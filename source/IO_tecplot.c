/** \file IO_tecplot.c
  * IO of the charge density, current density and fields with interpolation
  * 'on the fly' for field visualization.
  *
  * This set of routines is used to save interpolated into the same mesh nodes
  * all vector and scalar fields. Interpolation is done here to simplify all
  * visualization codes (they shouldn't guess about boundary conditions and
  * other numerical stuff). Data is saved here are just after the PIC step,
  * so all EM fields are defined on single time layer.
  *
  * Uses IO_unroll.c to pipeline everything through small buffer of
  * mc_tecCapacity size (basically few hundreds of Kb).
  *
  * \note Currents are saved separatedly because of a parallel decomposition.
  *       Current density data arrives only after the full update of the EM
  *       field (when E and H are again defined on separate time layers).
  *       So EM field and current density are saved for the same time layer
  *       t = n⋅τ but in two function calls.
  */

#include <unistd.h>

#include <mpi.h>

#include "IO_mesh.h"
#include "IO_names.h"
#include "IO_unroll.c"

#include "log.h"
#include "misc_cfgReader.h"

/// Size of the temporary buffer for tecplot interpolating/saving.
#define  mc_tecCapacity (128*1024*8)

// ---------------------------------------------------------------------------
/// Function for interpolation of the vector data and placing result to
/// given address (buffer). There will be one for E/J and one for H. Strides
/// along axises are precalc'd ans stored as global static variables, no need
/// to pass them through arguments.
// ---------------------------------------------------------------------------
typedef void (*interpFunc_t) (meshVec_RO_p mesh, long int pos, vec3D_t *resLocation);

static long int dx;		///< Stride along \b X for interpolators.
static long int dy;		///< Stride along \b Y for interpolators.
static long int dz;		///< Stride along \b Z for interpolators.

static int  tecplotNum = - 1;

///< Buffer to store interpolated chunk before saving.
static char tecBuffer[mc_tecCapacity];

static const char name_info[50]      = "binData/tec_%06d.txt";			///< Format string for global parameters of the record files' manifold.
static const char name_prtn_full[50] = "binData/tec_prtn_%06d.bin";		///< Format string for name of partition file (recNum).

static const char name_E_full[50]    = "binData/tec_E_%06d(%03d).bin";		///< Format string for name of electric field file (recNum, node).
static const char name_H_full[50]    = "binData/tec_H_%06d(%03d).bin";		///< Format string for name of magnetic field file (recNum, node).
static const char name_J_full[50]    = "binData/tec_J_%06d(%03d).bin";		///< Format string for name of current density field file (recNum, node).
static const char name_rho_full[50]  = "binData/tec_rho_%06d(%03d).bin";	///< Format string for name of charge density field file (recNum, node).

static const char name_E_map[50]     = "binData/tec_E_%06d(%%03d).bin";		///< Format string for generator of name of electric field file (recNum, node).
static const char name_H_map[50]     = "binData/tec_H_%06d(%%03d).bin";		///< Format string for generator of name of magnetic field file (recNum, node).
static const char name_J_map[50]     = "binData/tec_J_%06d(%%03d).bin";		///< Format string for generator of name of electric field file (recNum, node).
static const char name_rho_map[50]   = "binData/tec_rho_%06d(%%03d).bin";	///< Format string for generator of name of electric field file (recNum, node).

// ---------------------------------------------------------------------------
/// Returns number of tecplot meshes located on the disk.
// ---------------------------------------------------------------------------
int
tecIO_recordsTotal (void)
{
    int num = -1;

    FILE *fp = fopen ("binData/tecplot.N", "rt");
    if (fp) {
        num = cfg_readInt (fp);
        fclose (fp);
    }
    return num;
}

// ---------------------------------------------------------------------------
/// Sets tecplot number of the next record.
// ---------------------------------------------------------------------------
void
tecIO_setRecordNum (int num)
{
    // Direct-positioning way.
    if (num >= 0) {
        tecplotNum = num;
        return;
    }

    // Continue-from-the-last-possible way.
    tecplotNum = tecIO_recordsTotal ();
    if (tecplotNum < 0)
        tecplotNum = 0;
}

// ---------------------------------------------------------------------------
/// Gets global parameters of the given record set.
// ---------------------------------------------------------------------------
void
tecIO_parameters (int recordNum, double *time, int *cpuN, int *fileMapID)
{
    struct {
        int    cpuN, fMapID;
        double time;
    } pack;

    if (!cpu_here) {
        FILE *fp = cfg_open (IO_nameRec (name_info, recordNum), "rb", __func__);
        fread (&pack.time, sizeof (double), 1, fp);
        fread (&pack.cpuN, sizeof (int), 1, fp);
        fread (&pack.fMapID, sizeof (int), 1, fp);
        fclose (fp);
    }

    MPI_Bcast (&pack, sizeof (pack), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (time)       *time = pack.time;
    if (cpuN)       *cpuN = pack.cpuN;
    if (fileMapID)  *fileMapID = pack.fMapID;
}

// ---------------------------------------------------------------------------
/// Macros to address unrolled non-typed storage on top of span-module.
// ---------------------------------------------------------------------------
#define MV_STORG(mesh, pos) ((vec3D_t*) (((char*)mesh->storage) + (pos)))

// ---------------------------------------------------------------------------
/// Interpolator for E-type of mesh field (see 'em.h' for shifts).
// ---------------------------------------------------------------------------
static void
tecIO_interpE (meshVec_RO_p E, long int pos, vec3D_t *res)
{
    res->x = 0.5*(MV_STORG(E, pos)->x + MV_STORG(E, pos + dx)->x);
    res->y = 0.5*(MV_STORG(E, pos)->y + MV_STORG(E, pos + dy)->y);
    res->z = 0.5*(MV_STORG(E, pos)->z + MV_STORG(E, pos + dz)->z);
}

// ---------------------------------------------------------------------------
/// Interpolator for H-type of mesh field (see 'em.h' for shifts).
// ---------------------------------------------------------------------------
static void
tecIO_interpH (meshVec_RO_p H, long int pos, vec3D_t *res)
{
    res->x = 0.25*(MV_STORG(H, pos     )->x + MV_STORG(H, pos + dy     )->x +
                   MV_STORG(H, pos + dz)->x + MV_STORG(H, pos + dy + dz)->x);
    res->y = 0.25*(MV_STORG(H, pos     )->y + MV_STORG(H, pos + dx     )->y +
                   MV_STORG(H, pos + dz)->y + MV_STORG(H, pos + dx + dz)->y);
    res->z = 0.25*(MV_STORG(H, pos     )->z + MV_STORG(H, pos + dx     )->z +
                   MV_STORG(H, pos + dy)->z + MV_STORG(H, pos + dx + dy)->z);
}

#undef MV_STORG

// ---------------------------------------------------------------------------
/// Saves vector field with intermidiate interpolation/unrolling intro fixed
/// size buffer to save RAM (no tmp meshes engaged).
// ---------------------------------------------------------------------------
static void
tecIO_saveVec (const char *name, meshVec_RO_p mesh, const reg_t *reg,
               interpFunc_t combiner)
{
    const int size = sizeof (vec3D_t);
    spanSet_t spans;
    ENSURE (! mf_mesh_pointIsOutside (mesh, reg->min[0], reg->min[1], reg->min[2])
         && ! mf_mesh_pointIsOutside (mesh, reg->max[0], reg->max[1], reg->max[2]),
            "region to save is bigger than mesh");

    ENSURE (mc_tecCapacity > 30*size && size == mf_mesh_sizeofNode (mesh),
            "too small buffer or bad sizeof for tuned packer");

    // Saves valid region.
    FILE *fp = cfg_open (name, "wb", __func__);
    fwrite (&reg->min[0], sizeof (int), 1, fp);
    fwrite (&reg->min[1], sizeof (int), 1, fp);
    fwrite (&reg->min[2], sizeof (int), 1, fp);
    fwrite (&reg->max[0], sizeof (int), 1, fp);
    fwrite (&reg->max[1], sizeof (int), 1, fp);
    fwrite (&reg->max[2], sizeof (int), 1, fp);
    fwrite (mesh, sizeof (mesh_t), 1, fp);

    // Storage domain.
    reg_t domain = {{mesh->imin, mesh->jmin, mesh->kmin},
                    {mesh->imax, mesh->jmax, mesh->kmax}, 0};

    // Forms list of spans to process.
    span_init (&domain, reg, &spans, size);

    // Exports strides for combiner.
    dx = size*mc_have_x*mesh->width_yz;
    dy = size*mc_have_y*mesh->width_z;
    dz = size*mc_have_z;

    while (!span_allTouched (&spans)) {
        // Packs buffer.
        int pos = 0;
        do {
            long int start,
                     length;

            // Gets span ("start + length").
            length = span_iterate (&spans, mc_tecCapacity - pos, &start);

            // Interpolates into buffer.
            for (int l = 0 ; l < length ; l += size, pos += size)
                combiner (mesh, start + l, (vec3D_t*) (tecBuffer + pos));
        } while (mc_tecCapacity - pos - size > 0 && !span_allTouched (&spans));

        fwrite (tecBuffer, 1, pos, fp);
    }
    fclose (fp);
}

// ---------------------------------------------------------------------------
/// Saves meshes for tecplot.out.
// ---------------------------------------------------------------------------
void
tecIO_saveFields (double time, meshVec_RO_p E, meshVec_RO_p H)
{
    parameter_save ();
    int fileMap = fileMap_save ();

    // Saving region.
    reg_t reg = {{cpu_min[0], cpu_min[1], cpu_min[2]},
                 {cpu_max[0], cpu_max[1], cpu_max[2]}, 0};

    // Saves interpolated vector fields.
    tecIO_saveVec (IO_nameCpuRec (name_E_full, cpu_here, tecplotNum), E, &reg, tecIO_interpE);
    tecIO_saveVec (IO_nameCpuRec (name_H_full, cpu_here, tecplotNum), H, &reg, tecIO_interpH);

    if (cpu_here)
        return;

    // Saves header and descriptor.
    FILE *fp = cfg_open (IO_nameRec (name_info, tecplotNum), "wb", __func__);
    fwrite (&time, sizeof (double), 1, fp);
    fwrite (&cpu_total, sizeof (int), 1, fp);
    fwrite (&fileMap, sizeof (int), 1, fp);
    fclose (fp);

    // Saves partition.
    partition_save (IO_nameRec (name_prtn_full, tecplotNum));
}

// ---------------------------------------------------------------------------
/// Saves particle constributions for tecplot.out.
// ---------------------------------------------------------------------------
void
tecIO_saveCurrents (meshVec_RO_p J, meshDouble_RO_p rho)
{
    // Saving region.
    reg_t reg = {{cpu_min[0], cpu_min[1], cpu_min[2]},
                 {cpu_max[0], cpu_max[1], cpu_max[2]}, 0};

    // Saves interpolated vector field.
    tecIO_saveVec (IO_nameCpuRec (name_J_full, cpu_here, tecplotNum), J, &reg, tecIO_interpE);
    mesh_save (mcast_mesh_RO (rho), IO_nameCpuRec (name_rho_full, cpu_here, tecplotNum), cpu_min[0], cpu_min[1], cpu_min[2], cpu_max[0], cpu_max[1], cpu_max[2]);

    tecplotNum++;

    FILE *fp = cfg_open ("binData/tecplot.N", "wt", __func__);			// Updates total number of records.
    fprintf (fp, "@ %d  tecplots' records written.\n", tecplotNum);
    fclose (fp);
}

// ---------------------------------------------------------------------------
/// Loads meshes for tecplot.out.
// ---------------------------------------------------------------------------
void
tecIO_load (int recordNum, double *time, const reg_t *reg, meshVec_p E, meshVec_p H, meshVec_p J, meshDouble_p rho)
{
    int fMapID;
    fileMap_t fMap = mc_fileMap_init;

    // Loads parameters.
    tecIO_parameters (recordNum, time, NULL, &fMapID);
    fileMap_load (&fMap, fMapID);

    // Loads only requested meshes.
    if (E)	IO_loadMesh (IO_nameRec (name_E_map, recordNum),   mcast_mesh (E),   reg, &fMap);
    if (H)	IO_loadMesh (IO_nameRec (name_H_map, recordNum),   mcast_mesh (H),   reg, &fMap);
    if (J)	IO_loadMesh (IO_nameRec (name_J_map, recordNum),   mcast_mesh (J),   reg, &fMap);
    if (rho)	IO_loadMesh (IO_nameRec (name_rho_map, recordNum), mcast_mesh (rho), reg, &fMap);

    fileMap_free (&fMap);
}
