/** \file type_mesh.c
  * \brief Subroutines to create and work with meshes which allows negative indices and 1D/2D/3D access (governed by predefined constants
  * mc_have_x, mc_have_y, mc_have_z and macroses to avoid any memory/performance penalty). See type_mesh.h for macroses and type definitions.
  */

#include <unistd.h>
#include <stdlib.h>
#include <assert.h>

#include <mpi.h>

#include "type_mesh.h"

#include "log.h"
#include "misc_cfgReader.h"

#define mf_min(a, b) (((a) < (b)) ? (a) : (b))	///< \b Min of two numbers (typeless).
#define mf_max(a, b) (((a) > (b)) ? (a) : (b))	///< \b Max of two numbers (typeless).

#ifndef SGI
  // My Linux doesn't support files bigger 2 Gb.
  #define ftell64 ftell			///< Alias for non-SGI systems.
  #define fseek64 fseek			///< Alias for non-SGI systems.
#endif

// ---------------------------------------------------------------------------
/// If external define key MC_FILL_MESH_WITH_NAN is activated than fills mesh with qNaNs to
/// catch uninitialized nodes.
// ---------------------------------------------------------------------------
static void
mesh_NaNify (mesh_p mesh)
{
  assert (mesh && mesh->storage && mesh->size);
#ifdef MC_FILL_MESH_WITH_NAN												// Fills mesh with qNaN if asked.
  memset (mesh->storage, -1, mesh->size);
  SAY_DEBUG ("mesh_NaNify: fills mesh '%s' with NaNs.", (mesh->name) ? mesh->name : "?");
#endif
}

/* ======================================================================================= */
/*                    Memory (re)allocation and initialization part                        */
/* ======================================================================================= */

// ---------------------------------------------------------------------------
/// Sets sizes of the mesh-container.
// ---------------------------------------------------------------------------
static void
mesh_setLimits (mesh_p mesh, int imin, int jmin, int kmin, int imax, int jmax, int kmax)
{
  int sizeX, sizeY, sizeZ;

  imin *= mc_have_x;
  imax *= mc_have_x;
  jmin *= mc_have_y;
  jmax *= mc_have_y;
  kmin *= mc_have_z;
  kmax *= mc_have_z;

  mf_checkTypeID(mesh->type, "mesh_setLimits");

  ENSURE (imin <= imax && jmin <= jmax && kmin <= kmax, "bad arguments");

  mesh->imin = imin;
  mesh->jmin = jmin;
  mesh->kmin = kmin;
  mesh->imax = imax;
  mesh->jmax = jmax;
  mesh->kmax = kmax;

  /* Assignes axis limits and applies boundary degeneracy rules. */
  sizeX = imax - imin + 1;
  sizeY = jmax - jmin + 1;
  sizeZ = kmax - kmin + 1;

  /* Updates memory management parameters. */
  mesh->width_z = sizeZ;
  mesh->width_yz = sizeY*sizeZ;
  mesh->size = sizeX*sizeY*sizeZ*mf_mesh_sizeofNode (mesh);
}

// ---------------------------------------------------------------------------
/// Puts pointer on the origin of coordinate system inside the storage.
// ---------------------------------------------------------------------------
static void
mesh_setOrigin (mesh_p mesh)
{
  assert (mesh->type == mc_double || mesh->type == mc_vec3D_t);

  mesh->origin = mesh->storage;
  mesh->origin = (void *) mf_mesh_bytePointer(mesh, -mesh->imin, -mesh->jmin, -mesh->kmin);
}

// ---------------------------------------------------------------------------
/// \brief Allocates memory and shifts reference point to account for given axis limits.
/// \b Note: that is the only place you can set the type of the content (scalar/mc_double or vector/mc_vec3D_t).
// ---------------------------------------------------------------------------
void
mesh_allocate (mesh_p mesh, int imin, int jmin, int kmin, int imax, int jmax, int kmax, const char *name, int type)
{
  mf_checkTypeID(type, "mesh_allocate");										// Checks type of the mesh.

  strncpy (mesh->name, name, mc_meshNameSize - 1);									// Makes copy of the mesh name.
  mesh->name[mc_meshNameSize-1] = 0;

  mesh->type = type;													// Sets type of the mesh.
  mesh_setLimits (mesh, imin, jmin, kmin, imax, jmax, kmax);								// Sets imin, .., kmax, size, width_(y)z fields.

  mesh->storage = malloc (mesh->size);
  ENSURE (mesh->storage, "cannot malloc %.3f Kb", mesh->size/1024.0);

  mesh_NaNify (mesh);													// Optionally fills storage with NaNs.

  mesh_setOrigin (mesh);											// Envelops shift (-imin, -jmin, -kmin) into the origin pointer.
}

// ---------------------------------------------------------------------------
/// Frees allocated memory. Proper usage of this function avoids memory leaks :-))).
// ---------------------------------------------------------------------------
void
mesh_free (mesh_p mesh)
{
  if (mesh->storage)
    free (mesh->storage);
  mesh->storage = mesh->origin = NULL;
  mesh->type = -1;
}

// ---------------------------------------------------------------------------
/// \brief That is a func to call in case of reconfiguration of domain.
/// <b>Unlike realloc (..) content of the array is not saved but mixed up and destroyed</b>.
/// I do not use "realloc" because of mesh is warped and old coordinates do not work.
// ---------------------------------------------------------------------------
void
mesh_resize (mesh_p mesh, int imin, int jmin, int kmin, int imax, int jmax, int kmax)
{
  if (mesh->storage) 													// Frees memory.
    free (mesh->storage);

  mesh_setLimits (mesh, imin, jmin, kmin, imax, jmax, kmax);								// Calculates all sizes.
  mesh->storage = malloc (mesh->size);											// Allocates new chunk of memory.
  ENSURE (mesh->storage, "cannot malloc %.3f Kb", mesh->size/1024.0);

  mesh_NaNify (mesh);													// Optionally fills storage with NaNs.
  mesh_setOrigin (mesh);												// Tunes origin pointer.
}



/* ======================================================================================= */
/*              Parallel program management: mesh resizing, repartition, etc.              */
/* ======================================================================================= */

// ---------------------------------------------------------------------------
/// Packs mesh into 1D array.
// ---------------------------------------------------------------------------
double *
mesh_pack (mesh_p mesh, int i1, int j1, int k1, int i2, int j2, int k2, int *Size)
{
  int i, j, size, dJ, dK;
  double *data;

  mf_checkTypeID(mesh->type, "mesh_pack");							// All basic validity checks.

  if (i1 > i2 || j1 > j2 || k1 > k2)
  {
    *Size = 0;
    return NULL;
  }

  if (mf_mesh_pointIsOutside (mesh, i1, j1, k1) || mf_mesh_pointIsOutside (mesh, i2, j2, k2))
  {
    say ("mesh_pack: error detected.");
    say ("mesh_pack: region to pack (%d, %d, %d) - (%d, %d, %d).", i1, j1, k1, i2, j2, k2);
    mf_mesh_dumpInfo (mesh);
    DIE ("domain to pack is bigger than mesh");
  }

  dJ = j2 - j1 + 1;
  dK = k2 - k1 + 1;

  /* If overlapped region is not empty I save data. */
  size = mf_mesh_sizeofNode (mesh)/sizeof (double);
  ENSURE (mf_mesh_sizeofNode (mesh) % sizeof (double) == 0,
          "cannot pack node in doubles");

  *Size = (i2 - i1 + 1)*dJ*dK*size;

  data = (double *) malloc (*Size*sizeof (double));
  ENSURE (data, "cannot allocate memory for packed data");

  for (i = i1 ; i <= i2 ; i++)
    for (j = j1 ; j <= j2 ; j++)
      memcpy (data + size*dK*(j - j1 + dJ*(i - i1)), ((double *)(mesh->origin)) + size*(i*mc_have_x*mesh->width_yz + j*mc_have_y*mesh->width_z + k1*mc_have_z),
              size*dK*sizeof (double));

  return data;
}

// ---------------------------------------------------------------------------
/// Unpacks mesh from 1D array filled by mesh_pack.
// ---------------------------------------------------------------------------
void
mesh_unpack (mesh_p mesh, int i1, int j1, int k1, int i2, int j2, int k2, double *data, int Size)
{
  int i, j, size, dJ, dK;

  if (!data)														// Returns if buffer is empty.
    return;

  // Checks that type exists and all mesh->type is valid.
  mf_checkTypeID(mesh->type, "mesh_unpack");

  if (mf_mesh_pointIsOutside (mesh, i1, j1, k1) || mf_mesh_pointIsOutside (mesh, i2, j2, k2))
    DIE ("domain to pack is not covered by mesh");

  ENSURE (i1 < i2 && j1 < j2 && k1 < k2, "bad domain");

  size = mf_mesh_sizeofNode (mesh)/sizeof (double);
  ENSURE (mf_mesh_sizeofNode (mesh) % sizeof (double) == 0,
          "cannot unpack node as doubles");

  dJ = j2 - j1 + 1;
  dK = k2 - k1 + 1;
  ENSURE (Size == (i2 - i1 + 1)*dJ*dK*size, "packed data does not fit");

  for (i = i1 ; i <= i2 ; i++)
    for (j = j1 ; j <= j2 ; j++)
      memcpy (((double *)(mesh->origin)) + size*(i*mc_have_x*mesh->width_yz + j*mc_have_y*mesh->width_z + k1*mc_have_z), data + size*dK*(j - j1 + dJ*(i - i1)),
              size*dK*sizeof (double));
}

// ---------------------------------------------------------------------------
/// Changes size of the mesh and keeps data in an overlapped region.
// ---------------------------------------------------------------------------
void
mesh_reallocate (mesh_p mesh, int new_imin, int new_jmin, int new_kmin, int new_imax, int new_jmax, int new_kmax)
{
  int data_imin, data_jmin, data_kmin, data_imax, data_jmax, data_kmax, size, dI, dJ, dK;
  double *data = NULL;

  // Checks that type exists and all mesh->type is valid.
  mf_checkTypeID(mesh->type, "mesh_reallocate");

  new_imin *= mc_have_x;
  new_jmin *= mc_have_y;
  new_kmin *= mc_have_z;
  new_imax *= mc_have_x;
  new_jmax *= mc_have_y;
  new_kmax *= mc_have_z;

  /* Define size of the overlapped region which I need to save and copy into the new location. */
  data_imin = mf_max (mesh->imin, new_imin);
  data_jmin = mf_max (mesh->jmin, new_jmin);
  data_kmin = mf_max (mesh->kmin, new_kmin);
  data_imax = mf_min (mesh->imax, new_imax);
  data_jmax = mf_min (mesh->jmax, new_jmax);
  data_kmax = mf_min (mesh->kmax, new_kmax);

  dI = data_imax - data_imin + 1;
  dJ = data_jmax - data_jmin + 1;
  dK = data_kmax - data_kmin + 1;

  /* If overlapped region is not empty I save data. */
  if (dI > 0 && dJ > 0 && dK > 0)
    data = mesh_pack (mesh, data_imin, data_jmin, data_kmin, data_imax, data_jmax, data_kmax, &size);

  /* Resizes shape of mesh. */
  mesh_resize (mesh, new_imin, new_jmin, new_kmin, new_imax, new_jmax, new_kmax);
  mf_mesh_clean (mesh);

  if (dI > 0 && dJ > 0 && dK > 0)
    mesh_unpack (mesh, data_imin, data_jmin, data_kmin, data_imax, data_jmax, data_kmax, data, size);
}



/* ======================================================================================= */
/*                                         Disk IO.                                        */
/* ======================================================================================= */

// ---------------------------------------------------------------------------
/// \brief Saves entire mesh (geometry and content). If changed you must also correct the IO_loadSplittedMesh function and tecIO routines.
// ---------------------------------------------------------------------------
void
mesh_save (mesh_RO_p mesh, const char *name, int imin, int jmin, int kmin, int imax, int jmax, int kmax)
{
  if (mf_mesh_pointIsOutside (mesh, imin, jmin, kmin) || mf_mesh_pointIsOutside (mesh, imax, jmax, kmax))
    DIE ("region to save is bigger than mesh");

  imin *= mc_have_x;													// Removes duplicates due to degnerated axises.
  jmin *= mc_have_y;
  kmin *= mc_have_z;
  imax *= mc_have_x;
  jmax *= mc_have_y;
  kmax *= mc_have_z;

  FILE *fp = cfg_open (name, "wb", __func__);										// Opens file.
  fwrite (&imin, sizeof (int), 1, fp);											// Saves valid region.
  fwrite (&jmin, sizeof (int), 1, fp);
  fwrite (&kmin, sizeof (int), 1, fp);
  fwrite (&imax, sizeof (int), 1, fp);
  fwrite (&jmax, sizeof (int), 1, fp);
  fwrite (&kmax, sizeof (int), 1, fp);
  fwrite (mesh, sizeof (mesh_t), 1, fp);										// Saves structure fields.

  int span = kmax - kmin + 1;												// Defines size of the continuous chunk.
  for (int i = imin ; i <= imax ; i++)											// Saves body of the mesh.
    for (int j = jmin ; j <= jmax ; j++)
      fwrite (mf_mesh_bytePointer(mesh, i, j, kmin), mf_mesh_sizeofNode (mesh), span, fp);
  fclose (fp);														// Closes file.
}

// ---------------------------------------------------------------------------
/// \brief Loads entire mesh, geometry is reshaped to reproduce saved one on request, original name remains.
/// You can choose to extend mesh (to include old ghost regions) or to fit the mesh to the piece of
/// data shape.
// ---------------------------------------------------------------------------
void
mesh_load (mesh_p mesh, FILE *fp, int shrinkToValidPart)
{
  int i, j, span, type, imin, jmin, kmin, imax, jmax, kmax;
  void *storage;
  void *origin;
  char name[mc_meshNameSize];

  storage = mesh->storage;												// Saves buffers, pointers, name and type.
  origin = mesh->origin;
  type = mesh->type;													// Saves type.
  strncpy (name, mesh->name, mc_meshNameSize);										// Saves name.

  fread (&imin, sizeof (int), 1, fp);											// Loads valid region limits.
  fread (&jmin, sizeof (int), 1, fp);
  fread (&kmin, sizeof (int), 1, fp);
  fread (&imax, sizeof (int), 1, fp);
  fread (&jmax, sizeof (int), 1, fp);
  fread (&kmax, sizeof (int), 1, fp);
  fread (mesh, sizeof (mesh_t), 1, fp);											// Loads structure/names/etc.
  mesh->storage = storage;												// Restores pointers to catch allocated memory.
  mesh->origin  = origin;
  strncpy (mesh->name, name, mc_meshNameSize);										// Restores name.

  // Checks type and reshapes storage.
  ENSURE (mesh->type == type,
          "inconsistent type of mesh (loading of '%s' mesh, memory type %d, disk type %d)", mesh->name, type, mesh->type);
  if (shrinkToValidPart)
    mesh_resize (mesh, imin, jmin, kmin, imax, jmax, kmax);
  else
    mesh_resize (mesh, mesh->imin, mesh->jmin, mesh->kmin, mesh->imax, mesh->jmax, mesh->kmax);

  if (mf_mesh_pointIsOutside (mesh, imin, jmin, kmin) || mf_mesh_pointIsOutside (mesh, imax, jmax, kmax))
    DIE ("domain to load is bigger than mesh");

  span = kmax - kmin + 1;												// Defines size of the continuous chunk.
  for (i = imin ; i <= imax ; i++)											// Loads valid part of the mesh.
    for (j = jmin ; j <= jmax ; j++)
      fread (mf_mesh_bytePointer(mesh, i, j, kmin), mf_mesh_sizeofNode (mesh), span, fp);
}

// ---------------------------------------------------------------------------
/// \brief Uploads part of the saved mesh. Used for parallel repartitioned loads.
// ---------------------------------------------------------------------------
void
mesh_upload (mesh_p mesh, const reg_t *REG, FILE *fp)
{
  int i, j, imin, jmin, kmin, imax, jmax, kmax;
  reg_t reg = *REG;
  mesh_t src;

  fread (&imin, sizeof (int), 1, fp);											// Loads valid region limits.
  fread (&jmin, sizeof (int), 1, fp);
  fread (&kmin, sizeof (int), 1, fp);
  fread (&imax, sizeof (int), 1, fp);
  fread (&jmax, sizeof (int), 1, fp);
  fread (&kmax, sizeof (int), 1, fp);
  fread (&src, sizeof (mesh_t), 1, fp);											// Loads structure/names/etc.

  // Checks type.
  ENSURE (mesh->type == src.type,
          "inconsistent type of mesh (loading of '%s' mesh, memory type %d, disk type %d)", mesh->name, src.type, mesh->type);

  // Checks avaliable space.
  if (((mesh->imin > reg.min[0] || imin > reg.min[0] || mesh->imax < reg.max[0] || imax < reg.max[0]) && mc_have_x) ||
      ((mesh->jmin > reg.min[1] || jmin > reg.min[1] || mesh->jmax < reg.max[1] || jmax < reg.max[1]) && mc_have_y) ||
      ((mesh->kmin > reg.min[2] || kmin > reg.min[2] || mesh->kmax < reg.max[2] || kmax < reg.max[2]) && mc_have_z))
  {
    say ("%s: too big region to upload, dumping args ...", __func__);
    say ("  mesh: (%d, %d, %d) - (%d, %d, %d),", mesh->imin, mesh->jmin, mesh->kmin, mesh->imax, mesh->jmax, mesh->kmax);
    say ("  saved: (%d, %d, %d) - (%d, %d, %d),", imin, jmin, kmin, imax, jmax, kmax);
    say ("  region: (%d, %d, %d) - (%d, %d, %d),", reg.min[0], reg.min[1], reg.min[2], reg.max[0], reg.max[1], reg.max[2]);
    DIE ("domain to load is bigger than mesh");
  }

  // Shrinks size along degenerated axises to avoid repeated loading of the same piece.
  mf_reg_collapse (&reg);

  long long int posShift = ftell64 (fp);
  for (i = reg.min[0] ; i <= reg.max[0] ; i++)										// Loads region.
    for (j = reg.min[1] ; j <= reg.max[1] ; j++)
    {
      const int k = reg.min[2];
      fseek64 (fp, posShift + mf_mesh_sizeofNode (mesh)*(k - kmin + (kmax - kmin + 1)*(j - jmin + (jmax - jmin + 1)*(i - imin))), SEEK_SET);
      fread (mf_mesh_bytePointer(mesh, i, j, k), mf_mesh_sizeofNode (mesh), reg.max[2] - k + 1, fp);
    }
}
