/*
 * Input/Output of the all meshes --- charge density, currents, fields. Used to save data for TecPlot.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "type_mesh.h"
#include "type_marker.h"

#include "IO_sys.h"
#include "IO_mesh.h"
#include "IO_names.h"
#include "IO_fileMap.h"

#include "log.h"
#include "misc_MPItags.h"
#include "misc_partition.h"
#include "misc_cfgReader.h"

int  sysNum = - 1;

// ---------------------------------------------------------------------------
/// Gets global parameters of the given check-point.
// ---------------------------------------------------------------------------
void
sysIO_parameters (int recordNum, double *time, int *cpuN, int *fileMapID)
{
  struct
  {
    int    cpuN, fMapID;
    double time;
  } pack;

  if (!cpu_here)
  {
    FILE *fp = cfg_open (IO_nameRec (sysName_info, recordNum), "rt", __func__);
    pack.time = cfg_readDouble (fp);
    pack.cpuN = cfg_readInt (fp);
    pack.fMapID = cfg_readInt (fp);
    fclose (fp);
  }

  MPI_Bcast (&pack, sizeof (pack), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (time)       *time = pack.time;						// Returns requested parameters.
  if (cpuN)       *cpuN = pack.cpuN;
  if (fileMapID)  *fileMapID = pack.fMapID;
}

// ---------------------------------------------------------------------------
/// Loads number of system check-points record written at the moment.
// ---------------------------------------------------------------------------
int
sysIO_recordsTotal (void)
{
  int num = -1;

  if (!cpu_here)
  {
    FILE *fp = fopen ("binData/sys.N", "rt");
    if (fp)
    {
      num = cfg_readInt (fp);
      fclose (fp);
    }
  }

  MPI_Bcast (&num, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return num;
}

// ---------------------------------------------------------------------------
/// Sets record number to continue from.
// ---------------------------------------------------------------------------
void
sysIO_setRecordNum (int num)
{
  if (num >= 0)														// Direct-set branch.
  {
    sysNum = num;
    return;
  }

  sysNum = sysIO_recordsTotal ();											// Continue-from-possible branch.
  if (sysNum < 0)
    sysNum = 0;
}

// ---------------------------------------------------------------------------
/// Saves all system data (check-point) including particles.
// ---------------------------------------------------------------------------
void
sysIO_save (double time, meshVec_RO_p E, meshVec_RO_p H)
{
  parameter_save ();
  const int fileMap = fileMap_save ();

  // Saves meshes.
  mesh_save (mcast_mesh_RO (E), IO_nameCpuRec (sysName_E_full, cpu_here, sysNum), cpu_min[0], cpu_min[1], cpu_min[2], cpu_max[0], cpu_max[1], cpu_max[2]);
  mesh_save (mcast_mesh_RO (H), IO_nameCpuRec (sysName_H_full, cpu_here, sysNum), cpu_min[0], cpu_min[1], cpu_min[2], cpu_max[0], cpu_max[1], cpu_max[2]);

  partition_save (IO_nameRec (sysName_prtn_full, sysNum));			// Saves partition.
  ++sysNum;														// Switches to the next record.

  if (cpu_here)
    return;

  FILE *fp = cfg_open (IO_nameRec (sysName_info, sysNum - 1), "wt", __func__);						// Saves header and descriptor.
  fprintf (fp, "@ %.14e    time.\n", time);
  fprintf (fp, "@ %d       cpus total.\n", cpu_total);
  fprintf (fp, "@ %d       file map id.\n", fileMap);
  fprintf (fp, "That is checkpoint (record of the full system state) #%d.\n", sysNum - 1);
  fclose (fp);

  fp = cfg_open ("binData/sys.N", "wt", __func__);									// Updates total number of records.
  fprintf (fp, "@ %d\t  system check-point records written.\n", sysNum);
  fclose (fp);
}

// ---------------------------------------------------------------------------
/// Loads data from system checkpoint (set E or H to NULL to skip it).
// ---------------------------------------------------------------------------
void
sysIO_loadEM (int recordNum, double *time, const reg_t *reg, meshVec_p E, meshVec_p H)
{
  int fileMap;
  fileMap_t fMap = mc_fileMap_init;

  sysIO_parameters (recordNum, time, NULL, &fileMap);
  fileMap_load (&fMap, fileMap);											// Loads file map.

  if (E)	IO_loadMesh (IO_nameRec (sysName_E_map, recordNum), mcast_mesh (E), reg, &fMap);			// Loads only requested meshes.
  if (H)	IO_loadMesh (IO_nameRec (sysName_H_map, recordNum), mcast_mesh (H), reg, &fMap);

  fileMap_free (&fMap);
}
