/** \file spectr_dump.c
  * Saves raw data to the fragmented files for futher processing (see spectr.out diagnostic and spectr_process.h).
  */

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "type_mesh2.h"
#include "type_CFile.h"

#include "log.h"
#include "misc_cfgReader.h"

static int spectrNum = -1;			///< Number of the next record to export.

// ---------------------------------------------------------------------------
/// Resets number of the next record to 0.
// ---------------------------------------------------------------------------
void
spectr_startDump (void)
{
  assert (spectrNum < 0);												// Assumes module was unused.
  spectrNum = 0;
}

// ---------------------------------------------------------------------------
/// Sets number of the next record using number of records stored in the file <b>'binData/spectr.N'</b>.
// ---------------------------------------------------------------------------
void
spectr_continueDump (void)
{
  spectrNum = 0;
  FILE *fp = fopen ("binData/spectr.N", "rt");
  if (fp)
  {
    spectrNum = cfg_readInt (fp);
    fclose (fp);
  }
}

// ---------------------------------------------------------------------------
/// Saves vector mesh components into 3 different files for processing by \b spectr.out.
// ---------------------------------------------------------------------------
static void
spectr_dumpVec (CFile_t *file, const char *prefix, reg_t *dumpReg, reg_t *map, const vec3D_t *data)
{
  CF_openChunk (file, "%s.x", prefix);
  mf_scanRegion (dumpReg, i, j, k)
    CF_write (&mv_map(map, data, i, j, k).x, sizeof (double), 1, file);
  CF_closeChunk (file);

  CF_openChunk (file, "%s.y", prefix);
  mf_scanRegion (dumpReg, i, j, k)
    CF_write (&mv_map(map, data, i, j, k).y, sizeof (double), 1, file);
  CF_closeChunk (file);

  CF_openChunk (file, "%s.z", prefix);
  mf_scanRegion (dumpReg, i, j, k)
    CF_write (&mv_map(map, data, i, j, k).z, sizeof (double), 1, file);
  CF_closeChunk (file);
}

// ---------------------------------------------------------------------------
/// Saves meshes for spectr.out (\b warning: please pass only region [min, max) for each cpu, so there will be no duplication or boundary violation).
// ---------------------------------------------------------------------------
void
spectr_dump (double time, int cpu, int cpuN, reg_t *dumpReg, reg_t *map, const vec3D_t *dataE, const vec3D_t *dataH)
{
  mf_mesh_initUnroll (dumpReg);												// Checks and polishes mappers.
  mf_mesh_initUnroll (map);

  ENSURE (reg_isInside (dumpReg, map),
          "subregion to dump is outside of the array");

  CFile_t *file = CF_openWrite ("tmp/spectrDump_%d_%d", spectrNum, cpu);

  CF_openChunk (file, "global parameters");
  CF_print (file, "%d\t  cpus total.\n", cpu_total);
  CF_print (file, "%d\t  i-min.\n", dmn_min[0]*mc_have_x);
  CF_print (file, "%d\t  j-min.\n", dmn_min[1]*mc_have_y);
  CF_print (file, "%d\t  k-min.\n", dmn_min[2]*mc_have_z);
  CF_print (file, "%d\t  i-max.\n", dmn_max[0]*mc_have_x);
  CF_print (file, "%d\t  j-max.\n", dmn_max[1]*mc_have_y);
  CF_print (file, "%d\t  k-max.\n", dmn_max[2]*mc_have_z);
  CF_print (file, "%.14e\t  time.\n", time);
  CF_print (file, "%e\t  h1.\n", h1);
  CF_print (file, "%e\t  h2.\n", h2);
  CF_print (file, "%e\t  h3.\n", h3);
  CF_closeChunk (file);

  CF_openChunk (file, "cpu region");
  CF_write (dumpReg, sizeof (*dumpReg), 1, file);
  CF_closeChunk (file);

  spectr_dumpVec (file, "E", dumpReg, map, dataE);
  spectr_dumpVec (file, "H", dumpReg, map, dataH);
  CF_close (file);

  ++spectrNum;

  if (cpu_here)														// Master will update the counter file.
    return;

  FILE *fp = cfg_open ("binData/spectr.N", "wt", __func__);								// Updates total number of records.
  fprintf (fp, "@ %d  spectr dumps written.\n", spectrNum);
  fclose (fp);
}
