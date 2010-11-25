/** \file IO_mesh.c
  * Input/output of the given subregion of the distributed mesh.
  *
  * \warning For now support of the distributed file system is removed.
  *
  * Algorithm: loads file map, partitions region on this map and gets list of the regions (each entry in
  * list depends only on one file). Than all pieces are loaded.
  *
  * \note For distributed file systems it is better to use socket layer for non-blocking exchanges.
  *
  * \todo Simultaneous access to the same file causes hard-drive penalty => may be it is necessary to
  *       collect line of requests and provide access as single file/single reader. How to check: setup
  *       with one cpu and big domain and do loading with many-many cpus.
  */

#include "log.h"
#include "IO_mesh.h"
#include "misc_cfgReader.h"

// ---------------------------------------------------------------------------
/// Simple wrapper for \b sprintf to provide convinient way for nodeNum->name conversion.
// ---------------------------------------------------------------------------
static char *
name (const char *format, int nodeNum)
{
  static char name[100];
  ENSURE (nodeNum >= 0 && nodeNum < 1000,
          "bad parameter for node number format; len('%d') > 3", nodeNum);
  sprintf (name, format, nodeNum);
  return name;
}

// ---------------------------------------------------------------------------
/// Loads mesh distributed over many nodes and common disk space.
// ---------------------------------------------------------------------------
void
IO_loadMesh (const char *format, mesh_p mesh, const reg_t *reg, fileMap_t *fMap)
{
  regList_t pieces = mc_regList_init;											// List of pieces partitioned by map.
  reg_distributeOnMap (reg, &fMap->map, &pieces);

  for (int claim = 0 ; claim < pieces.N ; ++claim)									// Loads all readable sub-regions.
  {
    FILE *fp = cfg_open (name (format, pieces.list[claim].cpu), "rb", __func__);
    mesh_upload (mesh, pieces.list + claim, fp);
    fclose (fp);
  }
  regList_clean (&pieces);
}
