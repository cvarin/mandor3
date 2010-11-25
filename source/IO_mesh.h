/** \file IO_mesh.h
  * Input/output of the given subregion of the distributed mesh.
  */

#ifndef IO_mesh_header
#define IO_mesh_header				///< \internal Guard.

#include "type_mesh.h"
#include "IO_fileMap.h"

void IO_loadMesh (const char *format, mesh_p mesh, const reg_t *reg, fileMap_t *fMap);

#endif
