/** \file parr_ghostCells.h
  * Expensive explicit syncronization of the ghost-cells which overlap with other cpus.
  *
  * <b>Content of the local ghost cells must be updated by using em-module's ebc_<> group
  * of functions before calling this module.</b>
  *
  * Requires huge amount of memory, so make sure you call it before full framework for
  * main loop parallel exchange is allocated - this way one may reuse MPI_Alloc_mem memory.
  */

#ifndef PARR_GHOSTCELLS_HEADER
#define PARR_GHOSTCELLS_HEADER							///< Multiple include guard.

#include "type_mesh.h"

void ghostSync_sync (meshVec_p E, meshVec_p H);

#endif
