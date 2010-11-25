/** \file spectr_dump.h
  * Saves raw data to the fragmented files for futher processing using spectr.out diagnostic.
  */

#ifndef mc_spectr_dump_header
#define mc_spectr_dump_header						///< \internal Guard.

#include "type_mesh.h"

void spectr_startDump (void);
void spectr_continueDump (void);
void spectr_dump (double time, int cpu, int cpuN, reg_t *dumpReg, reg_t *map, const vec3D_t *dataE, const vec3D_t *dataH);

#endif
