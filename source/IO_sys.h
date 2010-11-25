#ifndef MC_IO_SYS_HEADER
#define MC_IO_SYS_HEADER

#include "type_mesh.h"
#include "misc_partition.h"

// XXX Simple saving (no subsystem required?)
int  sysNum;

int  sysIO_recordsTotal (void);
void sysIO_parameters (int recordNum, double *time, int *cpuN, int *fileMapID);

void sysIO_setRecordNum (int num);

void sysIO_save (double time, meshVec_RO_p E, meshVec_RO_p H);
void sysIO_loadEM (int recordNum, double *time, const reg_t *reg, meshVec_p E, meshVec_p H);

#endif
