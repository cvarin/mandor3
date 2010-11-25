/** \file IO_tecplot.h
  * IO of the charge density, current density and fields with interpolation
  * 'on the fly' for field visualization.
  */

#ifndef IO_TECPLOT_HEADER
#define IO_TECPLOT_HEADER

#include "type_mesh.h"
#include "misc_partition.h"

int  tecIO_recordsTotal (void);
void tecIO_parameters   (int recordNum, double *time, int *cpuN, int *fileMapID);

void tecIO_setRecordNum (int num);

void tecIO_saveFields   (double time, meshVec_RO_p E, meshVec_RO_p H);
void tecIO_saveCurrents (meshVec_RO_p J, meshDouble_RO_p rho);

void tecIO_load         (int recordNum,
                         double *time, const reg_t *reg,
                         meshVec_p E,
                         meshVec_p H,
                         meshVec_p J,
                         meshDouble_p rho);

#endif
