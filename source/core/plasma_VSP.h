/** \file plasma_VSP.h
  *
  * This set of routines does VSP transformation of the current and charge density as described by
  * <i>Vshivkov V.A., Romanov D.V., Snytnikov V.N., </i> <b>"The problem of spontaneous heating of model plasma in the method of particles"</b>,
  * Journal of Computational Technologies, Volume 4, N 3, 1999 year.
  *
  */

#ifndef MC_PLASMA_VSP_HEADER
#define MC_PLASMA_VSP_HEADER			///< \internal Guard.

#include "type_mesh.h"

void VSP_configure (double Alpha);

void VSP_vectorSweep (meshVec_t *J);
void VSP_scalarSweep (meshDouble_t *rho);

#endif
