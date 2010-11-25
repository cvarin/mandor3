/** \file em.h
  * Maxwell solver.
  *
  * Second order Yee scheme with central differences and staggered meshes,
  * conserves \f$ div \vec E \f$ and \f$ div \vec H \f$ with numerical
  * accuracy.
  *
  * Staggered mesh indexing: each component may have only <b>negative shift in
  * half of spatial step</b>, if there are any. Index to real world position
  * mapping:
  * \f[
  *   \begin{array}{llll}
  *   \vec H[i][j][k] = & \{ H_{x\ i,\     j-1/2,\ k-1/2}^{n-1/2};\quad
  *                     &    H_{y\ i-1/2,\ j,\     k-1/2}^{n-1/2};\quad
  *                     &    H_{z\ i-1/2,\ j-1/2,\ k    }^{n-1/2} \}\\[3mm]
  *   \vec E[i][j][k] = & \{ E_{x\ i-1/2,\ j,\     k    }^{n};
  *                     &    E_{y\ i,\     j-1/2,\ k    }^{n};
  *                     &    E_{z\ i,\     j,    \ k-1/2}^{n} \}
  *   \end{array}
  * \f]
  *
  * In other words, 'mv_fx(E, i + 1, j, k)' is 'Ex((i + 0.5)*h₁, j*h₂, k*h₃)'.
  *
  * Some implementation details are documented in em.c.
  */

#ifndef MC_EM_HEADER
#define MC_EM_HEADER

#include "type_mesh.h"

void em_init         (meshVec_p    E, meshVec_p H);
void em_HHalfStep    (meshVec_RO_p E, meshVec_p H);
void em_HStep        (meshVec_RO_p E, meshVec_p H);

void em_EStep_start  (meshVec_p E, meshVec_RO_p H);
void em_EStep_finish (meshVec_p E, meshVec_p    H, meshVec_RO_p J);

void em_energy       (meshVec_RO_p E, meshVec_RO_p H, double *W_E, double *W_M);

#endif
