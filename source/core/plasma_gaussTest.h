/** \file plasma_gaussTest.h
  * Few callbacks to test charge conservation and Gauss law.
  *
  * Computes charge denstity before and after timestep to compare dρ/dt against
  * div(j) as a first step (involves only plasma module).
  *
  * Than compares 'div (E)' against '4πρ' (plasma + EM module together).
  *
  * Eats performance.
  *
  * \attention To enable the test, one has to enable macro-variables
  *   ∘ MC_GAUSS_TEST	saves min/max errors for each time-step in 'XXX'
  *   ∘ MC_GAUSS_DEBUG	saves detailed balancing (each node) in HUGE file 'XXX'
  */

#ifndef MC_PLASMA_GAUSSTEST_HEADER
#define MC_PLASMA_GAUSSTEST_HEADER

#include "type_mesh.h"

void gauss_before (void);
void gauss_after (meshVec_RO_p J, meshVec_RO_p E);

#endif
