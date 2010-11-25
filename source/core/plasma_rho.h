/** \file plasma_rho.h
  * Charge density evaluation without storage detalization (linear arrays in,
  * charge density out).
  *
  * That is not a core functionality, because of a Maxwell solver uses current
  * density rather than charge density. Thus this module is neither speed nor
  * parallel exchange optimized. During computing we drop a chunks of plasma in
  * to evaluate charge density and than do parallel exchange to complete charge
  * density evaluation. Boundary conditions and VSP sweep are applied inside of
  * the parallel exchange routines as well. Make sure you have cleared 'rho'
  * mesh before starting to throw markers in.
  */

#ifndef MC_PLASMA_RHO_HEADER
#define MC_PLASMA_RHO_HEADER

#include "type_mesh.h"
#include "type_marker.h"

void plasmaRho_startExchange (meshDouble_p rho);
void plasmaRho_finishExchange (meshDouble_p rho);
void plasmaRho_add (meshDouble_p rho, marker_t *plasma, long int N);
void plasmaRho_setupConnections (void);

#endif
