/** \file plasma.h
  * Plasma storage and IO.
  *
  * Plasma is stored as set of linear arrays wrapped into one-way connected
  * list. All objects are allocated using one of the plasma creating tags
  * and can be referenced at any time using reverse counter (last object has
  * index '0', one before last '1' and so on).
  */

#ifndef MC_SETUP_PLASMA_HEADER
#define MC_SETUP_PLASMA_HEADER

#include "type_mesh.h"
#include "type_marker.h"

void      plasma_save      (const char *name);
void      plasma_newObject (void);
marker_t* plasma_getObject (int count, long int *N);
void      plasma_rho       (meshDouble_p rho);
marker_t* plasma_marker    (void);

#endif
