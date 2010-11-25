/** \file em_gaussSpot.h
  * Hard source of the radiation in form of X-plane with radiating gauss-spot.
  * (See setup/tag_gaussSpot.h).
  */

#ifndef em_gaussSpot_header
#define em_gaussSpot_header	///< \internal Guard.

#include "type_mesh.h"

void gaussSpot_init (void);
void gaussSpot_deactivateAll (void);
void gaussSpot_forceE (meshVec_p E, const double t);

#endif
