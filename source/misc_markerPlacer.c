/** \file misc_markerPlacer.c
  * \brief This module sorts particles and redistirbutes them to the owners' nodes. Uses fixed
  * (limited) size of the buffers to reduce memory pressure for bug 3D simulations. Call to this
  * functions guarantees that all particles are located on the proper nodes.
  *
  * <h3>Motivation.</h3>
  * This function is rather initialization/reconfiguration type of tools rather than high
  * performance routine. It is used after setup or loading to make sure that main step will start
  * work smoothly.
  *
  * <h3>Design ideas:</h3>
  * Parallel exchange is done by all2all module (see info in misc_all2all.h). What is done here is
  * sorting of the particles on local/scheduled and removing of the local particles from the
  * pipeline ASAP. Periodic boundary conditions are applied for local particles and particles
  * behind the absorbing boundary conditions are removed.
  */

#include "misc_PIC.h"
#include "log.h"
#include "misc_cfgReader.h"
#include "misc_markerPacker.h"

// ---------------------------------------------------------------------------
/// Macro-constructor of the proper switch statements for boundary condition processing.
// ---------------------------------------------------------------------------
#define MF_PLACER_BC(AXIS)									\
{												\
  if (f->AXIS < min[mc_ ## AXIS] && mc_have_ ## AXIS)						\
    switch (dmn.bound_## AXIS ## Min)								\
    {												\
      case BC_PERIODIC:										\
        f->AXIS += wrap[mc_ ## AXIS];								\
        break;											\
        											\
      case BC_MIRROR:										\
      case BC_OPEN:										\
        f->AXIS = 2*min[mc_ ## AXIS] - f->AXIS;							\
        break;											\
    }												\
    												\
  if (f->AXIS > max[mc_ ## AXIS] && mc_have_ ## AXIS)						\
    switch (dmn.bound_## AXIS ## Max)								\
    {												\
      case BC_PERIODIC:										\
        f->AXIS -= wrap[mc_ ## AXIS];								\
        break;											\
        											\
      case BC_MIRROR:										\
      case BC_OPEN:										\
        f->AXIS = 2*max[mc_ ## AXIS] - f->AXIS;							\
        break;											\
    }												\
}

// ---------------------------------------------------------------------------
/// Takes all particles and sends them to the owners. Local portal is slightly extended
/// to grab particles we can safely push ourself and remove boundary uncertainty.
// ---------------------------------------------------------------------------
void
placer_exchange (void)
{
  ENSURE (0, "Particles sorter is not reimplemented.");
  say ("placer_exchange: particles are forced to their domains.");
}
