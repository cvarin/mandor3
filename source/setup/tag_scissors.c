/** \file tag_scissors.c
  * Activates build in filtering of the markers.
  */

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"

// ---------------------------------------------------------------------------
/// \b CacheFilter entry point.
// ---------------------------------------------------------------------------
void
tag_scissors (FILE *fp)
{
    const double micron = units (mc_micron);
/*
    ENSURE (0, "I will implements scissoring on demand");
    // Reads bounding box.
    double x1 = cfg_readDouble (fp)*micron,
           x2 = cfg_readDouble (fp)*micron,
           y1 = cfg_readDouble (fp)*micron,
           y2 = cfg_readDouble (fp)*micron,
           z1 = cfg_readDouble (fp)*micron,
           z2 = cfg_readDouble (fp)*micron;

    plasma_set_filter (x1, x2, y1, y2, z1, z2);

    say ("tag_scissors:");
    say ("    All following particles will be filtered.");
    say ("    Filter (in [r₀]): [%e, %e] x [%e, %e] x [%e, %e]",
              x1, x2, y1, y2, z1, z2);
    say ("    Filter (in [μm]): [%e, %e] x [%e, %e] x [%e, %e]",
              x1/micron, x2/micron, y1/micron, y2/micron, z1/micron, z2/micron);*/
}
