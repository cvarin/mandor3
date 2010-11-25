/** \file tag_boundary.c
  * Sets boundary conditions.
  *
  * Example of the config file entry is
    <pre>
    [boundaryConditions]
    @ 0		xMin
    @ 0		xMax
    @ 0		yMin
    @ 0		yMax
    @ 0		zMin
    @ 0		zMax

      Type of boundary condition is:
        0 - periodic
        1 - mirror
        2 - Mur's first order absorbing boundary condition
    </pre>
  * \warning For 1D/2D simulations boundary conditions on the deactivated axises must be set to \b periodic.
  */

#include "misc_cfgReader.h"
#include "misc_parameters.h"

// ---------------------------------------------------------------------------
/// Reads boundary conditions and sets them.
// ---------------------------------------------------------------------------
void
tag_boundary (FILE *fp)
{
  int xMin, xMax, yMin, yMax, zMin, zMax;

  xMin = cfg_readInt (fp);
  xMax = cfg_readInt (fp);
  yMin = cfg_readInt (fp);
  yMax = cfg_readInt (fp);
  zMin = cfg_readInt (fp);
  zMax = cfg_readInt (fp);

  parameter_setupBounds (xMin, xMax, yMin, yMax, zMin, zMax);
}
