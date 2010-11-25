/** \file tag_vShift.c
  * \brief Shifts last chapter or entire plasma in velocity space. Goal is to create beams with arbitrary shaped DF in independent passes.
  *
  * Example of the config file entry is
    <pre>
    [velocity shift]
    @ 0.01          Delta Vx
    @ 0             Delta Vy
    @ 0             Delta Vz
    @ all           Component to shift (can be 'all' or 'last');

      Introduces perturbation of the mean velocity of the chosen component(s).
    </pre>
  */

#include <stdlib.h>
#include <string.h>

#include "log.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "setup/plasma.h"

// ---------------------------------------------------------------------------
/// Tags body.
// ---------------------------------------------------------------------------
void
tag_velocityShift (FILE *fp)
{
  double vx = cfg_readDouble (fp)*mc_have_x;
  double vy = cfg_readDouble (fp)*mc_have_y;
  double vz = cfg_readDouble (fp)*mc_have_z;

  const char *shift = cfg_readWord (fp);					// Gets type of shift.
  int globalShift = cfg_identifyWord (shift, "all", 1, "last", 0, mc_cfgTermGuesses);
  ENSURE (globalShift != -1, "Unknown parameter '%s'.", shift);

  ENSURE (globalShift == 0, "Global shift is not implemented (no need?).");

  long int N;
  marker_t *p = plasma_getObject (0, &N), *end = p + N;
  for ( ; p < end ; ++p)
  {
    p->vx += vx;
    p->vy += vy;
    p->vz += vz;
  }
}
