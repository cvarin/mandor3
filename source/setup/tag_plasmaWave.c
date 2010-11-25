/** \file tag_plasmaWave.c
  * \brief Adds perturbation to the speed of particles \f$ \displaystyle \vec v_j \to \vec v_j + A\cdot \frac qM \cdot \cos (\vec k\cdot \vec r_j), \f$
  * where \f$ \displaystyle \vec k = \left(\frac {2\pi m_x}{L_x}, \frac {2\pi m_y}{L_y}, \frac {2\pi m_z}{L_z}\right) \f$.
  *
  * Example of the config file entry is
    <pre>
    [plasmaWave]
    @ 1 	wave number along X.
    @ 0 	wave number along Y.
    @ 1 	wave number along Z.
    @ 0.1e-3	Amplitude of the perturbation [c].
    </pre>
  *
  */

#include <math.h>
#include <stdlib.h>

#include "type_marker.h"

#include "setup/main.h"
#include "setup/plasma.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

// ---------------------------------------------------------------------------
/// Adds cos-like perturbation to the speed of particles.
// ---------------------------------------------------------------------------
void
tag_plasmaWave (FILE *fp)
{
  const int mx = cfg_readInt (fp)*mc_have_x;
  const int my = cfg_readInt (fp)*mc_have_y;
  const int mz = cfg_readInt (fp)*mc_have_z;
  const double A = cfg_readDouble (fp);

  const double kx = 2*mc_pi*mx/Lx;
  const double ky = 2*mc_pi*my/Ly;
  const double kz = 2*mc_pi*mz/Lz;

  ENSURE (abs (mx) + abs (my) + abs (mz) > 0, "Bad wave numbers (%d, %d, %d).", mx, my, mz);

  if (memEstimateOnly)
    return;

  long int N;
  marker_t *p = plasma_getObject (0, &N);					// Gets array of markers.
  for ( ; N > 0 ; --N, ++p)
  {
    const double phase = p->x*kx + p->y*ky + p->z*kz;
    p->vx += A*kx*cos (phase)*p->qDivM;
    p->vy += A*ky*cos (phase)*p->qDivM;
    p->vz += A*kz*cos (phase)*p->qDivM;
  }

  say ("tag_plasmaWave:");
  say ("  - plasma wave is exited in last object");
  say ("  - amplitude of the speed perturbation is %e [c],", A);
  say ("  - wave mode is (%d, %d, %d), k = (%e, %e, %e).", mx, my, mz, kx, ky, kz);
}
