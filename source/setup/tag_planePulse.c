/*
 * Creates EM wave with given wave vector (initially exspressed in wave numbers mx, my, mz),
 * amplitude of the field, polarisation and plane electric field belongs to (defined by vectors n and k).
 */

#include <math.h>
#include <stdio.h>

#include "type_mesh.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

/*
 * Returns variable clamped to given limits
 */
static double
clampToInterval (double x, const double xMin, const double xMax)
{
  if (x < xMin)
    return xMin;

  if (x > xMax)
    return xMax;

  return x;
}

void
tag_planePulse (FILE *fp, meshVec_p E, meshVec_p H)
{
  int    i, j, k;
  int    axis, p;
  double qMin, qMax, E0;

  double e[3] = {0, 0, 0};
  double h[3] = {0, 0, 0};
  double clampMax, step;
  const char axisName[3] = "XYZ";
  const char names[3][10] = {"PLANE", "CIRCULAR", "CIRCULAR"};

  axis = cfg_readInt (fp);
  qMin = cfg_readDouble (fp)*units (mc_micron);
  qMax = cfg_readDouble (fp)*units (mc_micron);
  E0 = cfg_readDouble (fp);
  p  = cfg_readInt (fp);

  ENSURE (axis >= 0 && axis < 3, "bad axis '%d'", axis);

  ENSURE (p >= 0 && p < 3, "bad polarization '%d'", p);

  e[(axis + 1)%3] = fabs (E0);
  h[(axis + 2)%3] = E0;

  clampMax = (axis == 0)*Lx + (axis == 1)*Ly + (axis == 2)*Lz;
  step = (axis == 0)*h1 + (axis == 1)*h2 + (axis == 2)*h3;

  ENSURE (qMin >= 0 && qMax <= clampMax - 0 && (qMax - qMin) >= 10*step,
          "bad qMin(%le)/qMax(%le) parameters", qMin, qMax);

  for (i = E->imin ; i <= E->imax ; i++)						// That will work for periodic boundary conditions only.
    for (j = E->jmin ; j <= E->jmax ; j++)
      for (k = E->kmin ; k <= E->kmax ; k++)
      {
        double q = (axis == 0)*i*h1 + (axis == 1)*j*h2 + (axis == 2)*k*h3;
        double qE = 2*mc_pi*(clampToInterval (q, qMin, qMax) - qMin)/(qMax - qMin);
        double qH = 2*mc_pi*(clampToInterval (q - step/2, qMin, qMax) - qMin)/(qMax - qMin);

        mv_fx(E, i, j, k) += (1.0 - cos (qE))*(e[0]*sin (qE) + (p != 0)*h[0]*cos (qE));
        mv_fy(E, i, j, k) += (1.0 - cos (qE))*(e[1]*sin (qE) + (p != 0)*h[1]*cos (qE));
        mv_fz(E, i, j, k) += (1.0 - cos (qE))*(e[2]*sin (qE) + (p != 0)*h[2]*cos (qE));

        mv_fx(H, i, j, k) += (1.0 - cos (qE))*(h[0]*sin (qH) - (p != 0)*e[0]*cos (qH));
        mv_fy(H, i, j, k) += (1.0 - cos (qE))*(h[1]*sin (qH) - (p != 0)*e[1]*cos (qH));
        mv_fz(H, i, j, k) += (1.0 - cos (qE))*(h[2]*sin (qH) - (p != 0)*e[2]*cos (qH));
      }

  say ("tag_planePulse:");
  say ("  - electromagnetic plane pulse along %c-axis is excited,", axisName[axis]);
  say ("  - localized in %e [lambda] <= %c <= %e  [lambda] region,", qMin, axisName[axis], qMax);
  say ("  - amplitude = %e,", E0);
  say ("  - polarization is %s.", names[p]);
}
