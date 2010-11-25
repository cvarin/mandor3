/** \file tag_jBalancer.c
  * \brief Cancels mean current by addition of the velocity to all particles (kind of changing frame in
  * \b CLASSICAL limit - no relativistic transformations for EM field (but there are one for velocity)).
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "type_marker.h"

#include "log.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "setup/main.h"
#include "setup/plasma.h"

void
tag_currentBalancer (FILE *fp)
{
  int iters = cfg_readInt (fp);
  ENSURE (iters > 1, "At least 2 iterations is normally required.");

  say ("%s:", __func__);
  if (memEstimateOnly)
  {
    say ("  estimate mode, leaving ...");
    return;
  }

  double totalVx = 0, totalVy = 0, totalVz = 0;
  for (int i = 0 ; i < iters ; ++i)
  {
    double jx = 0, jy = 0, jz = 0, charge = 0;

    marker_t dummy, *p = &dummy, *end;						// 'p = &dummy' is to enter cycle below.
    for (int count = 0 ; p ; ++count)
    {
      long int N;								// Checks if next plasma object is avaliable.
      for (p = plasma_getObject (count, &N), end = p + N ; p < end ; ++p)
      {
        double factor = p->rho/sqrt (1.0 + p->vx*p->vx + p->vy*p->vy + p->vz*p->vz);
        jx += p->vx*factor;
        jy += p->vy*factor;
        jz += p->vz*factor;
        charge += p->rho;
      }
    }

    double Vx = jx/charge, Vy = jy/charge, Vz = jz/charge;
    double V = Vx*Vx + Vy*Vy + Vz*Vz;
    double gamma_V = sqrt (1 + V);
    V = sqrt (V)/gamma_V;

    p = &dummy;									// Initialization to enter scanning cycle.
    for (int count = 0 ; p ; ++count)
    {
      long int N;								// Checks if next plasma object is avaliable.
      for (p = plasma_getObject (count, &N), end = p + N ; p < end ; ++p)
      {
        double v = p->vx*p->vx + p->vy*p->vy + p->vz*p->vz;
        double gamma_v = sqrt (1 + v);
        v = sqrt (v)/gamma_v;							// Relativistic speed 'gamma*v' goes to speed 'v'.

        p->vx = (p->vx/gamma_v - Vx/gamma_V)/(1 + v*V);				// Relativistic summation of 'v' + 'V'.
        p->vy = (p->vy/gamma_v - Vy/gamma_V)/(1 + v*V);
        p->vz = (p->vz/gamma_v - Vz/gamma_V)/(1 + v*V);

        gamma_v = sqrt (1 + p->vx*p->vx + p->vy*p->vy + p->vz*p->vz);		// Turns speed back to 'relativistic speed'.
        p->vx *= gamma_v;
        p->vy *= gamma_v;
        p->vz *= gamma_v;
      }
    }

    SAY_DEBUG ("  o iter #%d: dV = (% .3e, % .3e, % .3e)", i, Vx, Vy, Vz);
    totalVx += Vx;
    totalVy += Vy;
    totalVz += Vz;
  }

  say ("  DV = (% e, % e, % e)", totalVx, totalVy, totalVz);
}
