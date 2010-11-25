/*
 * Zadaet DF f(v) ~ \delta (v_y)*\delta (v_x^2 + v_z^2 - v_0^2) metodom "quiet start" (prosto kol'tzo v prostranstve
 * skorostey, proshe nekuda. Parametri DF sleduushie:
 * ------------------------------------------------------------------------------------------
 * N - chislo chastitz v bazovom shablone;
 * nx, ny, nz - parametri razmesheniya chastitz v prostranstve;
 * V0 - radius in a speed-space.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "type_mesh.h"
#include "type_marker.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "main.h"
#include "setup/plasma.h"

static int    nx, ny, nz, N;
static double V0;

static double
tag_ringInit (double rho, double qDivM)
{
  double dx, dy, dz, number = 0;

  dx = h1/nx;
  dy = h2/ny;
  dz = h3/nz;

  plasma_newObject ();								// Starts new chapter.

  for (int i = dmn_min[0] ; i < (dmn_max[0] + 1 - mc_have_x)*nx ; ++i)
    for (int j = dmn_min[1] ; j < (dmn_max[1] + 1 - mc_have_y)*ny ; ++j)
      for (int k = dmn_min[2] ; k < (dmn_max[2] + 1 - mc_have_z)*nz ; ++k)
        for (int l = 0 ; l < N ; ++l)
        {
          marker_t *marker = plasma_marker ();					// Gets pointer on free marker.
          marker->x = (i + 0.5)*dx;
          marker->y = (j + 0.5)*dy;
          marker->z = (k + 0.5)*dz;
          marker->vx = V0*cos (2*mc_pi*l/N);
          marker->vy = 0;
          marker->vz = V0*sin (2*mc_pi*l/N);
          marker->rho = rho;
          marker->qDivM = qDivM;
          ++number;
        }

  return number;
}

double
tag_ringDF (FILE *fp)
{
  nx = cfg_readInt (fp);
  ny = cfg_readInt (fp);
  nz = cfg_readInt (fp);
  N  = cfg_readInt (fp);

  double omegaPeSquared = cfg_readDouble (fp);
  double qDivM = cfg_readDouble (fp);
  V0 = cfg_readDouble (fp);

  double chargeDensity = omegaPeSquared*4*mc_pi*mc_pi/qDivM;				// Calculates parameters based on input data.

  ENSURE (qDivM*chargeDensity > 0, "Opposite signs q/M (%.3e) and charge density (%.3e).", qDivM, chargeDensity);

  ENSURE (nx > 0 && ny > 0 && nz > 0 && N > 10,
             "Badly resolved distribution function: nx (%d), ny (%d), nz (%d), N (%d).", nx, ny, nz, N);

  ENSURE (cpu_total == 1, "Not parallel yet.");

  say ("tag_ringDF:");
  say ("  - ring-DF plasma component is added,");
  say ("  - V0 = %e, rho = %e (%e n_cr),", V0, chargeDensity, omegaPeSquared);
  say ("  - particles placement within a node is %d x %d x %d (%d particles per ring).", nx, ny, nz, N);

  if (memEstimateOnly)
    return 1.0*nx*ny*nz*N*sizeof (marker_t);

  return sizeof (marker_t)*tag_ringInit (chargeDensity/(double)(nx*ny*nz*N), qDivM);
}
