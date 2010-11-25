/** \file tag_twoStream.c
  * Two stream testing (see 'tag_twoStream.h' for details).
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "type_marker.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "main.h"
#include "setup/plasma.h"

// ---------------------------------------------------------------------------
/// Allocates two streams and fits speed to have classic growth rate as required.
// ---------------------------------------------------------------------------
double
tag_twoStream (FILE *fp)
{
  double gammaWE = cfg_readDouble (fp);						// Required growth rate.
  double qDivM = cfg_readDouble (fp);
  int Nx = cfg_readInt (fp);							// Number of particles.
  int Ny = cfg_readInt (fp);
  int Nz = cfg_readInt (fp);
  int mx = cfg_readInt (fp);
  int my = cfg_readInt (fp);
  int mz = cfg_readInt (fp);

  Nx = mc_have_x*(Nx - 1) + 1;							// Shuts down unused axices.
  Ny = mc_have_y*(Ny - 1) + 1;
  Nz = mc_have_z*(Nz - 1) + 1;
  mx *= mc_have_x;
  my *= mc_have_y;
  mz *= mc_have_z;

  int Nbeam = Nx*Ny*Nz;
  double rho = gammaWE*gammaWE/qDivM;

  double kx = 2*mc_pi*mx/Lx;							// Wave vector of the most unstable mode.
  double ky = 2*mc_pi*my/Ly;
  double kz = 2*mc_pi*mz/Lz;
  double kAbs = sqrt (kx*kx + ky*ky + kz*kz);

  double V = gammaWE*sqrt (3.0)/(2*kAbs);
  double vx = V*kx/kAbs;
  double vy = V*ky/kAbs;
  double vz = V*kz/kAbs;

  ENSURE (Nx > 0 && Ny > 0 && Nz > 0 && gammaWE > 0 && fabs (qDivM) <= 1, "Bad input parameters.");

  rho *= Lx*Ly*Lz/(h1*h2*h3*Nbeam);

  double dx = Lx/(double) Nx;
  double dy = Ly/(double) Ny;
  double dz = Lz/(double) Nz;

  if (memEstimateOnly)
    return 2.0*Nx*Ny*Nz*sizeof (marker_t);

  ENSURE (cpu_total == 1, "Not parallel yet.");

  plasma_newObject ();								// Starts plasma object.
  for (int i = 0 ; i < Nx ; ++i)
    for (int j = 0 ; j < Ny ; ++j)
      for (int k = 0 ; k < Nz ; ++k)
      {
        marker_t *marker = plasma_marker ();					// Pointer on uninitialized marker.
        marker->x = (i + 0.5)*dx;
        marker->y = (j + 0.5)*dy;
        marker->z = (k + 0.5)*dz;
        marker->vx = vx;
        marker->vy = vy;
        marker->vz = vz;
        marker->rho = rho;
        marker->qDivM = qDivM;

        marker = plasma_marker ();						// Pointer on uninitialized marker.
        marker->x = (i + 0.5)*dx;
        marker->y = (j + 0.5)*dy;
        marker->z = (k + 0.5)*dz;
        marker->vx = - vx;
        marker->vy = - vy;
        marker->vz = - vz;
        marker->rho = rho;
        marker->qDivM = qDivM;
      }

  rho /= Lx*Ly*Lz/(h1*h2*h3*Nbeam);
  say ("%s: ", __func__);
  say ("  - two cold beams are allocated (%d x %d x %d particles each)", Nx, Ny, Nz);
  say ("  - gamma_WE = %e [omega_0/2 pi]", gammaWE);
  say ("  - k = (%e, %e, %e) [omega_0/2 pi c]", kx, ky, kz);
  say ("  - mode = (%d, %d, %d)", mx, my, mz);
  say ("  - rho = %e, omega_pe/omega_0 = %e", 2*rho, sqrt (2*qDivM*rho)/(2*mc_pi));
  say ("  - V = (%e, %e, %e), |V| = %e", vx, vy, vz, V);

  return 2.0*Nbeam*sizeof (marker_t);
}
