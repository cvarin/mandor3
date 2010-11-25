/*
 * Zadaet tupoe volmushenie photoDF dlya neystoychivosti Weibelya. Parametri sleduushie:
 * ------------------------------------------------------------------------------------------
 * H0, dV0 --- amplitude of initial perturbation of field (and corresponding speed giving current J = rot H)
 * modeX --- mode of seed, phi0 --- phase shift of the seed.
 */

#include <stdio.h>
#include <assert.h>

#include "type_mesh.h"
#include "type_marker.h"

#include "log.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "setup/tag_photoelectrons.h"

#include "main.h"

/*
 * Adds perturbation to the photo-DF.
 * NOTE: PERTURBATION IS NOT A PROPER SOLUTION, JUST HELICAL FIELD!
 */
void
tag_seedWeibel (FILE *fp, meshVec_p H)
{
#if 0 // XXX
  int    p, i, j, k, modeX, ID;
  double chargeDensity, V0, H0, dV0, phi0, kx, omega2_pe;
#pragma set woff 1343										// To avoid warning about uninitialized const fields.
  markerIterator_t page;
#pragma reset woff 1343

  H0 = cfg_readDouble (fp);									// Reads parameters.
  modeX = cfg_readInt (fp);
  phi0 = mc_pi*cfg_readDouble (fp);

  cache_flush ();										// Sends all markers to storage (so we can access them).

  tag_photoDF_parameters (mc_photoDF_ID, &ID);							// Gets parameters of the "quiet start" created photoDF.
  tag_photoDF_parameters (mc_photoDF_V0, &V0);
  tag_photoDF_parameters (mc_photoDF_omega2_pe, &omega2_pe);

  markerPage_first (ID, &page);									// Calculates parameters based on input data.
  chargeDensity = omega2_pe*4*mc_pi*mc_pi/page.qDivM;
  kx = 2*mc_pi*modeX/Lx;
  dV0 = H0*kx/chargeDensity;

  for (i = H->imin ; i <= H->imax ; i++)							// Adds perturbation of the magnetic field.
    for (j = H->jmin ; j <= H->jmax ; j++)
      for (k = H->kmin ; k <= H->kmax ; k++)
        mv_fy(H, i, j, k) += H0*sin (kx*(i - 0.5)*h1 + phi0);

  for ( ; page.df ; markerPage_next (&page))							// Adds 'rot H' current to the mean velocity.
  {
    marker_t *f = page.df;
    for (p = 0 ; p < page.N ; p++)
      f[p].vz += dV0*cos (kx*f[p].x + phi0);
  }

  say ("tag_seedWeibel: ");
  say ("  - helical magnetic field is added");
  say ("  - rot H currents are exited");
  say ("  - H0 = %e, dV0 = %e, kx = %e (mx = %d).", H0, dV0, kx, modeX);
#endif
  assert (0);
}
