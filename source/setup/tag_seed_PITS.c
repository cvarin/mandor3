/** \file tag_seed_PITS.c
  * Excites unstable (PiTS) modes for photo-ionized plasma. For now only
  * \f$\vec k \parallel \vec e_z\f$ solutions are used.
  *
  * Creates perturbation of the electric field, initial positions and velocities
  * of the particles using linear solution for the 2-stream instability of the
  * photo-electrons distribution function (see tag_photoelectrons ()). Parameters
  * of the spectrum of modes are written to the 'output/tag_seed_PITS.dat' in
  * tecplot form.
  *
  * \todo XXX Comment tag parameters.
  */

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <mpi.h>

#include "misc_units.h"
#include "log.h"
#include "misc_cfgReader.h"

#include "main.h"
#include "tag_photoelectrons.h"

static int modeN = 0;			///< Number of modes.
static int modeM1 = 0;			///< First mode number (so <b>modeM1 <= m<sub>z</sub> < modeM1 + modeN</b>).
static double dEz = 0;			///< Amplitude of the electric field of the perturbations.
static double *modePhase = NULL;	///< Phase of the mode (<b>can be random and must be syncronized across cpus</b>).
static double *modeGamma = NULL;	///< Growth rate in \f$ \omega_{pe} \f$ (theoretical, used in linear solution expressions).

// ---------------------------------------------------------------------------
/// \brief LHS of the dispersion equation (18) from photoDF paper:
/// \f$\displaystyle  3-{k^2v_0^2 \over \omega^2_{pe}}+{\cal R}(\beta)+{3\,\beta^2\over \beta^2-1}=0  \f$.
/// In my case \f$\alpha = \omega_{pe}/kV_0,\ s = \gamma/kV_0\f$.
// ---------------------------------------------------------------------------
static double
dispFunc_photoES (double s, double alpha)
{
  return 3.0 - 1.0/(alpha*alpha) - 6*s*atan (1.0/s) + 3.0*s*s/(s*s + 1.0);
}

// ---------------------------------------------------------------------------
/// \brief Solver of the dispersion equation (18). \f$ \alpha = \omega_{pe}/kV_0, \ s = \gamma/kV_0 \f$,
/// dispersion equation is taken from EMI paper (eq. (18), \f$ \theta = 0 \f$). Only
/// unstable modes are excited (\f$ 0 \leqslant k_z \leqslant \sqrt 3 \cdot \omega_{pe}/V_0 \f$).
// ---------------------------------------------------------------------------
static double
dispFunc_s (double alpha)
{
  int    iter = 0;
  double sUp = 5, sDn = 0, fUp;

  if (fabs (alpha) <= 1.0/sqrt (3.0))											// Searches only in the region of instability.
    return 0;

  fUp = dispFunc_photoES (sUp, alpha);

  do {
    double sMd = 0.5*(sDn + sUp), fMd = dispFunc_photoES (sMd, alpha);

    if (fUp*fMd <= 0)
    {
      sDn = sMd;
    }
    else
    {
      sUp = sMd;
      fUp = fMd;
    }
  } while (fabs (sUp - sDn) > 1e-14*fabs (sUp + sDn) && iter++ < 50000);

  return 0.5*(sUp + sDn);
}

// ---------------------------------------------------------------------------
/// Initializes spectrum using tag_photoDF_parameters(int param, void *pntr) to get parameter of the photo-DF.
// ---------------------------------------------------------------------------
// XXX static
void
modes_init (int m1, int m2, double gamma_cutOff, double phase)
{
  int    mz;
  double V0, omega2_pe; 												// Parameters of the photo-DF.
  char name[50];

  tag_photoDF_parameters (mc_photoDF_V0, &V0);										// Takes photo-DF parameters from creator.
  tag_photoDF_parameters (mc_photoDF_omega2_pe, &omega2_pe);

  modeN = m2 - m1 + 1;
  modeM1 = m1;
  modePhase = (double*) calloc (modeN, sizeof (double));
  modeGamma = (double*) calloc (modeN, sizeof (double));

  const double unit_alpha = units (mc_r0)/(units (mc_t0)*units (mc_v0));						// [omega/kV] with k = 2\pi/\lambda.
  const double omega_pe = sqrt (omega2_pe);
  for (mz = m1 ; mz <= m2 ; mz++)
  {
    double kz, gammaDivOmegaPE;												// Parameters of the mode.
    kz = 2.0*mc_pi*mz/Lz;
    gammaDivOmegaPE = kz*V0/(unit_alpha*omega_pe)*dispFunc_s (unit_alpha*omega_pe/(kz*V0));

    if (gammaDivOmegaPE < gamma_cutOff)
      continue;

    modeGamma[mz-modeM1] = gammaDivOmegaPE;
    modePhase[mz-modeM1] = (phase < 0) ? 2.0*mc_pi*rand ()/(double) RAND_MAX : phase;
  }

  MPI_Bcast (modePhase, modeN, MPI_DOUBLE, 0, MPI_COMM_WORLD);								// Syncronizes global random parameters.

  sprintf (name, "output/tag_seed_PITS(%03d).dat", cpu_here);
  FILE *fp = cfg_open (name, "wt", "tag_seed2Stream_singleMode");
  fprintf (fp, "variables = mode, E<sub>0</sub>, <greek>f/p</greek>, <greek>g/w</greek><sub>pe</sub>, ");
  fprintf (fp, "<greek>g</greek><sub>WE</sub>, k<sub>z</sub>, \"2.0/(h3*kz)*sin (kz*h3/2.0) - 1.0\"\n");
  fprintf (fp, "zone t = \"L_z = %e, `w_p_e = %e\", f = point\n", Lz, sqrt (omega2_pe));
  for (mz = m1 ; mz <= m2 ; mz++)											// Prints spectrum of initial seed.
  {
    double kz, gamma;													// Parameters of the mode.
    kz = 2*mc_pi*mz/Lz;
    gamma = modeGamma[mz-modeM1];
    fprintf (fp, "%d %e %e %e %e %e %e\n", mz, dEz, modePhase[mz-modeM1]/mc_pi, gamma, 2.0*gamma*omega_pe, kz, 2.0/(h3*kz)*sin (kz*h3/2.0) - 1.0);
  }
  fclose (fp);
}

// ---------------------------------------------------------------------------
/// Function supposes that last cached chapter is photo-electrons with no perturbations and
/// adds perturbation with given spectral properties (to study electrostatic two-stream instability).
// ---------------------------------------------------------------------------
void
tag_seed_PITS (FILE *fp, meshVec_p E)
{
#if 0
  int    mz, mz_start, mz_stop, ID;
  double phase0, gamma_cutOff;
  markerIterator_t page;

  dEz = cfg_readDouble (fp);												// Reads parameters.
  mz_start = cfg_readInt (fp);
  mz_stop  = cfg_readInt (fp);
  gamma_cutOff = cfg_readDouble (fp);
  phase0 = cfg_readDouble (fp)*mc_pi;

  if (mz_start > mz_stop)
    error ("tag_seed2Stream: bad spectrum range\n.");

  modes_init (mz_start, mz_stop, gamma_cutOff, phase0);									// Initializes wave spectrum.

  say ("%s: ", __func__);											// Prints info.
  say ("  - spectrum of the unstable electrostatic modes is excited,");
  say ("  - detailes are written to 'output/tag_seed2Stream.dat',");
  say ("  - modes from %d to %d,", mz_start, mz_stop);
  if (10*abs (mz_start) >= dmn_max[2] - dmn_min[2] || 10*abs (mz_stop) >= dmn_max[2] - dmn_min[2])
    msg_warning (mf_here, "bad resolutions for modes >= %d.", (dmn_max[2] - dmn_min[2])/10);
  if (phase0 < 0)
    say ("  - phases are random,");
  else
    say ("  - phases are fixed to %f \\pi,", phase0/mc_pi);
  say ("  - amplitude %e,", dEz);
  say ("  - gamma-cutoff is %e.", gamma_cutOff);

  if (memEstimateOnly)
    return;

  tag_photoDF_parameters (mc_photoDF_ID, &ID);										// Gets ID to access particles.
  if (ID < 0)
    error ("tag_seed2stream: photo-DF component is not created yet.");

  cache_flush ();													// Flushes cache before doing any updates.
  double dzMax = -1, dzMin = +1;
  double dvMax = -1e10, dvMin = +1e10;
/*
  FILE *seedFile = cfg_open ("output/tag_seed_PITS.dat", "wt", __func__);
  fprintf (seedFile, "variables = z, v_z, `dz, `dv_z, `dz_n_u_m, `dv_z_n_u_m\n");*/
  for (markerPage_first (ID, &page) ; page.df ; markerPage_next (&page))						// Adds perturbation of DF.
  {
    int p;
    marker_t *f = page.df;
    for (p = 0 ; p < page.N ; ++p)
    {
      const double oldZ = f[p].z;
      const double oldVZ = f[p].vz;
      double newZ = oldZ, newVZ = oldVZ;
      for (mz = mz_start ; mz <= mz_stop ; ++mz)									// Adds perturbation for each mode.
      {
        double kz, gamma, phase;											// Parameters of the mode.
        double kzz, kzvz, gamma2, dv, dz;										// Shortcuts for derivatives.

        if (!mz)													// There are no zero mode.
          continue;

        kz = 2*mc_pi*mz/Lz;												// Extracts mode form spectrum arrays.
        gamma = modeGamma[mz-modeM1];
        phase = modePhase[mz-modeM1];

        if (gamma < gamma_cutOff)
          continue;

        kzz = oldZ*kz;													// Adds perturbation of velocity and position.
        kzvz = oldVZ*kz/(2*mc_pi);
        gamma2 = gamma*gamma;
        dv = - page.qDivM*dEz*(kzvz*cos (kzz + phase) - gamma*sin (kzz + phase))/(2*mc_pi*(gamma2 + kzvz*kzvz));
        dz = page.qDivM*dEz*((gamma2 - kzvz*kzvz)*sin (kzz + phase) - 2*gamma*kzvz*cos (kzz + phase))/(4*mc_pi*mc_pi*(gamma2 + kzvz*kzvz)*(gamma2 + kzvz*kzvz));

        newZ += dz;
        newVZ += dv;
      }

      dzMax = (dzMax < newZ - oldZ) ? newZ - oldZ : dzMax;
      dzMin = (dzMin > newZ - oldZ) ? newZ - oldZ : dzMin;
      dvMax = (dvMax < newVZ - oldVZ) ? newVZ - oldVZ : dvMax;
      dvMin = (dvMin > newVZ - oldVZ) ? newVZ - oldVZ : dvMin;

/*      double z = f[p].z, vz = f[p].vz;
*/
      f[p].z = fmod (newZ + 4*Lz, Lz);
      f[p].vz = newVZ;
/*
      fprintf (seedFile, "%e %e %e %e %e %e\n", z, vz, f[p].z - z, f[p].vz - vz, f[p].z - z, f[p].vz - vz);*/
    }
  }
/*  fclose (seedFile);
*/
  for (mz = mz_start ; mz <= mz_stop ; mz++)										// Adds perturbation of electric field.
  {
    int i, j, k;
    double kz, gamma, phase;												// Parameters of the mode.

    if (!mz)														// There are no zero mode.
      continue;

    kz = 2*mc_pi*mz/Lz;													// Extracts mode from spectrum arrays.
    gamma = modeGamma[mz-modeM1];
    phase = modePhase[mz-modeM1];

    if (gamma < gamma_cutOff)
      continue;

    for (i = E->imin ; i <= E->imax ; i++)										// Adds perturbation of electric field.
      for (j = E->jmin ; j <= E->jmax ; j++)
        for (k = E->kmin ; k <= E->kmax ; k++)
        {
          const double kzz = (k - 0.5)*kz*h3;
          mv_fz(E, i, j, k) += dEz*sin (kzz + phase);
        }
  }

  say ("  - %.3e <= dz <= %.3e (h3 = %.3e),", dzMin, dzMax, h3);
  say ("  - %.3e <= dvz <= %.3e,", dvMin, dvMax);
#endif
  assert (0);
}
