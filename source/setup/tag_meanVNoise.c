/** \file tag_meanVNoise.c
  * \brief Creates noise perturbation of the mean velocity. Used as cheap way to introduce initial electrostatic (ES) perturbation when
  * exact expression for corresponding ES wave is complicated or unknown (typically it is used on the first stages of the new plasma
  * distribution function research).
  *
  * Example of the config file entry is
    <pre>
    [<V> perturbation]
    @ SpecDomain:shell          Spectral region (can be shell or (after implementation) region).
    > 1                         Parameters of the spectral region (for shell it is rMin, rMax in wave-numbers).
    > 6
    @ Direction:chaotic         Direction of the velocity (can be 'chaotic' or 'aligned').
    @ 0.01e0                    Vx
    @ 0.00e0                    Vy
    @ 0.00e0                    Vz
    @ Touch:last                Plasma component to perturb (can be 'all' or 'last').

      Introduces perturbation of the mean velocity of the chosen component(s) with given spectral parameters of the noise.
    </pre>
  *
  * \todo Choose random number generator type (clib or ranmar).
  */

#include <stdlib.h>
#include <assert.h>

#include <mpi.h>

#include "log.h"
#include "misc_PIC.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"
#include "misc_markerPlacer.h"

#include "setup/main.h"			// For memEstimateOnly flag.

// ---------------------------------------------------------------------------
/// Function to initialize vector.
// ---------------------------------------------------------------------------
typedef void (*dV_func) (vec3D_t *v);

int    noise_specDomain;		///< Type of the spectral domain of the noise (0 - spherical shell, 1 - region and so on).
int    noise_minR;			///< \b Min radius of the shell.
int    noise_maxR;			///< \b Max radius of the shell.

int    noise_vDir;			///< Type of the velocity direction (chaotic/aligned/curlFree).
double noise_dV[3];			///< Amplitude of the velocity perturbation.

int    noise_touchAll;			///< Flag to trigger processing of the all particles or last chapter only.

// XXX static meshVec_t     dV;		///< \f$ \delta \vec V (\vec r) = \sum_{\vec k} \delta \vec V_{\vec k}\cdot \cos (\vec k\cdot\vec r + \phi_{\vec k}) \f$.
// XXX static meshVec_t     spectrV;		///< \f$ \delta \vec V_{\vec k} \f$.
// XXX static meshDouble_t  spectrPhi;		///< \f$ \phi_{\vec k} \f$.

// ---------------------------------------------------------------------------
/// Generates random number in [-1, 1] interval.
// ---------------------------------------------------------------------------
static double
dV_rand (void)
{
  return rand ()*2.0/(double) RAND_MAX - 1;
}

// ---------------------------------------------------------------------------
/// Generates random vector, aligned to noise.dV.
// ---------------------------------------------------------------------------
// XXX static
void
dV_randAligned (vec3D_t *v)
{
  const double ampl = dV_rand ();
  v->x = ampl*noise_dV[0];
  v->y = ampl*noise_dV[1];
  v->z = ampl*noise_dV[2];
}

// ---------------------------------------------------------------------------
/// Generates random vector, aligned to noise_dV.
// ---------------------------------------------------------------------------
// XXX static
void
dV_randChaotic (vec3D_t *v)
{
  v->x = dV_rand ()*noise_dV[0];
  v->y = dV_rand ()*noise_dV[1];
  v->z = dV_rand ()*noise_dV[2];
}

// ---------------------------------------------------------------------------
/// Reads parameters of the mean velocity noise and excites perturbation requested.
// ---------------------------------------------------------------------------
void
tag_meanVNoise (FILE *fp)
{
#if 0
  say ("tag_meanVNoise:");

  const char *arg = cfg_readWord (fp);											// Reads type of the spectral domain.
  switch (noise_specDomain = cfg_identifyWord (arg, "SpecDomain:shell", 0, "SpecDomain:reg", 1, mc_cfgTermGuesses))
  {
    case 0:
      noise_minR = cfg_readOptInt (fp);
      noise_maxR = cfg_readOptInt (fp);
      if (noise_minR <= 0 || noise_maxR <= noise_minR)
        error ("tag_meanVNoise: bad k-shell; %d <= |k| <= %d.", noise_minR, noise_maxR);
      mesh_allocate (mcast_mesh (&spectrV), -noise_maxR, -noise_maxR, -noise_maxR, noise_maxR, noise_maxR, noise_maxR, "meanDV:dV_k", mc_vec3D_t);
      mesh_allocate (mcast_mesh (&spectrPhi), -noise_maxR, -noise_maxR, -noise_maxR, noise_maxR, noise_maxR, noise_maxR, "meanDV:phi_k", mc_double);
      say ("  - spectrum domain is shell-type, rMin = %d, rMax = %d", noise_minR, noise_maxR);
    break;

    default:
      error ("tag_meanVNoise: unimplemented spectral domain type '%s'.", arg);
    break;
  }

  arg = cfg_readWord (fp);												// Reads type of orientation field.
  if (-1 == (noise_vDir = cfg_identifyWord (arg, "Direction:chaotic", 1, "Direction:aligned", 0, mc_cfgTermGuesses)))
    error ("tag_meanVNoise: unimplemented type of direction field '%s'.", arg);
  say ("  - speed vector alignment is %s", (noise_vDir == 1) ? "chaotic" : "aligned");
  dV_func initFunc[2] = {dV_randAligned, dV_randChaotic};								// Function switch.
  dV_func func = initFunc[noise_vDir];

  noise_dV[0] = cfg_readDouble (fp);											// Reads velocity perturbation amplitude.
  noise_dV[1] = cfg_readDouble (fp);
  noise_dV[2] = cfg_readDouble (fp);
  say ("  - amplitude of the speed vector is (%e, %e, %e)", noise_dV[0], noise_dV[1], noise_dV[2]);

  arg = cfg_readWord (fp);												// Reads type of orientation field.
  if (-1 == (noise_touchAll = cfg_identifyWord (arg, "Touch:last", 0, "Touch:all", 1, mc_cfgTermGuesses)))
    error ("tag_meanVNoise: unimplemented set of particles to process '%s'.", arg);

  if (! cpu_here)
  {
    FILE *fp = cfg_open ("output/meanVNoise.dat", "wt", __func__);
    fprintf (fp, "variables = kk, kj, ki, dVx, dVy, dVz, phase\nzone t=\"dV_k\", i = %d, j = %d, k = %d, f = point\n",
             spectrV.kmax - spectrV.kmin + 1, spectrV.jmax - spectrV.jmin + 1, spectrV.imax - spectrV.imin + 1);
    for (int i = spectrV.imin ; i <= spectrV.imax ; ++i)								// Inits spectrum harmonics.
      for (int j = spectrV.jmin ; j <= spectrV.jmax ; ++j)
        for (int k = spectrV.kmin ; k <= spectrV.kmax ; ++k)
        {
          int iFlag = 2*(i > 0) + (i == 0);										// Sign triplets.
          int jFlag = 2*(j > 0) + (j == 0);
          int kFlag = 2*(k > 0) + (k == 0);

          int r2 = i*i + j*j + k*k;
          if (iFlag*9 + jFlag*3 + kFlag > 9 + 3 + 1 && (noise_specDomain == 0 && r2 < noise_maxR*noise_maxR + 0.1 && r2 > noise_minR*noise_minR - 0.1))
          {
            func (&mv_v (&spectrV, i, j, k));
            mv_f (&spectrPhi, i, j, k) = dV_rand ()*mc_pi;
          }
          else
          {
            mv_fx (&spectrV, i, j, k) = mv_fy (&spectrV, i, j, k) = mv_fz (&spectrV, i, j, k) = 0;			// This \vec k is aliased => dropping it.
            mv_f (&spectrPhi, i, j, k) = 0;
          }
          fprintf (fp, "%d %d %d %.5e %.5e %.5e %.5e\n", k, j, i, mv_fx (&spectrV, i, j, k), mv_fy (&spectrV, i, j, k), mv_fz (&spectrV, i, j, k),
                   mv_f (&spectrPhi, i, j, k));
        }
    fclose (fp);
  }

  MPI_Bcast (spectrV.storage, spectrV.size, MPI_BYTE, 0, MPI_COMM_WORLD);						// Syncs noise spectrums over cluster.
  MPI_Bcast (spectrPhi.storage, spectrPhi.size, MPI_BYTE, 0, MPI_COMM_WORLD);

  if (memEstimateOnly)													// Returns in memory usage estimate mode.
  {
    say ("  - estimate mode, leaving ...");
    return;
  }

  say_doing ("  - allocating mesh ...");
  mesh_allocate (mcast_mesh (&dV), cpu_min[0] - mc_have_x, cpu_min[1] - mc_have_y, cpu_min[2] - mc_have_z,			// Allocates mesh for mean velocity field.
                                   cpu_max[0] + mc_have_x, cpu_max[1] + mc_have_y, cpu_max[2] + mc_have_z, "meanDV:dV", mc_vec3D_t);

  const double ax = 2*mc_have_x*mc_pi/(dmn_max[0] - dmn_min[0] + 1 - mc_have_x);						// Base wave vector components.
  const double ay = 2*mc_have_y*mc_pi/(dmn_max[1] - dmn_min[1] + 1 - mc_have_y);
  const double az = 2*mc_have_z*mc_pi/(dmn_max[2] - dmn_min[2] + 1 - mc_have_z);
  for (int i = dV.imin ; i <= dV.imax ; ++i)										// Adds all modes of perturbation.
    for (int j = dV.jmin ; j <= dV.jmax ; ++j)
      for (int k = dV.kmin ; k <= dV.kmax ; ++k)
      {
        say_doing ("  - creating dV field mesh at %d, %d, %d ...", i, j, k);
        double x = 0, y = 0, z = 0;
        for (int ki = spectrV.imin ; ki <= spectrV.imax ; ++ki)
          for (int kj = spectrV.jmin ; kj <= spectrV.jmax ; ++kj)
            for (int kk = spectrV.kmin ; kk <= spectrV.kmax ; ++kk)
            {
              const double wave = cos (i*ki*ax + j*kj*ay + k*kk*az + mv_f (&spectrPhi, ki, kj, kk));			/// \todo Replace by array precalc.
              x += mv_fx (&spectrV, ki, kj, kk)*wave;
              y += mv_fy (&spectrV, ki, kj, kk)*wave;
              z += mv_fz (&spectrV, ki, kj, kk)*wave;
            }
        mv_fx (&dV, i, j, k) = x;
        mv_fy (&dV, i, j, k) = y;
        mv_fz (&dV, i, j, k) = z;
      }

  say_doing ("  - redistributing particles ...");
  placer_exchange ();													// Sorts particles over cluster.

  int particles = 0;
  int lastID = cache_lastChapter ();											// Gets last chapter ID.
  for (int ID = markerChapter_first () ; ID >= 0 ; markerChapter_next (&ID))						// Lists all chapters.
  {
    if (!noise_touchAll && lastID != ID)										// Choses last/all chapter to update.
      continue;

    #pragma set woff 1343												// To avoid warning on SGI.
    markerIterator_t page;
    #pragma reset woff 1343
    static int pageN = 0;
    say_doing ("  - updating page %d ...", ++pageN);
    for (markerPage_first (ID, &page) ; page.df ; markerPage_next (&page))
    {
      marker_t *f = page.df;
      for (int p = 0 ; p < page.N ; ++p)
      {
        double sx = f[p].x/h1, sy = f[p].y/h2, sz = f[p].z/h3;
        int i = mf_PIC_double2int(sx), j = mf_PIC_double2int(sy), k = mf_PIC_double2int(sz);
        sx -= i;
        sy -= j;
        sz -= k;
        f[p].vx += mf_PIC_readX (&dV, i, sx, j, sy, k, sz);
        f[p].vy += mf_PIC_readY (&dV, i, sx, j, sy, k, sz);
        f[p].vz += mf_PIC_readZ (&dV, i, sx, j, sy, k, sz);
      }
      particles += page.N;
    }
  }

  say_doing ("  - %d particles is kicked (%s).", particles, (noise_touchAll) ? "all particles" : "last chapter");
#endif
  assert (0);
}
