/** \file plasma_gaussTest.c
  * Functions to pick charge density before and after and test conservation and Gauss's laws.
  *
  * \attention It is DEBUGGING functions, they are simple (to drop mistakes) and they are SLOW.
  */

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "profiler.h"
#include "log.h"
#include "misc_partition.h"
#include "misc_cfgReader.h"

#include "core/plasma.h"

static int          gaussStarted = 0;					///< Flag to show that test is configured and running.
static mesh_t       gaussRho1_, gaussRho2_;				///< Actual storages, accessed through pointers below.
static meshDouble_p gaussRho1 = mcast_meshDouble (&gaussRho1_), 	///< Properly "typed" meshes.
                    gaussRho2 = mcast_meshDouble (&gaussRho2_);		//   XXX remove aliases after uniform meshes implementation.

static FILE *gaussRes = NULL, *gaussDump = NULL;			///< Files for max error and exact error distribution.

// Full "DEBUG" mode must include "TEST" output as well.
#ifdef MC_GAUSS_DEBUG
  #define MC_GAUSS_TEST
#endif

// ---------------------------------------------------------------------------
/// Closes files on exit (submitted to 'atexit(3)').
// ---------------------------------------------------------------------------
static void
gauss_stop (void)
{
  if (gaussRes)									// Closes dump files.
    fclose (gaussRes);
  if (gaussDump)
    fclose (gaussDump);
  gaussRes = gaussDump = NULL;
}

// ---------------------------------------------------------------------------
/// Prepares everything for testing (allocated meshes, etc).
// ---------------------------------------------------------------------------
static void
gauss_start (void)
{
   if (gaussStarted)								// No double initialization allowed.
      return;

   mesh_allocate (&gaussRho1_, cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2,
                               cpu_max[0] + 2, cpu_max[1] + 2, cpu_max[2] + 2, "plasma_rho1", mc_double);
   mesh_allocate (&gaussRho2_, cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2,
                               cpu_max[0] + 2, cpu_max[1] + 2, cpu_max[2] + 2, "plasma_rho2", mc_double);
   memset (gaussRho1_.storage, -1, gaussRho1_.size);				// NaNs to catch missed initialization.
   memset (gaussRho2_.storage, -1, gaussRho2_.size);

   ENSURE (atexit (gauss_stop) == 0,
           "cannot submit automatic test file finalizator");

   char *msg = "charge conservation test is active:\n"
               "  - run is slower due to extra blocking parallel exchanges"
#ifdef MC_GAUSS_DEBUG
               "\n  - HUGE dump is saved EACH TIME STEP."
#endif
   ;
  SAY_WARNING (msg);

#ifdef MC_GAUSS_DEBUG
  char name[40];
  sprintf (name, "output/gauss_dump.%d.dat", cpu_here);
  gaussDump = cfg_open (name, "wt", __func__);					// Saves domain size into the head of dump.
  fwrite (cpu_min, sizeof (int), 6, gaussDump);
#endif

  // Computes charge density (usually @ gauss_before we assume that 'rho2' contains valid data at the end of old timestep).
  plasma_rho (gaussRho2);

  if (!cpu_here)								// Only master maintains the report file.
  {
    gaussRes = cfg_open ("output/gauss_test.dat", "wt", __func__);
    fprintf (gaussRes, "variables = t, deltaRho_J, deltaRho_E\n");
  }

  gaussStarted = 1;								// Starts test.
}

// ---------------------------------------------------------------------------
/// First step - remembers charge density before time-step and current density calculation.
// ---------------------------------------------------------------------------
void
gauss_before (void)
{
#if !defined (MC_GAUSS_TEST)
  return;
#endif

  gauss_start ();								// Opens all files and allocates all meshes.

  profiler_endBegin (mc_prof_divJ_prep);					// Saves initial charge density.
  mf_mesh_copy (gaussRho2, gaussRho1);
}

// ---------------------------------------------------------------------------
/// Last test - compares charge density difference with current fluxes over cell's faces and 'div E' against '4*pi*rho'.
// ---------------------------------------------------------------------------
void
gauss_after (meshVec_RO_p J, meshVec_RO_p E)
{
  if (!gaussStarted)								// Quits if test is not running.
    return;

  profiler_endBegin (mc_prof_divJ);

  const int *mn = cpu_min, *mx = cpu_max;				// General sanity test.
  ENSURE (!mf_mesh_pointIsOutside(gaussRho1, mn[0], mn[1], mn[2]) &&
             !mf_mesh_pointIsOutside(gaussRho1, mx[0] + 1, mx[1] + 1, mx[2] + 1),
             "Bad mesh domain.");

  plasma_rho (gaussRho2);							// Computes final charge density.

  double delta[2] = {-1, -1}, final[2];						// Finds maximal deviation in {J, E}.
  for (int i = mn[0] ; i <= mx[0] ; ++i)
    for (int j = mn[1] ; j <= mx[1] ; ++j)
      for (int k = mn[2] ; k <= mx[2] ; ++k)
      {
        double tmp[2];
        tmp[0] = (mv_f(gaussRho2, i, j, k) - mv_f(gaussRho1, i, j, k))/tau +
                 (mv_fx(J, i+1, j, k) - mv_fx(J, i, j, k))/h1 +
                 (mv_fy(J, i, j+1, k) - mv_fy(J, i, j, k))/h2 +
                 (mv_fz(J, i, j, k+1) - mv_fz(J, i, j, k))/h3;
        tmp[1] = (mv_fx(E, i+1, j, k) - mv_fx(E, i, j, k))/h1 +
                 (mv_fy(E, i, j+1, k) - mv_fy(E, i, j, k))/h2 +
                 (mv_fz(E, i, j, k+1) - mv_fz(E, i, j, k))/h3 -
                 mv_f(gaussRho2, i, j, k);
        delta[0] = (fabs (tmp[0]) > delta[0]) ? fabs (tmp[0]) : delta[0];
        delta[1] = (fabs (tmp[1]) > delta[1]) ? fabs (tmp[1]) : delta[1];
#ifdef MC_GAUSS_DEBUG
        double data[5] = {tmp[0], tmp[1], mv_f(gaussRho1, i, j, k), mv_f(gaussRho2, i, j, k), tmp[1] + mv_f(gaussRho2, i, j, k)};
        fwrite (data, sizeof (double), 5, gaussDump);
#endif
      }

  MPI_Reduce (delta, final, 2, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);		// Collects final max over all cpus.
  if (!cpu_here)
    fprintf (gaussRes, "%e %e %e\n", Time, final[0], final[1]);
}
