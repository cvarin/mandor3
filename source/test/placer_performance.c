/** \file placer_performance.c
  * \brief Simple test of the different implementations of the particles exchange problem. Domain is populated by
  * particles and that sorter is invoked to place them according to their positions and current partition.
  *
  * This routine for now is called on loading and start/initializations.
  */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "type_marker.h"

#include "log.h"
#include "misc_partition.h"
#include "misc_cfgReader.h"
#include "misc_markerPlacer.h"
#include "misc_markerPacker.h"

#ifdef test_main
  ERROR: ONLY ONE TEST MAY BE PLUGGED!
#endif

#define test_main() test_placer()			///< To replace call at main_test.c to our entry function.
#define mc_testBufferSize	2000			///< Size of the test buffer (bytes).
#define mc_testFieldsNumber	11			///< Number of parameters we generate on top of the total DF to ensure conservation of everything.
#define mc_testRealRand		1			///< Set it to generate different sequence of random positions on every start.

// ---------------------------------------------------------------------------
/// Attempt to find the source of the SIGBUS crush in the a2a_sortCache().
// ---------------------------------------------------------------------------
static void
test_sigbusAligment (void)
{
  char cache[100];
  char *pos = cache + sizeof (markerPackHeader_t);
  marker_t *POS = (marker_t*) pos;
  POS->x = 10;
  say ("Y = %e", POS->y);
}

// ---------------------------------------------------------------------------
/// Fills domain with markers (returns ID).
// ---------------------------------------------------------------------------
static int
test_fill (int N)
{
  marker_cleanAll ();													// Removes all markers.
  int ID = marker_addChapter ();											// Gets chapter(ID) for markers.
  marker_t tmp[mc_markerPageN];

  srand (10);														// Resets repeatable way (for test).
#ifdef mc_testRealRand
  srand (MPI_Wtime ()*3.20451e5);											// Resets truly randomly.
#endif

  for (int i = 0 ; i < N ; i += mc_markerPageN)
  {
    int j = 0;
    for ( ; i + j < N && j < mc_markerPageN ; ++j)
    {
      tmp[j].x = (rand ()*Lx)/(double) RAND_MAX;
      tmp[j].y = (rand ()*Ly)/(double) RAND_MAX;
      tmp[j].z = (rand ()*Lz)/(double) RAND_MAX;
      tmp[j].vx = rand ()/(double) RAND_MAX - 0.5;
      tmp[j].vy = rand ()/(double) RAND_MAX - 0.5;
      tmp[j].vz = rand ()/(double) RAND_MAX - 0.5;
      tmp[j].rho = rand ()/RAND_MAX - 0.5;
    }
    marker_addExtMarkers (ID, j, 10, tmp);
  }

  marker_syncTypes ();

  return ID;
}

// ---------------------------------------------------------------------------
/// Configures domain using file 'test_placer.cfg' (returns number of particles, the rest goes to misc_partition.c, etc).
// ---------------------------------------------------------------------------
static int
test_configure (void)
{
  FILE *fp = cfg_open ("test_placer.cfg", "rt", __func__);								// Reads configuration.
  int N = cfg_readInt (fp);
  int imax = cfg_readInt (fp);
  int jmax = cfg_readInt (fp);
  int kmax = cfg_readInt (fp);
  int bc_xmin = cfg_readInt (fp);
  int bc_ymin = cfg_readInt (fp);
  int bc_zmin = cfg_readInt (fp);
  int bc_xmax = cfg_readInt (fp);
  int bc_ymax = cfg_readInt (fp);
  int bc_zmax = cfg_readInt (fp);
  fclose (fp);

  imax *= mc_have_x;													// Applies degeneration state.
  jmax *= mc_have_y;
  kmax *= mc_have_z;
  bc_xmin *= mc_have_x;
  bc_ymin *= mc_have_y;
  bc_zmin *= mc_have_z;
  bc_xmax *= mc_have_x;
  bc_ymax *= mc_have_y;
  bc_zmax *= mc_have_z;

  parameter_setTau (1e-4);												// Sets global parameters of the run.
  parameter_setSteps (1.0, 1.0, 1.0);
  parameter_setSize (imax, jmax, kmax);
  parameter_setBounds (bc_xmin, bc_xmax, bc_ymin, bc_ymax, bc_zmin, bc_zmax);
  parameter_info ();

  partition_init ();													// Partitions the domain.
  partition_show ();													// Shows partition.

  return N;
}

// ---------------------------------------------------------------------------
/// Collects information about particles to check that all particles are redistributed, but DF is absolutely the same.
/// Parameters are N, Q, W, P, R, R^2.
// ---------------------------------------------------------------------------
static void
test_particlesDF (double param[mc_testFieldsNumber])
{
  double tmp[mc_testFieldsNumber];

  memset (tmp, 0, mc_testFieldsNumber*sizeof (double));									// Cleans all accumulators.

  for (int ID = markerChapter_first () ; ID >= 0 ; markerChapter_next (&ID))
  {
    #pragma set woff 1343												// To avoid warning on SGI.
    markerIterator_t page;
    #pragma reset woff 1343

    for (markerPage_first (ID, &page) ; page.df ; markerPage_next (&page))						// Checks all pages of chapter.
    {
      marker_t *f = page.df;
      for (int p = 0 ; p < page.N ; p++)
      {
        double x = f[p].x, y = f[p].y, z = f[p].z;
        double vx = f[p].vx, vy = f[p].vy, vz = f[p].vz;
        tmp[0] += 1.0;
        tmp[1] += f[p].rho;
        tmp[2] += vx*vx + vy*vy + vz*vz;
        tmp[3] += vx;
        tmp[4] += vy;
        tmp[5] += vz;
        tmp[6] += x;
        tmp[7] += y;
        tmp[8] += z;
        tmp[9] += x*x + y*y + z*z;
        tmp[10] += f[p].rho/page.qDivM;
#if 10 >= mc_testFieldsNumber
  Compile time check of the range: too small array size.
#endif
      }
    }
  }

  MPI_Allreduce (tmp, param, mc_testFieldsNumber, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);					// Collects global parameters.
}

// ---------------------------------------------------------------------------
/// Checks that we can pack/unpack easily using misc_markerPacker.c library.
// ---------------------------------------------------------------------------
void
test_packer (void)
{
  char storage[mc_testBufferSize];											// Local buffer.
  markerPack_t buffer;

  buffer_init (&buffer, storage, mc_testBufferSize);
  int transferID = marker_addChapter ();

  double before[mc_testFieldsNumber];
  test_particlesDF (before);												// Hashes start parameters.

  say ("Transferring...");
  double startTime = MPI_Wtime ();

  for (int ID = markerChapter_first () ; ID >= 0 ; markerChapter_next (&ID))
  {
    #pragma set woff 1343												// To avoid warning on SGI.
    markerIterator_t page;
    #pragma reset woff 1343

    if (ID == transferID)												// Skips received particles.
      continue;

    for (markerPage_first (ID, &page) ; page.df ; markerPage_next (&page))						// Checks all pages of chapter.
    {
      marker_t *f = page.df;

      buffer_setExtType (&buffer, page.qDivM);										// Tells packer what is q/M.
      for (int p = 0 ; p < page.N ; ++p)
      {
        if (buffer_addMarker (&buffer, f + p))										// Packs particle.
        {
          buffer_unpackExt2chapter (&buffer, transferID);								// Unpacks data.
          buffer_rewind (&buffer);
        }
      }
      markerPage_setN (&page, 0);											// Updates number of particles left on page.
    }

    marker_defragChapter (ID);												// Defrags chapter (basically, removes pages).
  }

  buffer_finalize (&buffer, 0);												// Unpacks left-overs.
  buffer_unpackExt2chapter (&buffer, transferID);

  marker_syncTypes ();													// Syncronizes markers.

  double finishTime = MPI_Wtime (), after[mc_testFieldsNumber];
  test_particlesDF (after);												// Hashes finish parameters.

  say ("Packing/unpacking of %d particles per cpu (%d cpus total) takes %.3e second(s).", (int) before[0], cpu_total, finishTime - startTime);

  say ("before    after    delta");
  for (int i = 0 ; i < mc_testFieldsNumber ; ++i)
    say ("%+.2e %+.2e, %+.2e", before[i], after[i], before[i] - after[i]);
}

// ---------------------------------------------------------------------------
/// Entry point.
// ---------------------------------------------------------------------------
static void
test_placer (void)
{
  int N;

  test_sigbusAligment ();

  say ("Start-up configuration...");
  N = test_configure ();												// Configures domain.

  say ("Fills domain with %d particles...", N);									// Fills domain by particles.
  test_fill (N);

  say ("Testing pack/unpack speed and reliability on %d particles...", N);						// Checks pack/unpack routines.
  test_packer ();

  double before[mc_testFieldsNumber];
  test_particlesDF (before);												// Hashes start parameters.

  double startTime = MPI_Wtime ();
  say ("Testing new placer (based on misc_socket.c)");
  placer_exchange ();

  double finishTime = MPI_Wtime (), after[mc_testFieldsNumber];
  test_particlesDF (after);												// Hashes finish parameters.

  say ("Sorting of %d particles per cpu (%d cpus total) takes %.3e second(s).", N, cpu_total, finishTime - startTime);

  say ("before      after       delta        %%     conclusion");
  for (int i = 0 ; i < mc_testFieldsNumber ; ++i)
    say ("%+.2e %+.2e %+.4e %+.2e%%    %s", before[i], after[i], before[i] - after[i], (before[i] - after[i])*100.0/before[i],
              (fabs ((before[i] - after[i])*100.0/before[i]) > 1e-13*sqrt (N)) ? "bad" : "good");
}
