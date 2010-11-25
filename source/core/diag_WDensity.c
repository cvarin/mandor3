/** \file diag_WDensity.c
  * Outputs averaged energies (thermal, electric field, magnetic field and so on) into text file.
  *
  * This is one of the basic diagnostics to monitor energy distribution over the degrees of freedom system have.
  * Energy is dumped into local binary file and on exit all files are merged into single text file prepared for TecPlot.
  *
  * \note In order to simplify the usage of the \b log axises in both OpenDX and Tecplot initial point is not included.
  */

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "log.h"
#include "misc_MPItags.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

static char  wName[40];				///< File name, initialized in wDensity_prepare().
static FILE *wDatFile = NULL;			///< Output filestream.
static int   wContinue = -1;			///< Flag of the data writing mode: new file or continue old file.

// ---------------------------------------------------------------------------
/// Parallel merging of all data (done on exit).
// ---------------------------------------------------------------------------
static void
wDensity_finalize (void)
{
  int newFile, N, bcastN;
  double *data, *res, t1, t2;

  fclose (wDatFile);												// Reopens buffer file for reading.
  wDatFile = cfg_open (wName, "rb", "wDensity_finalize");

  fseek (wDatFile, 0, SEEK_END);										// Gets number of records buffered.
  N = ftell (wDatFile);
  ENSURE (N % (sizeof (double)*6) == 0,
          "bad size of the buffer file (%d != %f*6*sizeof(double)).",
          N, N/(6.0*sizeof (double)));
  N /= sizeof (double)*6;

  if (!N)													// Exits if content is flushed already.
  {
    fclose (wDatFile);
    return;
  }

  rewind (wDatFile);												// Gets to the start of the file.

  data = (double *) malloc (N*6*sizeof (double));
  ENSURE (data, "out of memory (%f Kb)", N*6*sizeof (double)/1024.0);
  fread (data, sizeof (double), 6*N, wDatFile);									// Transfer the content of the file into buffer.

  fclose (wDatFile);												// Removes temporary file.
  wDatFile = NULL;
  remove (wName);

  // XXX refactor all!
  t1 = data[0];													// Simple verification of the consistency.
  t2 = data[(N-1)*6];
  bcastN = N;
  MPI_Bcast (&bcastN, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&t1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&t2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ENSURE (N == bcastN,                          "inconsistent data size");
  ENSURE (t1 == data[0] && t2 == data[(N-1)*6], "inconsistent time signature");

  res = (double *) malloc (N*6*sizeof (double));
  ENSURE (res, "out of memory (%.2f Kb)", N*6*sizeof (double)/1024.0);

  MPI_Reduce (data, res, N*6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);						// Accumulates data.

  if (!cpu_here)												// Merges data on the master node.
  {
    int i;
    FILE *fp;

    fp = fopen ("binData/diag.dat", "rt");									// Defines mode: append or new file.
    newFile = (!fp) || (!wContinue);
    if (fp)
      fclose (fp);

    if (newFile)												// Opens file to send the output to.
    {
      fp = cfg_open ("binData/diag.dat", "wt", "wDensity_finalize");
      fprintf (fp, "variables = t, W_E, W_M, W_T_x, W_T_y, W_T_z, W_p_e_r_p, W_k_i_n, W_E_M, W_`S\nzone f=point\n");
    }
    else
    {
      fp = cfg_open ("binData/diag.dat", "at", "wDensity_finalize");
    }

    for (i = 0 ; i < N ; i++)
    {
      double time = res[i*6]/cpu_total, W_E = res[i*6+1], W_H = res[i*6+2], WTx = res[i*6+3], WTy = res[i*6+4], WTz = res[i*6+5];

      if (time < 0.5*tau)
        continue;

      fprintf (fp, "%e %e %e %e %e %e %e %e %e %e\n", time, W_E, W_H, WTx, WTy, WTz,
               WTy + WTz, WTx + WTy + WTz, W_E + W_H, W_E + W_H + WTx + WTy + WTz);
    }
    fclose (fp);
  }

  free (data);
  free (res);
}

// ---------------------------------------------------------------------------
/// Opens buffer file, sets name and so on.
// ---------------------------------------------------------------------------
void
wDensity_prepare (int startNew)
{
  static int submitted = 0;
  sprintf (wName, "binData/W(t).%03d.tmp", cpu_here);
  wDatFile = cfg_open (wName, "wb", "wDensity_prepare");

  if (!submitted)
  {
    ENSURE (!atexit (wDensity_finalize), "cannot submit 'wDensity_finalize'");
    submitted = 1;
  }

  wContinue = (startNew == 0);
}

// ---------------------------------------------------------------------------
/// Writes new data point into file (OS will cache it automatically, no need to manually cache it).
// ---------------------------------------------------------------------------
void
wDensity_addPoint (double W_E, double W_M, double WTx, double WTy, double WTz)
{
  double data[6] = {Time, W_E, W_M, WTx, WTy, WTz};

  ENSURE (wDatFile, "diagnostic is not initialized");

  fwrite (data, sizeof (double), 6, wDatFile);

  ENSURE (!isnan (W_E + W_M + WTx + WTy + WTz),
          "nan detected: W_E = %e, W_H = %e, WTx = %e, WTy = %e, WTz = %e",
          W_E, W_M, WTx, WTy, WTz);
}

// ---------------------------------------------------------------------------
/// Merges all buffers and adds result to the data file.
// ---------------------------------------------------------------------------
void
wDensity_flushData (void)
{
  ENSURE (wDatFile, "diagnostic is not initialized");
  wDensity_finalize ();
  wDensity_prepare (0);
}
