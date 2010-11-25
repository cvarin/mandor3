/** \file spectr_process.c
  * \brief Processes data saved using spectr_dump.h routines (see spectr_process.h).
  */

#include <math.h>
#include <stdlib.h>
#include <complex.h>					// Complex should be included _before_ fftw to enjoy c99 complex types inside of the fftw.

#include <mpi.h>
#include <fftw3.h>

#include "IO_mesh.h"
#include "IO_names.h"

#include "type_mesh2.h"
#include "type_CFile.h"

#include "fft.h"

#include "log.h"
#include "misc_cfgReader.h"

static double  time;						///< Time of the record.
static double  dx, dy, dz;					///< Spatial steps.
static int     NKx, NKy, NKz;					///< Spectral region size is [-NKx, NKx] x .. x [-NKz, NKz].
static reg_t   fftMap = mc_reg_initBad;				///< Simulation domain with origin placed to the (0,0,0) node: all ranges are in [0, Nx), .., [0, Nz).
static reg_t   window = mc_reg_initBad;				///< Subdomain in computational domain from where to extract data for fft.

static fftw_plan planX = NULL, planY = NULL, planZ = NULL;	///< Plans to do in-place fft transformation in 'tmp1D' array.
static fftw_complex *source = NULL;				///< Temporary array for dumped data.
static fftw_complex *tmp1D = NULL;				///< Temporary 1D array for fftw transformation.

// ---------------------------------------------------------------------------
/// Updates all plans and temporary arrays if configuration is changed.
// ---------------------------------------------------------------------------
static void
spectr_updatePlans (void)
{
  static int oldMaxX = -1, oldMaxY = -1, oldMaxZ = -1;									// To remember if plans are created already.
  static int oldVolume = -1;
  if (oldMaxX == fftMap.max[0] && oldMaxY == fftMap.max[1] && oldMaxZ == fftMap.max[2] && 				// Avoids repetitions (planning is slow).
      source && oldVolume == reg_volume (&window))
  {
    return;
  }

  if (tmp1D)														// Cleans all.
  {
    fftw_free (source);
    fftw_free (tmp1D);
    fftw_destroy_plan (planX);
    fftw_destroy_plan (planY);
    fftw_destroy_plan (planZ);
  }

  tmp1D = source = NULL;
  planX = planY = planZ = NULL;
  oldMaxX = oldMaxY = oldMaxZ = -72623563;

  int N = (fftMap.max[1] > fftMap.max[0]) ? fftMap.max[1] : fftMap.max[0];						// Gets size along the longest direction.
  N = (fftMap.max[2] > N) ? fftMap.max[2] : N;
  ++N;															// +1 to get number of nodes from max index.

  tmp1D = fftw_malloc (N*sizeof (fftw_complex));									// Allocates tmp array.
  assert (tmp1D);

  planX = fftw_plan_dft_1d (fftMap.max[0] + 1, tmp1D, tmp1D, FFTW_FORWARD, FFTW_MEASURE);				// Designs plans.
  planY = fftw_plan_dft_1d (fftMap.max[1] + 1, tmp1D, tmp1D, FFTW_FORWARD, FFTW_MEASURE);
  planZ = fftw_plan_dft_1d (fftMap.max[2] + 1, tmp1D, tmp1D, FFTW_FORWARD, FFTW_MEASURE);
  assert (planX && planY && planZ);

  oldMaxX = fftMap.max[0];												// Remembers state.
  oldMaxY = fftMap.max[1];
  oldMaxZ = fftMap.max[2];

  source = fftw_malloc (sizeof (fftw_complex)*reg_volume (&window));							// Allocates memory for source data.
  assert (source);

  oldVolume = reg_volume (&window);
}

// ---------------------------------------------------------------------------
/// Loads all pieces of mesh to update entire \b domain (\b warning - mesh is initialized with \b NaN to help to catch unloaded regions).
// ---------------------------------------------------------------------------
static int
spectr_loadDump (const char *format, const char *component, int cpuN)
{
  memset (source, -1, sizeof (*source)*reg_volume (&window));								// To catch uninitialized vars.

  for (int cpu = 0 ; cpu < cpuN ; ++cpu)
  {
    CFile_t *file = CF_openRead (format, cpu);										// Erased files are checked earlier.
    assert (file);

    reg_t localMap;													// Reads region occupied by piece.
    CF_findChunk (file, "cpu region");
    CF_read (&localMap, sizeof (localMap), 1, file);

    reg_t overlap = window;												// If cpu piece is not in window - drop it.
    if (reg_overlap (&overlap, &localMap, NULL))
    {
      CF_close (file);
      continue;
    }

    mf_reg_collapse (&overlap);												// Prepares domain of loading.
    mf_mesh_initUnroll (&localMap);											// Prepares mapper.

    CF_findChunk (file, "%s", component);										// Opens proper chunk.
    const int span = (overlap.max[2] - overlap.min[2] + 1);								// Defines size of the continuous chunk.
    for (int i = overlap.min[0] ; i <= overlap.max[0] ; ++i)								// Loads piece.
      for (int j = overlap.min[1] ; j <= overlap.max[1] ; ++j)
      {
        double *tmp = (double*) tmp1D;											// Uses 1D buffer for loading.
        CF_seek (file, mf_mesh2_offset (&localMap, i, j, overlap.min[2])*sizeof (double), SEEK_SET);			// Positions pointer on the line's start.
        CF_read (tmp, sizeof (double), span, file);
        for (int k = overlap.min[2] ; k <= overlap.max[2] ; ++k)
          mv_map (&window, source, i, j, k) = *(tmp++);
      }
    CF_close (file);
  }

  return 0;
}

// ---------------------------------------------------------------------------
/// Converts temporary 3D complex array - in-place, only subregion of the k-space is constructed, data are rearranged to the origin of array.
// ---------------------------------------------------------------------------
static void
spectr_processDump (const char *formatIn, CFile_t *outputFile, const char *component, int cpuN)
{
  spectr_updatePlans ();
  spectr_loadDump (formatIn, component, cpuN);										// Loads data.

  const int Nx = fftMap.max[0] + 1, Ny = fftMap.max[1] + 1, Nz = fftMap.max[2] + 1;					// Local aliases for arrray sizes.

  for (int i = 0 ; i < Nx ; ++i)											// Transforms along Z.
    for (int j = 0 ; j < Ny ; ++j)
    {
      for (int k = 0 ; k < Nz ; ++k)											// Prepares data for fft.
        tmp1D[k] = mv_map (&fftMap, source, i, j, k);

      fftw_execute (planZ);												// Makes fft.

      int nz = 0;
      for (int k = Nz - NKz ; k < Nz ; ++k, ++nz)									// Gets first continuous part.
        mv_map(&fftMap, source, i, j, nz) = tmp1D[k];

      for (int k = 0 ; k <= NKz ; ++k, ++nz)										// Gets second continuous part.
        mv_map(&fftMap, source, i, j, nz) = tmp1D[k];
    }

  for (int i = 0 ; i < Nx ; ++i)											// Transforms along Y.
    for (int k = 0 ; k <= 2*NKz ; ++k)
    {
      for (int j = 0 ; j < Ny ; ++j)											// Prepares data for fft.
        tmp1D[j] = mv_map(&fftMap, source, i, j, k);

      fftw_execute (planY);												// Makes fft.

      int ny = 0;
      for (int j = Ny - NKy ; j < Ny ; ++j, ++ny)									// Gets first continuous part.
        mv_map(&fftMap, source, i, ny, k) = tmp1D[j];

      for (int j = 0 ; j <= NKy ; ++j, ++ny)										// Gets second continuous part.
        mv_map(&fftMap, source, i, ny, k) = tmp1D[j];
    }

  for (int j = 0 ; j <= 2*NKy ; ++j)											// Transforms along X.
    for (int k = 0 ; k <= 2*NKz ; ++k)
    {
      for (int i = 0 ; i < Nx ; ++i)											// Prepares data for fft.
        tmp1D[i] = mv_map(&fftMap, source, i, j, k);

      fftw_execute (planX);												// Makes fft.

      int nx = 0;
      for (int i = Nx - NKx ; i < Nx ; ++i, ++nx)									// Gets first continuous part.
        mv_map(&fftMap, source, nx, j, k) = tmp1D[i];

      for (int i = 0 ; i <= NKx ; ++i, ++nx)										// Gets second continuous part.
        mv_map(&fftMap, source, nx, j, k) = tmp1D[i];
    }

  const double scale = 4.0/(1.0*Nx*Nx*Ny*Ny*Nz*Nz);									// Normalization factor for after fftw.
  mv_map(&fftMap, source, NKx, NKy, NKz) = 0.5*mv_map(&fftMap, source, NKx, NKy, NKz);					// Zero mode is mirrored to itself.

  if (!CF_probeChunk (outputFile, "spectral domain size"))								// Saves region size.
  {
    CF_openChunk (outputFile, "spectral domain size");
    CF_write (&NKx, sizeof (int), 1, outputFile);
    CF_write (&NKy, sizeof (int), 1, outputFile);
    CF_write (&NKz, sizeof (int), 1, outputFile);
    CF_write (&time, sizeof (double), 1, outputFile);
    double L = (window.max[0] - window.min[0] + 1)*dx;									// Writes Lx of the transformed window.
    CF_write (&L, sizeof (double), 1, outputFile);
    L = (window.max[1] - window.min[1] + 1)*dy;										// Writes Ly of the transformed window.
    CF_write (&L, sizeof (double), 1, outputFile);
    L = (window.max[2] - window.min[2] + 1)*dz;										// Writes Lz of the transformed window.
    CF_write (&L, sizeof (double), 1, outputFile);
    CF_closeChunk (outputFile);
  }

  CF_openChunk (outputFile, "%s", component);										// Saves data in compact float format.
  for (int i = 0 ; i <= 2*NKx ; ++i)
    for (int j = 0 ; j <= 2*NKy ; ++j)
      for (int k = 0 ; k <= 2*NKz ; ++k)
      {
        // Some thoughts. Waves are assumed to be cos or sin shape so after FT they contribute to two places (+k and -k) due to Eiler formulas:
        // sin(a) = (exp(ia) - exp(-ia))/2i, cos(a) = (exp(ia) + exp(-ia))/2. For real input signal that means that amplitudes of inputs are
        // the same (only sign is different) and amplitude of the mode can be computed using direct square of either +k or -k contribution.
        complex double cmplx = mv_map (&fftMap, source, i, j, k);
        float tmp = (creal (cmplx)*creal (cmplx) + cimag (cmplx)*cimag (cmplx))*scale;
        CF_write (&tmp, sizeof (float), 1, outputFile);
      }
  CF_closeChunk (outputFile);
}

// ---------------------------------------------------------------------------
/// Prepares global parameters and returns \b cpuN (it is left as argument to ensure proper order
/// of function calls - you have to obtain cpuN to pass it down the chain).
// ---------------------------------------------------------------------------
static int
spectr_prepareDump (int recordNum, int kx, int ky, int kz, const reg_t *subdomain)
{
  CFile_t *file = CF_openRead ("tmp/spectrDump_%d_%d", recordNum, 0);							// Loads globals from master cpu file.
  if (!file)
    return -1;

  int found = (CF_findChunk (file, "global parameters") == 0);
  assert (found);

  int cpuN;
  char garbage[1024];
  garbage[1023] = 0;													// Simple check of array boundary violation.
  CF_scan (file, "%d%[^\n]", &cpuN, garbage);			assert (!garbage[1023]);
  CF_scan (file, "%d%[^\n]", &window.min[0], garbage);		assert (!garbage[1023]);
  CF_scan (file, "%d%[^\n]", &window.min[1], garbage);		assert (!garbage[1023]);
  CF_scan (file, "%d%[^\n]", &window.min[2], garbage);		assert (!garbage[1023]);
  CF_scan (file, "%d%[^\n]", &window.max[0], garbage);		assert (!garbage[1023]);
  CF_scan (file, "%d%[^\n]", &window.max[1], garbage);		assert (!garbage[1023]);
  CF_scan (file, "%d%[^\n]", &window.max[2], garbage);		assert (!garbage[1023]);
  CF_scan (file, "%le%[^\n]", &time, garbage);			assert (!garbage[1023]);
  CF_scan (file, "%le%[^\n]", &dx, garbage);			assert (!garbage[1023]);
  CF_scan (file, "%le%[^\n]", &dy, garbage);			assert (!garbage[1023]);
  CF_scan (file, "%le%[^\n]", &dz, garbage);			assert (!garbage[1023]);
  CF_close (file);

  mf_mesh_initUnroll (&window);												// Resets and checks window.

  int badSubdomain = reg_overlap (&window, subdomain, NULL);
  assert (!badSubdomain);

  --window.max[0];													// Excludes top boundary (it isn't saved too).
  --window.max[1];
  --window.max[2];

  fftMap = window;
  for (int axis = 0 ; axis < 3 ; ++axis)
  {
    fftMap.min[axis] -= window.min[axis];
    fftMap.max[axis] -= window.min[axis];
  }
  mf_mesh_initUnroll (&fftMap);
  mf_mesh_initUnroll (&window);

  NKx = mc_have_x*((kx >= (fftMap.max[0] + 1)/2 || kx < 0) ? (fftMap.max[0] + 1)/2 - 1 : kx);				// Clamps size of the spectral domain.
  NKy = mc_have_y*((ky >= (fftMap.max[1] + 1)/2 || ky < 0) ? (fftMap.max[1] + 1)/2 - 1 : ky);
  NKz = mc_have_z*((kz >= (fftMap.max[2] + 1)/2 || kz < 0) ? (fftMap.max[2] + 1)/2 - 1 : kz);

  return cpuN;
}

// ---------------------------------------------------------------------------
/// Converts EM-field distribution file into spectral energy density file.
// ---------------------------------------------------------------------------
void
spectr_convert (int recordNum, int kx, int ky, int kz, const reg_t *subdomain, int removeDumps)
{
  int cpuN = spectr_prepareDump (recordNum, kx, ky, kz, subdomain);							// Checks if file was processed and removed.
  if (cpuN < 0)
    return;

  int done = 0;
  CFile_t *outputFile;
  if ((outputFile = CF_openRead ("binData/spectr_%d", recordNum)))							// Checks if file is created already.
  {
    int nKx, nKy, nKz;													// Gets parameters of the
    int found = (CF_findChunk (outputFile, "spectral domain size") == 0);
    assert (found);
    CF_read (&nKx, sizeof (int), 1, outputFile);
    CF_read (&nKy, sizeof (int), 1, outputFile);
    CF_read (&nKz, sizeof (int), 1, outputFile);
    CF_close (outputFile);

    if (nKx == NKx && nKy == NKy && nKz == NKz)										// Returns if fft is done for the same sizes.
      done = 1;
  }

  char format[300];													// Creates spectral energy density files.
  sprintf (format,  "tmp/spectrDump_%d_%%d", recordNum);
  if (!done)
  {
    outputFile = CF_openWrite ("binData/spectr_%d", recordNum);
    assert (outputFile);

    spectr_processDump (format, outputFile, "E.x", cpuN);
    spectr_processDump (format, outputFile, "E.y", cpuN);
    spectr_processDump (format, outputFile, "E.z", cpuN);
    spectr_processDump (format, outputFile, "H.x", cpuN);
    spectr_processDump (format, outputFile, "H.y", cpuN);
    spectr_processDump (format, outputFile, "H.z", cpuN);
    CF_close (outputFile);
  }

  if (removeDumps)													// Removes dump data files if requested.
    for (int cpu = 0 ; cpu < cpuN ; ++cpu)
    {
      char name[200];
      sprintf (name, format, cpu);
      remove (name);
    }
}
