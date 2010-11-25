/** \file main.c
  * \brief This diagnostic reads binary data from probe file(s) and
  *        translates it into tecplot file.
  * \todo Remove MPI from here.
  */

/** \mainpage Probe diagnostic of "Mandor" package.
  * This executable just reads binary files written by simulation
  * module and produces single TecPlot data file.
  *
  * Implementation details are documented in the diag_probe.h file.
  */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>						// Complex should be _before_ fftw.

#include <fftw3.h>

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "core/diag_probe.h"

static double uE0, uMicron, uFemtosecond;			///< Units of the base parameters.

static double windowWidth = 0;					///< The only parameter of the \b window wavelet - width.

static int           sourceN = 0;				///< Number of elements in allocated fftw tmp array.
static fftw_complex *source = NULL;				///< Buffer for fftw tmp array for doing transformations.
static int           spectrN = 0;				///< Number of elements in allocated fftw tmp array.
static fftw_complex *spectr = NULL;				///< Buffer for fftw tmp array for doing transformations.
static fftw_plan     fftwPlan = NULL;				///< Plan to use for transformations.

// ---------------------------------------------------------------------------
/// Updates plan used to do fftw transformation and reallocates all arrays (use N <= 0 to deallocate everything back).
// ---------------------------------------------------------------------------
static void
fft_updatePlan (const int srcN, const int specN)
{
  if (spectrN == specN || sourceN == srcN)
    return;

  if (fftwPlan)
  {
    fftw_destroy_plan (fftwPlan);
    fftw_free (source);
    fftw_free (spectr);
  }

  spectrN = sourceN = 0;
  spectr = source = NULL;
  fftwPlan = NULL;

  if (srcN <= 0 || specN <=0)
    return;

  sourceN = srcN;
  source = fftw_malloc (sizeof (fftw_complex)*sourceN);

  spectrN = specN;
  spectr = fftw_malloc (sizeof (fftw_complex)*spectrN);

  fftwPlan = fftw_plan_dft_1d (spectrN, source, spectr, FFTW_FORWARD, FFTW_ESTIMATE);

  if (!spectr || !source || !fftwPlan)
    error ("fft_updatePlan: memory/plan allocation error.");
}

// ---------------------------------------------------------------------------
/// Sets up all parameters for the \b window wavelet transformation.
// ---------------------------------------------------------------------------
static void
window_setup (FILE *fp)
{
  windowWidth = cfg_readOptDouble (fp);
  windowWidth *= (windowWidth < 0) ? -1 : uFemtosecond;
  say ("FFT window width is set to %e [t0] = %e [femtosecond]", windowWidth, windowWidth/uFemtosecond);
}

// ---------------------------------------------------------------------------
/// Sets up all parameters for the \b window wavelet transformation.
// ---------------------------------------------------------------------------
static void
window_transform (FILE *fp, double tau, int N, probeSample_t *array)
{
  const int windN = (windowWidth/tau < N) ? windowWidth/tau : N;

  if (windN < 1)
    error ("window_transform: bad parameters width (%e) and tau (%e).", windowWidth, tau);

  fft_updatePlan (windN, windN);											// Updates transformation plan.
  fprintf (fp, "zone t = \"window/%.3f[T]\", i = %d, j = %d\n", windN*tau, windN/2, (N - windN + 1)/10 + 1);		// Creates tecplot header.

  for (int j = 0 ; j <= N - windN ; j += 10)										// Scans samples.
  {
    for (int i = 0 ; i < windN ; ++i)											/// \todo Optimize fill.
      source[i] = array[i+j].Ey;

    fftw_execute (fftwPlan);

    for (int i = 0 ; i < windN/2 ; ++i)
      fprintf (fp, "%d %d %.3e\n", i, j, creal (spectr[i])*creal (spectr[i]) + cimag (source[i])*cimag (source[i]));
  }
}

// ---------------------------------------------------------------------------
/// Entry point.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
  parameter_enterMPI (argc, argv, 0, NULL);										// Initializes log-file and MPI parameters.
  msg_setRefreshTime (0.1);

  units_load ();													// Loads units to dimensionalize.
  uE0 = units (mc_E0)/units (mc_A0);
  uMicron = units (mc_micron);
  uFemtosecond = units (mc_femtosecond);

  FILE *fp = cfg_open ("diag_probes.cfg", "rt", __func__);								// Gets input parameters.
  const char *waveletType = cfg_readWord (fp);										// Reads wavelet type.
  const int wavelet = cfg_identifyWord (waveletType, "window", 0, mc_cfgTermGuesses);
  switch (wavelet)
  {
    case 0:
      window_setup (fp);
    break;

    case -1:
    default:
      error ("main: cannot recognize the type of the wavelet '%s'.", waveletType);
  }
  fclose (fp);

  probeGlobals_t globals;												// Gets global parameters.
  fp = cfg_open ("binData/probe_globals.bin", "rb", __func__);
  fread (&globals, sizeof (probeGlobals_t), 1, fp);
  fclose (fp);

  FILE *fpw = cfg_open ("output/wavelets.dat", "wt", __func__);								// Converts data into tecplot file format.
  fprintf (fpw, "variables = i, j, W\n");

  fp = cfg_open ("output/probes.dat", "wt", "probe2tec.out");
  fprintf (fp, "variables = t, E_x, E_y, E_z, H_x, H_y, H_z, E^2, H^2\n");

  int N;
  probeHeader_t header;
  probeSample_t *data = NULL;
  while ((N = probe_loadData (&header, &data)) != -1)
  {
    fprintf (fp, "zone t=\"(%.2f,%.2f,%.2f)[micron]\", f=point\n", header.i*globals.h1/uMicron, header.j*globals.h2/uMicron, header.k*globals.h3/uMicron);

    for (int t = 0 ; t < N ; ++t)											// Writes time series.
    {
      data[t].Ex *= uE0;
      data[t].Ey *= uE0;
      data[t].Ez *= uE0;
      data[t].Hx *= uE0;
      data[t].Hy *= uE0;
      data[t].Hz *= uE0;
      fprintf (fp, "%.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e %.9e\n",
        data[t].time, data[t].Ex, data[t].Ey, data[t].Ez, data[t].Hx, data[t].Hy, data[t].Hz,
        data[t].Ex*data[t].Ex + data[t].Ey*data[t].Ey + data[t].Ez*data[t].Ez,
	data[t].Hx*data[t].Hx + data[t].Hy*data[t].Hy + data[t].Hz*data[t].Hz);
    }

    say ("Probe @ (%.2f,%.2f,%.2f)[micron] / [%d, %d, %d] node is prepared.",
              header.i*globals.h1/uMicron, header.j*globals.h2/uMicron, header.k*globals.h3/uMicron, header.i, header.j, header.k);

    switch (wavelet)
    {
      case 0:
        window_transform (fpw, globals.tau, N, data);
      break;

      case -1:
      default:
        error ("main: cannot recognize the type of the wavelet '%s'.", waveletType);
    }
  }
  fclose (fp);
  fclose (fpw);

  say ("\nDone.\n");

  return 0;
}
