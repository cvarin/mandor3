/** \file test_FFT.c
  * Subroutines to test 1D/2D/3D FFT transformation for given N.
  *
  * |attention After version <b>/339</b> custom implementation of the fast fourier transformation was replaced
  * by the library <a href="http://www.fftw.org">fftw3</a>.
  */

#include <math.h>
#include <complex.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include <fftw3.h>

#include "log.h"
#include "misc_cfgReader.h"

// /**
//   * Checks that energy density for mode -k goes to mode N - k.
//   */
// double
// test_FFT3D (int Nx, int Ny, int Nz, int tryN, double *execTime)
// {
//   if (Nx <= 0 || Ny <= 0 || Nz <= 0)
//     error ("test_FFT: bad parameters Nx (%d), Ny (%d) or Nz (%d).", Nx, Ny, Nz);
//
//   complex double *data1 = (complex double*) malloc (sizeof (complex double)*Nx*Ny*Nz);
//   complex double *data2 = (complex double*) malloc (sizeof (complex double)*Nx*Ny*Nz);
//
//   for (int i = 0 ; i < Nx ; ++i)
//     for (int j = 0 ; j < Ny ; ++j)
//       for (int k = 0 ; k < Nz ; ++k)
//       {
//         data1[i*Ny*Nz + j*Nz + k] = (rand ()/(double) RAND_MAX - 0.5) + I*(rand ()/(double) RAND_MAX - 0.5);
//         data2[i*Ny*Nz + j*Nz + k] = data1[i*Ny*Nz + j*Nz + k];
//       }
//
//   double dt = MPI_Wtime ();
//   for (int i = 0 ; i < tryN ; ++i)
//   {
//     fft_transform3D (data1, Nx, Ny, Nz, -1);
//     fft_transform3D (data1, Nx, Ny, Nz, +1);
//   }
//   dt -= MPI_Wtime ();
//
//   double delta = 0;
//   for (int i = 0 ; i < Nx ; ++i)
//     for (int j = 0 ; j < Ny ; ++j)
//       for (int k = 0 ; k < Nz ; ++k)
//       {
//         double tmp = cabs (data1[i*Ny*Nz + j*Nz + k] - data2[i*Ny*Nz + j*Nz + k]);
//         delta = (delta < tmp) ? tmp : delta;
//       }
//
//   free (data1);
//   free (data2);
//
//   *execTime = - (double) dt/(double) tryN;
//   return delta/sqrt (tryN);
// }

// ---------------------------------------------------------------------------
/// Tests fftw3 library the same way local library was tested.
// ---------------------------------------------------------------------------
double
test_FFTW2 (int N, int tryN, double *execTime)
{
  if (N <= 0)
    error ("test_FFT: bad parameter N (%d).", N);

  fftw_complex *dataFFT = fftw_malloc (sizeof (fftw_complex)*N);
  fftw_complex *dataReference = fftw_malloc (sizeof (fftw_complex)*N);

  fftw_plan pForward = fftw_plan_dft_1d (N, dataFFT, dataFFT, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan pBackward = fftw_plan_dft_1d (N, dataFFT, dataFFT, FFTW_BACKWARD, FFTW_ESTIMATE);

  for (int i = 0 ; i < N ; ++i)
  {
    dataFFT[i] = (rand ()/(double) RAND_MAX - 0.5) + I*(rand ()/(double) RAND_MAX - 0.5);
    dataReference[i] = dataFFT[i];
  }

  const double normalizer = 1.0/N;
  double dt = MPI_Wtime ();
  for (int i = 0 ; i < tryN ; ++i)
  {
    fftw_execute (pForward);
    fftw_execute (pBackward);
    for (int i = 0 ; i < N ; ++i)
      dataFFT[i] = dataFFT[i]*normalizer;
  }
  dt -= MPI_Wtime ();

  double delta = 0;
  for (int i = 0 ; i < N ; i++)
  {
    double tmp = sqrt ((creal (dataFFT[i]) - creal (dataReference[i]))*(creal (dataFFT[i]) - creal (dataReference[i])) +
                       (cimag (dataFFT[i]) - cimag (dataReference[i]))*(cimag (dataFFT[i]) - cimag (dataReference[i])));
    delta = (delta < tmp) ? tmp : delta;
  }

  fftw_destroy_plan (pForward);
  fftw_destroy_plan (pBackward);
  fftw_free (dataFFT);
  fftw_free (dataReference);

  *execTime = - (double)dt/(double) tryN;										// Returns FFT time.
  return delta/(double) tryN;
}
