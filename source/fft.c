/** \file fft.c
  * \brief Fft of reg-array (see fft.h).
  */

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>		// Complex should be included _before_ fftw to enjoy c99 complex types inside of the fftw.
#include <fftw3.h>

#include "reg_unroller.h"

static fftw_plan     fft_X = NULL, fft_Y = NULL, fft_Z = NULL;	///< Plans to do in-place fft transformation in 'tmp1D' array.
static fftw_complex *fft_source = NULL;				///< Temporary array for 3D-data converted to complex.
static fftw_complex *fft_tmp = NULL;				///< Temporary 1D array for 1D fftw transformation.

// ---------------------------------------------------------------------------
/// Updates all plans and temporary arrays if configuration is changed.
// ---------------------------------------------------------------------------
static void
fft_initFFTW3 (const reg_t *dmn)
{
   static vec3i_t oldRegSizes = {{-1, -1, -1}};										// To remember if plans are created already.
   if (MF_VEC_EQ (oldRegSizes.r, dmn->max) && fft_source)					 			// Avoids repetitions (planning is expensive).
      return;

   // Cleans all.
   if (fft_tmp) {
      fftw_free (fft_source);
      fftw_free (fft_tmp);
      fftw_destroy_plan (fft_X);
      fftw_destroy_plan (fft_Y);
      fftw_destroy_plan (fft_Z);
   }

   fft_tmp = fft_source = NULL;
   fft_X = fft_Y = fft_Z = NULL;

   int N = (dmn->max[1] > dmn->max[0]) ? dmn->max[1] : dmn->max[0];							// Gets biggest size of domain.
   N = (dmn->max[2] > N) ? dmn->max[2] : N;
   ++N;															// Get 'number of nodes' (not 'max index').

   fft_tmp = fftw_malloc (N*sizeof (fftw_complex));									// Allocates tmp array.
   assert (fft_tmp);

   fft_X = fftw_plan_dft_1d (dmn->max[0] + 1, fft_tmp, fft_tmp, FFTW_FORWARD, FFTW_MEASURE);				// Makes fftw plans.
   fft_Y = fftw_plan_dft_1d (dmn->max[1] + 1, fft_tmp, fft_tmp, FFTW_FORWARD, FFTW_MEASURE);
   fft_Z = fftw_plan_dft_1d (dmn->max[2] + 1, fft_tmp, fft_tmp, FFTW_FORWARD, FFTW_MEASURE);
   assert (fft_X && fft_Y && fft_Z);

   MF_VEC_COPY (oldRegSizes.r, dmn->max);										// Remembers state.

   fft_source = fftw_malloc (sizeof (fftw_complex)*reg_volume (dmn));							// Allocates memory for source data.
   assert (fft_source);
}

// ---------------------------------------------------------------------------
/// Converts temporary 3D complex fft_source array (in-place, only subregion of the k-space is processed,
/// data are unrolled in original array).
// ---------------------------------------------------------------------------
static void
fft_doFFTW3 (const reg_t *dmn, const reg_t *spectr)
{
  assert (fft_source && fft_tmp);

  const int Nx = dmn->max[0] + 1, Ny = dmn->max[1] + 1, Nz = dmn->max[2] + 1;						// Local aliases for source sizes.
  const int NKx = spectr->max[0], NKy = spectr->max[1], NKz = spectr->max[2];						// Local aliases for spectr sizes.

  long int _s_[3], _o_;													// Gets unrolling coefficients for data.
  MF_UNROLL (_s_, _o_, dmn);
  assert (_o_ == 0);

  for (int i = 0 ; i < Nx ; ++i)											// Z-transforms: all i, all j.
    for (int j = 0 ; j < Ny ; ++j)
    {
      for (int k = 0 ; k < Nz ; ++k)											// Prepares data for fft.
        fft_tmp[k] = fft_source[MF_AIM(i, j, k, _s_, _o_)];

      fftw_execute (fft_Z);												// Makes fft along Z.

      int nz = 0;
      for (int k = Nz - NKz ; k < Nz ; ++k, ++nz)									// [-NKz, 0) part of spectrum.
        fft_source[MF_AIM(i, j, nz, _s_, _o_)] = fft_tmp[k];

      for (int k = 0 ; k <= NKz ; ++k, ++nz)										// [0, NKz] part of spectrum.
        fft_source[MF_AIM(i, j, nz, _s_, _o_)] = fft_tmp[k];
    }

  for (int i = 0 ; i < Nx ; ++i)											// Y-transforms: all i, kz in [-NKz, NKz].
    for (int k = 0 ; k <= 2*NKz ; ++k)
    {
      for (int j = 0 ; j < Ny ; ++j)											// Prepares data for fft.
        fft_tmp[j] = fft_source[MF_AIM(i, j, k, _s_, _o_)];

      fftw_execute (fft_Y);												// Makes fft along Y.

      int ny = 0;
      for (int j = Ny - NKy ; j < Ny ; ++j, ++ny)									// [-NKy, 0) part of spectrum.
        fft_source[MF_AIM(i, ny, k, _s_, _o_)] = fft_tmp[j];

      for (int j = 0 ; j <= NKy ; ++j, ++ny)										// [0, NKy] part of spectrum.
        fft_source[MF_AIM(i, ny, k, _s_, _o_)] = fft_tmp[j];
    }

  for (int j = 0 ; j <= 2*NKy ; ++j)											// X-transforms: ky,kz in subregion.
    for (int k = 0 ; k <= 2*NKz ; ++k)
    {
      for (int i = 0 ; i < Nx ; ++i)											// Prepares data for fft.
        fft_tmp[i] = fft_source[MF_AIM(i, j, k, _s_, _o_)];

      fftw_execute (fft_X);												// Makes fft along X.

      int nx = 0;
      for (int i = Nx - NKx ; i < Nx ; ++i, ++nx)									// [-NKx, 0) part of spectrum.
        fft_source[MF_AIM(nx, j, k, _s_, _o_)] = fft_tmp[i];

      for (int i = 0 ; i <= NKx ; ++i, ++nx)										// [0, NKx] part of spectrum.
        fft_source[MF_AIM(nx, j, k, _s_, _o_)] = fft_tmp[i];
    }

  fft_source[MF_AIM(NKx, NKy, NKz, _s_, _o_)] = fft_source[MF_AIM(NKx, NKy, NKz, _s_, _o_)]*0.5;			// Zero mode is mirrored to itself.
}

// ---------------------------------------------------------------------------
/// \brief Evaluates spectral energy density \f$|f_{\vec k}|^2\f$ where
/// \f$ f_{\vec r} = \sum_{\vec k} f_{\vec k}\cdot \exp(i\vec k\cdot\vec r)\f$. \b domain and \b data give
/// source, <b>spectr->max</b> holds desired spectral domain, <b>spectr->min</b> is ignored). Returns spectral
/// energy density in \b data and actual spectral domain in \b spectr (\b spectr can be used to access data).
// ---------------------------------------------------------------------------
void
fft (const reg_t *domain, double *data, reg_t *spectr)
{
  reg_t dmn = *domain;
  for (int axis = 0 ; axis < 3 ; ++axis)										// Updates local control structures.
  {
    dmn.max[axis] -= dmn.min[axis];											// Removes biase - only size matters.
    dmn.min[axis] -= 0;

    int Nk = spectr->max[axis], maxNk = dmn.max[axis]/2;								// Desired and max possible mode number.
    spectr->max[axis] = (Nk >= maxNk || Nk < 0) ? maxNk : Nk;
    spectr->min[axis] = - spectr->max[axis];
  }

  fft_initFFTW3 (&dmn);													// Allocates/prepares all.

  for (int i = 0 ; i < reg_volume (&dmn) ; ++i)										// Converts source to complex double.
    fft_source[i] = data[i];

  fft_doFFTW3 (&dmn, spectr);												// Does fft, result is in fft_source.

  long int _s1_[3], _o1_, _s2_[3], _o2_;
  MF_UNROLL (_s1_, _o1_, &dmn);												// Gets unrolling coefficients for data.
  MF_UNROLL (_s2_, _o2_, spectr);											// Gets unrolling coefficients for spectr.

  vec3i_t maxK;														// Alias for spectr->max.
  MF_VEC_COPY (maxK.r, spectr->max);
  const double scale = 4.0/(1.0*reg_volume (&dmn)*reg_volume (&dmn));							// Normalization factor.
  for (int i = 0 ; i <= 2*maxK.v.x ; ++i)										// Packs result back in 'data'.
    for (int j = 0 ; j <= 2*maxK.v.y ; ++j)
      for (int k = 0 ; k <= 2*maxK.v.z ; ++k)
      {
        complex double cmplx = fft_source[MF_AIM(i, j, k, _s1_, _o1_)];
        data[MF_AIM(i - maxK.v.x, j - maxK.v.y, k - maxK.v.z, _s2_, _o2_)] =
          (creal (cmplx)*creal (cmplx) + cimag (cmplx)*cimag (cmplx))*scale;
      }
}

// =================================================================================================== //

#if 0

// ---------------------------------------------------------------------------
/// Fills array with given mode, amplitude and phase.
// ---------------------------------------------------------------------------
static void
trainer_fill (const reg_t *dmn, double *data, int mx, int my, int mz, double A, double phase)
{
  const double kx = 2*acos (-1)/(dmn->max[0] - dmn->min[0] + 1);
  const double ky = 2*acos (-1)/(dmn->max[1] - dmn->min[1] + 1);
  const double kz = 2*acos (-1)/(dmn->max[2] - dmn->min[2] + 1);

  long int _s_[3], _o_;
  MF_UNROLL (_s_, _o_, dmn);												// Gets unrolling coefficients for source data.
  for (int i = dmn->min[0] ; i <= dmn->max[0] ; ++i)
    for (int j = dmn->min[1] ; j <= dmn->max[1] ; ++j)
      for (int k = dmn->min[2] ; k <= dmn->max[2] ; ++k)
        MV_CELL(data, i, j, k, _s_, _o_) += A*cos (i*kx*mx + j*ky*my + k*kz*mz + phase*acos(-1)/180.0);
}

// ---------------------------------------------------------------------------
/// \brief Tester of this library. Uncomment and recompile to use.
// ---------------------------------------------------------------------------
int
main (void)
{
  reg_t domain = {{0, 0, 0}, {10, 11, 12}};
  mf_reg_collapse (&domain);

  long int _s_[3], _o_;
  MF_UNROLL (_s_, _o_, &domain);											// Gets unrolling coefficients for source data.

  double *data = (double *) malloc (reg_volume (&domain)*sizeof (double));						// Allocates array.
  assert (data);

  while (1)
  {
    printf ("\n\n\nDimensionality: %dD\n", mc_have_x + mc_have_y + mc_have_z);
    printf ("Domain: %s\n", reg_printRanges (&domain));

    int mx, my, mz;
    printf ("Input wave numbers (mx, my, mz):\n");
    scanf ("%d %d %d", &mx, &my, &mz);

    double A, phase;
    printf ("Input amplitude and phase[degrees]:\n");
    scanf ("%le %le", &A, &phase);

    memset (data, 0, sizeof (double)*reg_volume (&domain));								// Clears array.
    trainer_fill (&domain, data, mx, my, mz, sqrt (A), phase);								// Fills mode.

    reg_t spectr = {{-1, -1, -1}, {-1, -1, -1}};									// All spectr.
    fft (&domain, data, &spectr);

    printf ("\n\nSpectr: %s\nModes with energy > 1e-12:\n", reg_printRanges (&spectr));
    MF_UNROLL (_s_, _o_, &spectr);											// Gets unrolling coefficients for fft'd data.
    for (int i = spectr.min[0] ; i <= spectr.max[0] ; ++i)
      for (int j = spectr.min[1] ; j <= spectr.max[1] ; ++j)
        for (int k = spectr.min[2] ; k <= spectr.max[2] ; ++k)
          if (MV_CELL(data, i, j, k, _s_, _o_) > 1e-12)
            printf ("mode(%d, %d, %d) -> %e\n", i, j, k, MV_CELL(data, i, j, k, _s_, _o_));
  }
}

#endif
