/** \file fft.h
  * FFT transformation of the 3D region using temporary 1D array. In-place,
  * returns spectral energy density \f$|f_{\vec k}|^2\f$, where:
  * \f[
  *   f_{\vec r} = \sum\limits_{\vec k} f_{\vec k}\cdot \exp(i\vec k\cdot\vec r)
  * \f]
  * Possible to compute and return only subregion of the spectral domain:
  * usefull to ignore high-frequency regions. Core job is done by
  * <a href="http://www.fftw.org">fftw3</a>.
  *
  * \attention Trainer is commented at the end of fft.c - use it to test
  *            library independently.
  *
  * \attention Top boundary is \b NOT included or passed here because of for
  *            periodical boundary conditions it is excessive information which
  *            is not required by fft.
  *
  * \attention FFT transforms f[i], i = 0, .., N-1 into f[k], k = 0, .., N-1.
  *            It is lossless and f[-k] is aliased into f[N-k]. In my case I
  *            want to use the same array fft.c::fft_source and spectr is
  *            assumed to be in [-Nk, Nk] region mapped by reg_t *spectr. That
  *            means that for odd N we loose one mode with the highest wave
  *            number (not critical for visualization - this mode is highly
  *            numerical).
  *
  * \todo Input signal is real, so for energy of the modes one can cut off half
  *       of the spectral domain, e.g.
  *       \f$[-N_{Kx},N_{Kx}]\times[-N_{Ky},N_{Ky}]\times[-N_{Kz},N_{Kz}]\f$ ->
  *       \f$[-N_{Kx},N_{Kx}]\times[0,N_{Ky}]\times[-N_{Kz},N_{Kz}]\f$
  *       (something like that, any axis).
  */

#ifndef MC_FFT_HEADER
#define MC_FFT_HEADER		///< \internal Guard.

void fft (const reg_t *domain, double *data, reg_t *spectr);

#endif
