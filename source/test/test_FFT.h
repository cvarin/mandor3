#ifndef test_FFT_header
#define test_FFT_header

void   test_FFTdomainWrap (int N, int k);
double test_FFT3D (int Nx, int Ny, int Nz, int tryN, double *dt);

double test_FFTW2 (int N, int tryN, double *execTime);

#endif
