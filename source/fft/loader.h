/** \file spectr_process.h
  * \brief Loads data saved using functions in spectr_dump.h and performs fourier transformation of the data saving
  * energy density for each field component. Its possible to reduce data size by setting high-frequency cut-off to save
  * space, especially for 3D runs, where size scales as third power of the node numbers. Cut-off can be used to remove
  * "dirty" region with \f$ \lambda < 5 h \f$, which cannot be resolved using mesh.
  *
  * <h3>Mapping note.</h3>
  * I use reg_t structure to hold region sizes. It lets me access elements of the unrolled arrays and compute unrolling
  * parameters using macroses mf_initUnroller() and mv_map(). I have \b two reg_t descriptor for the same computational
  * domain. One is for loading and holds all sizes the way they are defined in main computational module and the way
  * they are saved - that is \b spectr_process.c::domain. Another (spectr_process.c::fftMap) threats the same array as
  * unrolled 3D array with the same sizes along axises and all indices started from zero. This way is natural for fftw
  * and local processing.
  *
  * \attention Top boundary is \b NOT included here because of for periodical boundary conditions it is repeated and it
  * is not required for fft.
  *
  * \warning In function spectr_loadDump() mesh is initialized with \b NaN to detect early unloaded pieces. Coverage test
  * or some other type of consistency analysis may be much faster.
  */

#ifndef MC_SPECTR_PROCESS_HEADER
#define MC_SPECTR_PROCESS_HEADER						///< \internal Guard.

void spectr_convert (int recordNum, int NKx, int NKy, int NKz, const reg_t *subdomain, int removeDumps);

#endif
