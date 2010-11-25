/** \file reg_unroller.h
  * \brief Computes memory coverage for 2D/3D arrays unrolled into 1D array.
  *
  * This routine is core element of the modelling engine. Unrolling of the array into 1D
  * lump of memory is a way of becoming independent on the dimensionality at compile time.
  * Main idea is that frame.h will hold predefined macro-variables mc_have_x, mc_have_y
  * and mc_have_z which are set to '0' or '1'. Using this macro-variables in the addressing
  * macros will help to completely exclude unused dimensionality index at compile time.
  *
  * Because of different sets of axis may produce badly optimized loops special routine
  * can compute memory in 1D array we should process. Effectively it means replacement of
  * loops in form
  <pre>
  for (int i = reg.min[0] ; i <= reg.max[0] ; ++i)
    for (int j = reg.min[1] ; j <= reg.max[1] ; ++j)
      for (int k = reg.min[2] ; k <= reg.max[2] ; ++k)
        array[offset(i,j,k)] = ...;
  </pre>
  * by
  <pre>
  vec3i_t skip, N;
  long long int offset = reg_getMemCoverage (domain, reg, NULL, &skip, &N);
  for (int i = 0 ; i < N.v.x ; ++i, offset += skip.v.x)
    for (int j = 0 ; j < N.v.y ; ++j, offset += skip.v.y)
      for (int k = 0 ; k < N.v.z ; ++k, offset += skip.v.z)
        array[offset] = ...;
  </pre>
  *
  * One should note that if memory piece is continuous than efficiently only inner loop is not
  * trivial (have N > 1) and two outer loops do nothing.
  *
  * Other macroses are used to find element in 1D array (storage) which corresponds to
  * particular indices '(i, j, k)'. It is done by computing unrolling coefficients (strides)
  * and using them. Resulting strides are kept separate from region structure (one can use more
  * fields inside of reg_t) on purpose - if you do not have this arrays that means that strides
  * are not found yet and there will be no case you use uninitialized strides by accident.
  */

#ifndef MF_UNROLL

#include "frame.h"
#include "type_reg.h"
#include "type_vector.h"

/**
  * Sets strides to map cells withing 3D region into 1D unrolled array with (i, j, k) element
  * of 'reg' going to 'i*strides[0] + j*strides[1] + k*strides[2] + origin' ('origin' sets
  * unrolled address of 'reg->min[]' to the first cell of 1D array with address 0).
  */
#define MF_UNROLL(strides, origin, reg)												\
{																\
  (strides)[2] = mc_have_z;													\
  (strides)[1] = mc_have_y*(((reg)->max[2] - (reg)->min[2])*mc_have_z + 1);							\
  (strides)[0] = mc_have_x*(((reg)->max[1] - (reg)->min[1])*mc_have_y + 1)*(((reg)->max[2] - (reg)->min[2])*mc_have_z + 1);	\
  (origin) = - ((strides)[0]*(reg)->min[0] + (strides)[1]*(reg)->min[1] + (strides)[2]*(reg)->min[2]);				\
}

/**
  * Uses prepared unrolling parameters to get offset of '(i, j, k)' cell.
  */
#define MF_AIM(i, j, k, strides, origin)		((origin) + (i)*mc_have_x*(strides)[0] + (j)*mc_have_y*(strides)[1] + (k)*mc_have_z)

long long int reg_getMemCoverage (const reg_t *domain, const reg_t *reg, const vec3i_t *d, vec3i_t *skip, vec3i_t *N);

#endif
