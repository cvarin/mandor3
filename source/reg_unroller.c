/** \file reg_unroller.c
  * \brief Computes memory coverage for 2D/3D arrays unrolled into 1D array.
  *
  * Please see reg_unroller.h for details.
  */

#include <assert.h>

#include "reg_unroller.h"

// ---------------------------------------------------------------------------
/// \brief Takes region/subregion and returns memory coverage information (see reg_unroll.h header).
/// Input parameter 'd' is optional steps for loops (NULL means normal steps of size '1').
// ---------------------------------------------------------------------------
long long int
reg_getMemCoverage (const reg_t *domain, const reg_t *reg, const vec3i_t *d, vec3i_t *skip, vec3i_t *N)
{
  assert (reg_isInside (reg, domain));

  static const vec3i_t unit = {{1, 1, 1}};										// To avoid temporary for unit steps.
  d = (d) ? d : &unit;
  assert (d->v.x > 0 && d->v.y > 0 && d->v.z > 0);

  long long int _s_[3], _o_;
  MF_UNROLL(_s_, _o_, domain);												// Sets unrolling strides and origin.

  N->v.x = ((reg->max[0] - reg->min[0])/d->v.x)*mc_have_x + 1;								// Nx
  N->v.y = ((reg->max[1] - reg->min[1])/d->v.y)*mc_have_y + 1;								// Ny
  N->v.z = ((reg->max[2] - reg->min[2])/d->v.z)*mc_have_z + 1;								// Nz

  skip->v.x = MF_AIM(d->v.x, reg->min[1], 0, _s_, _o_) - MF_AIM(0, reg->min[1] + N->v.y*d->v.y, 0, _s_, _o_);		// Skip to continue new X line.
  skip->v.y = MF_AIM(0, d->v.y, reg->min[2], _s_, _o_) - MF_AIM(0, 0, reg->min[2] + N->v.z*d->v.z, _s_, _o_);		// Skip to continue new Y line.
  skip->v.z = MF_AIM(0, 0, d->v.z, _s_, _o_) - MF_AIM(0, 0, 0, _s_, _o_);						// Skip for single step along Z.

  for (int iter = 0 ; iter < 2 ; ++iter)										// Bubble-sorting-type variant.
    for (int l = 2 ; l > 0 ; --l)
    {
      if (skip->r[l-1] == 0)												// No additional skip => merge loops.
      {
        N->r[l] *= N->r[l-1];
        N->r[l-1] = 1;
      }
      if (N->r[l] == 1 && skip->r[l] == 0)										// Removes trivial iteration.
      {
        N->r[l] = N->r[l-1];
        N->r[l-1] = 1;
        skip->r[l] = skip->r[l-1];
        skip->r[l-1] = 0;
      }
    }

  return MF_AIM(reg->min[0], reg->min[1], reg->min[2], _s_, _o_);					 		// Offset of the beginning (first span).
}
