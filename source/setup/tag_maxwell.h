/** \file tag_maxwell.h
  * Creates three-themperature maxwellian distribution function.
  *
  * Creates maxwellian DF using "quiet start" from Denavit paper. Basic
  * pattern is copied several times after applying some symmetrization.
  * 'ux', 'uy', 'uz' are thermal speeds:
  *
  *     f(vx, vy, vz) ~ exp (- vx^2/ux^2 - .. - vz^2/uz^2)).
  *
  * Example of the config file entry is
    <pre>
    [Maxwell]
    @ 5.0       n/n_cr
    @ -1.0      q/m
    @ 0.01      ux
    @ 0.01      uy
    @ 0.01      uz
    @ 1         tuneEnergy flag
    @ 1         uniformWeight
    </pre>
  *
  * \note Quiet start uses bit reversal generation of sampled coordinate and for X, Y and Z
  * radix of the bit reversal procedure is different => for small N quality of sampling is
  * different for different directions (probably have to sample using spherical velocity
  * coordinates).
  *
  * \todo Energy normalization - define v_min, v_max, take integral
  * \f$ \displaystyle \int\limits_{v_{min}}^{v_{max}} f(v)\ dv \f$ and
  * normalize on it rather on value \f$ \displaystyle \int\limits_{-\infty}^{+\infty} f(v)\ dv \f$.
  *
  * \todo Check debay scale resolution (to avoid self-heating runs).
  *
  * \todo Use units like in tag_photoelectrons ().
  *
  * \todo \b Rotational transformations symmetry question is tricky - we need it to ensure
  * symmetry in appications which already have some cylindrical symmetry in mind. In this
  * case it may be necessary to disregard statistical independency of all transversal components
  * and initialize them not like \f$ f(v_x)\cdot f(v_y)\cdot f(v_z) \f$ but rather like
  * \f$ f(v)\cdot f(\theta)\cdot f(\phi) \f$.
  */

#ifndef MC_TAG_MAXWELL_HEADER
#define MC_TAG_MAXWELL_HEADER

#include <stdio.h>

double tag_maxwell (FILE *fp);

#endif
