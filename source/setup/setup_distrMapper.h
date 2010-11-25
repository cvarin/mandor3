/** \file setup_distrMapper.h
  * Headers and data structures for mapper (see setup_distrMapper.c).
  *
  * The job of mapper is to receive a distribution function f(x), to build integral
  * distribution function "F(x) = integral (f(t)*dt, [-oo, x])". After that for each
  * Y in [0, 1] mapper can find corresponding X so F(X) = Y as shown below:
  * <pre>
  *   ^ F(x)
  * 1 +                     ----------
  *   |                   /
  * Y +------------------+
  *   |                 /|
  * 0 +     -----------  |
  *   |                  |
  *   +------------------+---------------------> x
  *                      X
  * </pre>
  *
  * Function F(x) is stored in const-step-array and integrated numerically.
  * Search for X is done using binary division search and smooth interpolation.
  */

#ifndef MC_MISC_FUNCMAPPER_HEADER
#define MC_MISC_FUNCMAPPER_HEADER

/**
  * Structured data for integral mapping of one function to another (see setup_distrMapper.c, setup_distrMapper.h).
  */
typedef struct
{
  int     N;		///< Number of points in data array (inclusive, 0 <= i <= N).
  double *func;         ///< \f$ \displaystyle \int_{-\infty}^{x_i} func(x')dx' \f$.
  double *arg;		///< \f$ x_i \f$.
} mapper_t;

mapper_t *mapper_create (double (*func)(double x), double xMin, double xMax, int N);
double    mapper_invoke (double y, const mapper_t *mapper);
void      mapper_destroy (mapper_t *mapper);

#endif
