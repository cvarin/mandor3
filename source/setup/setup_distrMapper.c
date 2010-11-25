/** \file setup_distrMapper.c
  * A group of functions to work with integral mapping of one function to another.
  *
  * Typical usage is to map uniformly distributed random variable to one with given
  * distribution function for quiet start. Structures are described in setup_distrMapper.h.
  *
  * Main notes: input parameter is distribution function \f$ f(x) \f$ which is
  * transformed into integral distribution function
  *   \f$ \displaystyle F(x) = \frac 1\alpha \cdot \int\limits_{-\infty}^{x} f(\zeta)\,d\zeta \f$.
  * Normalisation coefficient \f$ \alpha \f$ is easy to find by taking integral and
  * putting \f$ \alpha = F(+\infty) \f$.
  *
  * Integral is taken by approximating function using 4 points and two quadratic
  * splines (and adding linear interpolation to blend left and right splines). At the
  * boundary only one of splines is used. Points are uniformly spaces so expressions
  * are quite simplified. All integrals are written for x > 0, so for x < 0 sign goes
  * to the 'b' coefficient(s) of \f$ y(x) = a\cdot (x - x_1)^2 + b\cdot (x - x_1) + c\f$.
  *
  * Interpolation of the random point is done using quadratic spline again but this
  * time points are not uniformly spaced. At the boundary interpolation falls back
  * to linear interpolation.
  */

#include <math.h>
#include <stdlib.h>

#include "log.h"

#include "setup_distrMapper.h"

#define MC_EPSILON	(1e-7)		///< Accuracy for range check inside mapper_invoke ().

// ---------------------------------------------------------------------------
/// Allocates new mapper structure.
// ---------------------------------------------------------------------------
static mapper_t *
mapper_alloc (int N)
{
  mapper_t *mapper = (mapper_t*) malloc (sizeof (mapper_t));

  ENSURE (mapper, "Cannot allocate %d bytes.", sizeof (mapper_t));		// Allocates memory.

  mapper->N = N;
  mapper->func = (double*) malloc (sizeof (double)*(N + 1));
  mapper->arg  = (double*) malloc (sizeof (double)*(N + 1));

  ENSURE (mapper->func && mapper->arg, "Cannot allocate memory to store func data (%f Kb).",
             2*(N + 1)*sizeof (double)/1024.0);

  return mapper;
}

// ---------------------------------------------------------------------------
/// Defines 'a', 'b', 'c' that \f$ y(x) = a\cdot (x - x_1)^2 + b\cdot (x - x_1) + c \f$
/// and \f$ y(x_i) = y_i \f$. Spatial step is constant and equal to 'h'.
// ---------------------------------------------------------------------------
#define mf_quadSpline(h, y0, y1, y2, a, b, c)										\
{															\
  (a) = ((y2) - 2*(y1) + (y0))/(2*(h)*(h));										\
  (b) = ((y2) - (y0))/(2*(h));												\
  (c) = y1;														\
}

// ---------------------------------------------------------------------------
/// Evaluates \f$ \displaystyle \int\limits_0^h (a\cdot x^2 + b\cdot x + c) \cdot dx \f$.
// ---------------------------------------------------------------------------
#define mf_integralAll(h, a, b, c) ((a)*(h)*(h)*(h)/3.0 + 0.5*(b)*(h)*(h) + (c)*(h))

// ---------------------------------------------------------------------------
/// Evaluates \f$ \displaystyle \int\limits_0^h (a\cdot x^2 + b\cdot x + c)\cdot \frac xh \cdot dx \f$.
// ---------------------------------------------------------------------------
#define mf_integralRight(h, a, b, c) ((a)*(h)*(h)*(h)/4.0 + (b)*(h)*(h)/3.0 + 0.5*(c)*(h))

// ---------------------------------------------------------------------------
/// Evaluates \f$ \displaystyle \int\limits_0^h (a\cdot x^2 + b\cdot x + c)\cdot \left(1 - \frac xh\right) \cdot dx \f$.
// ---------------------------------------------------------------------------
#define mf_integralLeft(h, a, b, c) (mf_integralAll(h, a, b, c) - mf_integralRight(h, a, b, c))

// ---------------------------------------------------------------------------
/// Builds integral mapping data for function \b func.
// ---------------------------------------------------------------------------
mapper_t *
mapper_create (double (*func)(double x), double xMin, double xMax, int N)
{
  ENSURE (N > 5 && xMin < xMax - 1e-10, 					// 4 points required for spline.
             "Bad parameters; N(%d), xMin (%e), xMax (%e).", N, xMin, xMax);

  mapper_t *mapper = mapper_alloc (N);						// Links mapper and allocates memory.

  const double dx = (xMax - xMin)/(double) N;

  double x = xMin + dx;
  double y0 = func (x - dx), y1 = func (x), y2 = func (x + dx), sum;
  double a, b, c;								// Left quad spline coefficients.

  mapper->arg[0] = xMin;							// Left boundary data: start point.
  mapper->func[0] = 0;
  mf_quadSpline (dx, y0, y1, y2, a, b, c);					// Defines quad spline.
  mapper->arg[1] = x;								// Left boundary data: quad-spline integral.
  sum = mapper->func[1] = mf_integralAll(dx, a, -b, c);				// '-b' due to negative sign of displacement.

  x += dx;
  for (int i = 2 ; i < N ; ++i)							// All inner points.
  {
    // Inside loop x corresponds to the upper limit of integral (to y2 point).
    double y3 = func (x + dx), A, B, C;						// Right quad spline coefficients.
    mf_quadSpline (dx, y1, y2, y3, A, B, C);					// Defines right quad spline.

    sum += mf_integralLeft (dx, a, b, c) + mf_integralLeft (dx, A, -B, C);	// Makes integration from x1 to x = x2.

    mapper->func[i] = sum;							// Saves data point.
    mapper->arg[i] = x;

    x += dx;									// Shifts stencil to the right by one step.
    y0 = y1;
    y1 = y2;
    y2 = y3;
    a = A;
    b = B;
    c = C;
  }

  sum += mf_integralAll (dx, a, b, c);						// Right quad spline integral.
  mapper->func[N] = sum;
  mapper->arg[N] = x;

  sum = 1.0/mapper->func[N];
  for (int i = 1 ; i <= N ; i++)						// Normalisation.
    mapper->func[i] *= sum;

  return mapper;
}

// ---------------------------------------------------------------------------
/// Returns \f$ x: F(x) = y \f$ for given mapper. Root is found by binary tree search and smooth interpolation.
// ---------------------------------------------------------------------------
double
mapper_invoke (double y, const mapper_t *mapper)
{
  int iL = 0, iR = mapper->N;

  ENSURE (y >= 0.0 - MC_EPSILON && y <= 1.0 + MC_EPSILON, "Bad mapping ordinate Y (%e).\n", y);

  while (iR - iL > 1)								// Binary divisions search.
  {
    int i = (iL + iR) >> 1;
    if (y < mapper->func[i])
      iR = i;
    else
      iL = i;
  }

  if (iL == 0 || iL == mapper->N - 1)						// Linear fall-back for an edge points.
  {
    double delta = (y - mapper->func[iL])/(mapper->func[iL+1] - mapper->func[iL]);
    return mapper->arg[iL] + delta*(mapper->arg[iL+1] - mapper->arg[iL]);
  }

  double x0 = mapper->func[iL-1], x1 = mapper->func[iL], x2 = mapper->func[iL+1];
  double y0 = mapper->arg[iL-1], y1 = mapper->arg[iL], y2 = mapper->arg[iL+1];
  double a = ((y2 - y1)/(x2 - x1) - (y0 - y1)/(x0 - x1))/(x2 - x0);
  double b = (y2 - y1)/(x2 - x1) - a*(x2 - x1);
  double dy = y - mapper->func[iL], left = a*dy*dy + b*dy + mapper->arg[iL];	// Left quadratic spline guess.

  x0 = mapper->func[iL];
  x1 = mapper->func[iL+1];
  x2 = mapper->func[iL+2];
  y0 = mapper->arg[iL];
  y1 = mapper->arg[iL+1];
  y2 = mapper->arg[iL+2];
  a = ((y2 - y1)/(x2 - x1) - (y0 - y1)/(x0 - x1))/(x2 - x0);
  b = (y2 - y1)/(x2 - x1) - a*(x2 - x1);
  dy = y - mapper->func[iL+1];
  double right = a*dy*dy + b*dy + mapper->arg[iL+1];				// Right quadratic spline guess.

  // Final interpolation.
  return - dy/(mapper->func[iL+1] - mapper->func[iL])*left + (1.0 + dy/(mapper->func[iL+1] - mapper->func[iL]))*right;
}

// ---------------------------------------------------------------------------
/// Destroys mapper to free allocated memory.
// ---------------------------------------------------------------------------
void
mapper_destroy (mapper_t *mapper)
{
  if (mapper)
  {
    free (mapper->arg);
    free (mapper->func);
    free (mapper);
  }
}
