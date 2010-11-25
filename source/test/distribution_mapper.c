/** \file distribution_mapper.c
  * \brief Test of the distribution mapper (see setup_distrMapper.c, setup_distrMapper.h).
  *
  * Creates mappers for simple functions \f$ f(x) = x \f$ and \f$ f(x) = exp (- x^2) \f$ and checks accuracy, etc. Results are shown in
  * <em>output/test_mapper.dat</em>, <em>output/test_mapper2.dat</em> files accessible through <em>lay/test_mapper*.lay</em> layouts.
  */

#include <math.h>
#include <stdlib.h>

#include "log.h"
#include "misc_cfgReader.h"

#include "setup/setup_distrMapper.h"

// ---------------------------------------------------------------------------
/// Function to test accuracy of mapper.
// ---------------------------------------------------------------------------
static double
func_lin (double x)
{
  return x;
}

// ---------------------------------------------------------------------------
/// Simple tester of the mapper (uses linear distribution function to generate mapper).
// ---------------------------------------------------------------------------
void
mapper_testLinear (void)
{
  int i;
  FILE *fp;
  double y, dy;
  mapper_t *mapper;

  say ("Testing of the distribution function mapper using linear function with known inverse mapping.");
  say ("See file output/test_mapper.dat for results.");

  mapper = mapper_create (func_lin, 0, sqrt (2.0), 20);

  fp = cfg_open ("output/test_mapper.dat", "wt", __func__);
  fprintf (fp, "variables = y invGood invNum delta\nzone t=\"test\", f = point\n");
  for (y = 0, dy = 0.01 ; y <= 1.0 ; y += dy)
  {
    double x1 = mapper_invoke (y, mapper);
    fprintf (fp, "%e %e %e %e\n", y, sqrt (2*y), x1, x1 - sqrt (2*y));
  }

  fprintf (fp, "zone t=\"integral\", f = point\n");
  for (i = 0 ; i < mapper->N ; i++)
    fprintf (fp, "%e %e %e %e\n", mapper->arg[i], mapper->func[i], 0.5*mapper->arg[i]*mapper->arg[i], mapper->func[i] - 0.5*mapper->arg[i]*mapper->arg[i]);
  fclose (fp);

  mapper_destroy (mapper);
}

// ---------------------------------------------------------------------------
/// Function to test accuracy of mapper.
// ---------------------------------------------------------------------------
static double
func_maxwell (double x)
{
  return exp (-x*x);
}

// ---------------------------------------------------------------------------
/// More advanced tester of the mapper (uses maxwellian distribution function to generate mapper).
// ---------------------------------------------------------------------------
void
mapper_testMaxwell (void)
{
  int i;
  FILE *fp;
  double y, dy;
  mapper_t *mapper;

  say ("Testing of the distribution function mapper using Maxwell function with known inverse mapping.");
  say ("See file output/test_mapper2.dat for results.");

  mapper = mapper_create (func_maxwell, - 10.0, 10.0, 5000);								// Creates mapper.

  const double norm = 1.0/(erf (mapper->arg[mapper->N]) - erf (mapper->arg[0]));					// Normalization for analitic result.
  fp = cfg_open ("output/test_mapper2.dat", "wt", __func__);
  fprintf (fp, "variables = x yGood yNum delta\nzone t=\"test\", f = point\n");
  for (y = 0.0, dy = 0.01 ; y <= 1.0 + 1e-7 ; y += dy)
  {
    double x = mapper_invoke (y, mapper);
    double analitic = norm*(erf(x) - erf (mapper->arg[0]));
    fprintf (fp, "%e %e %e %e\n", x, y, analitic, y - analitic);
  }

  fprintf (fp, "zone t=\"integral\", f = point\n");
  for (i = 0 ; i <= mapper->N ; i++)
  {
    double analitic = norm*(erf (mapper->arg[i]) - erf (mapper->arg[0]));
    fprintf (fp, "%e %e %e %e\n", mapper->arg[i], analitic, mapper->func[i], mapper->func[i] - analitic);
  }
  fclose (fp);

  mapper_destroy (mapper);
}
