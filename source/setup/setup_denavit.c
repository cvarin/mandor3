/** \file setup_denavit.c
  * Quiet start sampler of a phase space as described by Denavit and Walsh.
  *
  * Fills unit cube in 4D space using method of Hammersley from paper
  *
  *   J. Denavit, J.M. Walsh, 'Nonrandom initialization of particle codes',
  *   Comments Plasma Phys. Cont. Fusion, v.6, 1981.
  *
  * Than I map the unit intervals on proper distributions (see 'setup_distMapper.h').
  */

#include <math.h>
#include <stdio.h>

#include "log.h"

// XXX - refactor file

// ---------------------------------------------------------------------------
/// Finds representation of 'num' under given 'radix' so 'num = a0*radix^0 + a1*radix^1 + a2*radix^2 + ...'
/// and returns bit-reversed number 'res = a0/radix^0 + a1/radix^1 + a2/radix^2 + ...'.
// ---------------------------------------------------------------------------
static double
denavit_bitReverse (int num, int radix)
{
   double res         = 0,
          invRadix    = 1.0/(double) radix,
          digitWeight = invRadix;
   while (num > 0) {
      int digit = num % radix;
      num /= radix;
      res += digit*digitWeight;
      digitWeight *= invRadix;
   }
   return res;
}

// ---------------------------------------------------------------------------
/// Returns sample number 'i' in length 'N' sequence.
// ---------------------------------------------------------------------------
void
denavit_createQuartet (int i, int N, double *reg, double *phi2, double *phi3, double *phi5)
{
   ENSURE (i >= 0 && i < N, "sequence number is out of range: i = %d, N = %d",
                            i, N);

   // Uniformly generates 4 'uncorrelated' randoms in [0,1] range.
   *reg = (i + 0.5)/(double) N;
   *phi2 = denavit_bitReverse (i + 1, 2);
   *phi3 = denavit_bitReverse (i + 1, 3);
   *phi5 = denavit_bitReverse (i + 1, 5);
}

// ---------------------------------------------------------------------------
/// This function exports data in data-file (to compare pictures from here with Denavit's paper).
// ---------------------------------------------------------------------------
void
denavit_test (void)
{
   int   N;
   printf ("Testing of the uniform cube coverint (sampling).\nInput number of points:");
   scanf ("%d", &N);

   FILE *fp = fopen ("output/denavit.dat", "wt");
   fprintf (fp, "variables = N, x, y, z, t\nzone t=\"N=%d\"", N);
   for (int i = 0 ; i < N ; ++i)
   {
      double x, y, z, t;
      denavit_createQuartet (i, N, &x, &y, &z, &t);
      fprintf (fp, "%d %e %e %e %e\n", i, x, y, z, t);
   }

   fclose (fp);
}
