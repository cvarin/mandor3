/** \file misc_ranmar.c
  * This random number generator originally appeared in "Toward a Universal
  * Random Number Generator" by George Marsaglia and Arif Zaman,
  * Florida State University Report: FSU-SCRI-87-50 (1987).
  *
  * It was later modified by F. James and published in "A Review of Pseudo
  * random Number Generators".
  *
  * THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
  *
  * (However, a newly discovered technique can yield a period of 10^600, but
  * that is still in the development stage.)
  *
  * It passes ALL of the tests for random number generators and has a period
  * of 2^144, is completely portable (gives bit identical results on all
  * machines with at least 24-bit mantissas in the floating point
  * representation).
  *
  * The algorithm is a combination of a Fibonacci sequence (with lags of 97
  * and 33, and operation "subtraction plus one, modulo one") and an
  * "arithmetic sequence" (using subtraction).
  *
  * On a Vax 11/780, this random number generator can produce a number in
  * 13 microseconds.
  */

#include <stdio.h>

#include "log.h"
#include "misc_ranmar.h"

static int    set = 0;
static double u[97], c, cd, cm;
static int    i97, j97;

void
ranmar_test (void)
{
   // Initialization (these are the seeds needed to produce the test case results).
   ranmar_init (1802, 9373);

   // Generates 20000 random numbers.
   double temp[100];
   for (int i = 0 ; i < 200 ; i++)
      ranmar (temp, 100);

   say ("If the random number generator works, the next six random\n"
        "numbers should be:\n"
        "  6533892.0  14220222.0   7275067.0\n"
        "  6172232.0   8354498.0  10633180.0\n"
   );

   ranmar (temp, 6);
   for (int i = 0 ; i < 6 ; i++)
      say ("%12.1f", 4096.0*4096.0*temp[i]);
}

/*
 * This is the initialization routine for the random number generator RANMAR()
 * NOTE: The seed variables can have values between:    0 <= IJ <= 31328
 *                                                      0 <= KL <= 30081
 * The random number sequences created by these two seeds are of sufficient
 * length to complete an entire calculation with. For example, if several
 * different groups are working on different parts of the same calculation,
 * each group could be assigned its own IJ seed. This would leave each group
 * with 30000 choices for the second seed. That is to say, this random
 * number generator can create 900 million different subsequences -- with
 * each subsequence having a length of approximately 10^30.
 *
 * Use IJ = 1802 & KL = 9373 to test the random number generator. The
 * subroutine RANMAR should be used to generate 20000 random numbers.
 * Then display the next six random numbers generated multiplied by 4096*4096
 * If the random number generator is working properly, the random numbers
 * should be:
 *           6533892.0  14220222.0  7275067.0
 *           6172232.0  8354498.0   10633180.0
 */
void
ranmar_init (int ij, int kl)
{
   if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081)
      DIE("the first random number seed must have a value between 0 and 31328, "
         "the second seed must have a value between 0 and 30081");

   int i = (ij/177)%177 + 2,
       j =  ij%177      + 2,
       k = (kl/169)%178 + 1,
       l =  kl%169;

   for (int ii = 0 ; ii < 97 ; ii++)
   {
      double sum         = 0,
             mantissaBit = 0.5;
      for (int bit = 0 ; bit < 24 ; bit++)
      {
         int m = (((i*j)%179)*k)%179;
         i = j;
         j = k;
         k = m;
         l = (53*l + 1)%169;
         if ((l*m)%64 >= 32)
            sum += mantissaBit;
         mantissaBit *= 0.5;
      }
      u[ii] = sum;
   }

   c  =   362436.0/16777216.0;
   cd =  7654321.0/16777216.0;
   cm = 16777213.0/16777216.0;

   i97 = 96;
   j97 = 32;

   set = 1;
}

/*
 * This is the random number generator proposed by George Marsaglia in
 * Florida State University Report: FSU-SCRI-87-50
 * It was slightly modified by F. James to produce an array of pseudorandom
 * numbers.
 */
void
ranmar (double *rvec, int len)
{
  int ivec;
  double uni;

  if (!set)
    error ("ranmar: random number generator is not initialized (use ranmar_init).");

  for (ivec = 0 ; ivec < len ; ivec++)
  {
    uni = u[i97] - u[j97];
    uni += (uni < 0.0);
    u[i97] = uni;
    i97 = (i97 - 1 + 97)%97;
    j97 = (j97 - 1 + 97)%97;
    c -= cd;
    c += cm*(c < 0.0);
    uni -= c;
    uni += (uni < 0.0);
    rvec[ivec] = uni;
  }
}
