/** \file realToInt.c
  * \brief Simple test of the real to integer conversion (basic operation at particles inputs evaluation).
  */

#include <math.h>

#include <mpi.h>

#include "real2int.c"

#include "log.h"

#ifdef test_main
  ERROR: ONLY ONE TEST MAY BE PLUGGED AT THE TIME AND test_main is defined already!
#endif

#define test_main() test_real2int()		///< Replaces entry point by ours.
#define mc_counts (1000)			///< Number of convertions to perform.

// ---------------------------------------------------------------------------
/// Performs 100000 convertions.
// ---------------------------------------------------------------------------
static void
test_real2int (void)
{
  double startTime = MPI_Wtime (), delta = startTime - ceil (startTime), dSum = 1.193845719345;
  long int sum = 0;

  for (int l = 0 ; l < 10000 ; ++l)
  {
    for (int i = 0 ; i < mc_counts ; ++i)
    {
      sum += double2int (dSum);
//       sum += (int) (dSum + 1000) - 1000;
//       sum += floor (dSum);
      dSum += delta;
      delta *= 0.997;
    }
    dSum = dSum - ceil (dSum);
    sum = sum % 10;
  }

  say ("%d000 iterations have taken %f sec (%d/%f).", mc_counts, MPI_Wtime () - startTime, sum, dSum);
}
