/*
 * Mesh module --- subroutines to check limits, save, load, copy and clean a mesh-variable.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>
#include "type_mesh.h"

#include "log.h"


// Meshes for the test.
static mesh_t       test1, test2;
static meshDouble_p testDouble = mcast_meshDouble (&test1);
static meshVec_p    testVec = mcast_meshVec (&test2);

/*
 * Returns random number in [min, max].
 */
static int
mesh_randomInt (int min, int max)
{
  ENSURE (max >= min, "bad range: min (%d) > max (%d)", min, max);
  return min + rand ()%(max + 1 - min);
}

/*
 * Returns random number in [-1, 1].
 */
static double
mesh_randomDouble (void)
{
   return 2*(double)(rand ())/((double) RAND_MAX) - 1;
}

/*
 * This function tests that resize routine works properly.
 * NOTE: I NEED TO PASS THE SAME LIMITS TO ALL REALLOCATION CALLS TO DOUBLE AND VECTOR ARRAYS BECAUSE OF
 *       I USE ONE COMMON LOOP TO CHECK THE RESULT. OTHERWISE I BREAK BOUNDARY OF ARRAY.
 *
 */
void
mesh_testResize (void)
{
   int imin, jmin, kmin, imax, jmax, kmax;
   int i, j, k, di, dj, dk;

   // Arrays to save data in control points to compare.
   int I[5], J[5], K[5];
   vec3D_t saveVec[5][5][5];
   double saveDouble[5][5][5];

   // Sets initial size.
   imin = mesh_randomInt (5, 10)*mc_have_x;
   jmin = mesh_randomInt (5, 10)*mc_have_y;
   kmin = mesh_randomInt (5, 10)*mc_have_z;
   imax = (imin + mesh_randomInt (5, 10))*mc_have_x;
   jmax = (jmin + mesh_randomInt (5, 10))*mc_have_y;
   kmax = (kmin + mesh_randomInt (5, 10))*mc_have_z;

   // Allocates meshes.
   mesh_allocate (&test1, imin, jmin, kmin, imax, jmax, kmax, "test_Double", mc_double);
   mesh_allocate (&test2, imin, jmin, kmin, imax, jmax, kmax, "test_Vec", mc_vec3D_t);

   ENSURE (!mf_mesh_cmpSize (testDouble, testVec),  "bad initialization");

   // Initializes by random value.
   for (i = testVec->imin ; i <= testVec->imax ; i++)
   for (j = testVec->jmin ; j <= testVec->jmax ; j++)
   for (k = testVec->kmin ; k <= testVec->kmax ; k++) {
      mv_f (testDouble, i, j, k) = mesh_randomDouble ();
      mv_fx (testVec, i, j, k) = mesh_randomDouble ();
      mv_fy (testVec, i, j, k) = mesh_randomDouble ();
      mv_fz (testVec, i, j, k) = mesh_randomDouble ();
   }

   // Picks the control point.
   for (i = 0 ; i < 5 ; i++) {
      I[i] = mesh_randomInt (testVec->imin, testVec->imax);
      J[i] = mesh_randomInt (testVec->jmin, testVec->jmax);
      K[i] = mesh_randomInt (testVec->kmin, testVec->kmax);
   }

   say ("mesh_testResize: grid is\n  {%d, %d, %d, %d, %d}\n  {%d, %d, %d, %d, %d}\n  {%d, %d, %d, %d, %d}.",
        I[0], I[1], I[2], I[3], I[4], J[0], J[1], J[2], J[3], J[4], K[0], K[1], K[2], K[3], K[4]);

   // Saves the values.
   for (i = 0 ; i < 5 ; i++)
   for (j = 0 ; j < 5 ; j++)
   for (k = 0 ; k < 5 ; k++) {
      saveDouble[i][j][k] =  mv_f (testDouble, I[i], J[j], K[k]);
      saveVec[i][j][k] =  mv_v (testVec, I[i], J[j], K[k]);
   }

   // Expands the meshes.
   di = mesh_randomInt (0, 4)*mc_have_x;
   dj = mesh_randomInt (0, 4)*mc_have_y;
   dk = mesh_randomInt (0, 4)*mc_have_z;

   // Expands the array.
   mesh_reallocate (&test1, imin - di, jmin - dj, kmin - dk, imax + di*2, jmax + dj/2, kmax + dk);
   mesh_reallocate (&test2, imin - di, jmin - dj, kmin - dk, imax + di*2, jmax + dj/2, kmax + dk);

   // Verifies the values.
   for (i = 0 ; i < 5 ; i++)
   for (j = 0 ; j < 5 ; j++)
   for (k = 0 ; k < 5 ; k++) {
      ENSURE (saveDouble[i][j][k] ==  mv_f (testDouble, I[i], J[j], K[k]),
              "bad mv_f");
      ENSURE (saveVec[i][j][k].x ==  mv_fx (testVec, I[i], J[j], K[k]),
              "bad mv_fx");
      ENSURE (saveVec[i][j][k].y ==  mv_fy (testVec, I[i], J[j], K[k]),
              "bad mv_fy");
      ENSURE (saveVec[i][j][k].z ==  mv_fz (testVec, I[i], J[j], K[k]),
              "bad mv_fz");
   }

   say ("mesh_testResize: expand test is passed.");

   // Initializes by random value new array.
   for (i = testVec->imin ; i <= testVec->imax ; i++)
   for (j = testVec->jmin ; j <= testVec->jmax ; j++)
   for (k = testVec->kmin ; k <= testVec->kmax ; k++) {
      mv_f (testDouble, i, j, k) = mesh_randomDouble ();
      mv_fx (testVec, i, j, k) = mesh_randomDouble ();
      mv_fy (testVec, i, j, k) = mesh_randomDouble ();
      mv_fz (testVec, i, j, k) = mesh_randomDouble ();
   }

   // Saves the values.
   for (i = 0 ; i < 5 ; i++)
   for (j = 0 ; j < 5 ; j++)
   for (k = 0 ; k < 5 ; k++) {
      saveDouble[i][j][k] =  mv_f (testDouble, I[i], J[j], K[k]);
      saveVec[i][j][k] =  mv_v (testVec, I[i], J[j], K[k]);
   }

   // Shrinks the meshes back.
   mesh_reallocate (&test1, imin, jmin, kmin, imax, jmax, kmax);
   mesh_reallocate (&test2, imin, jmin, kmin, imax, jmax, kmax);

   // Verifies the values.
   for (i = 0 ; i < 5 ; i++)
   for (j = 0 ; j < 5 ; j++)
   for (k = 0 ; k < 5 ; k++) {
      ENSURE (saveDouble[i][j][k] ==  mv_f (testDouble, I[i], J[j], K[k]),
              "bad mv_f");
      ENSURE (saveVec[i][j][k].x ==  mv_fx (testVec, I[i], J[j], K[k]),
              "bad mv_fx");
      ENSURE (saveVec[i][j][k].y ==  mv_fy (testVec, I[i], J[j], K[k]),
              "bad mv_fy");
      ENSURE (saveVec[i][j][k].z ==  mv_fz (testVec, I[i], J[j], K[k]),
              "bad mv_fz");
   }

   say ("mesh_testResize: shrink test is passed.");

   mesh_free (&test1);
   mesh_free (&test2);
}
