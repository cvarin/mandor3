/** \file shared.template.c
  * Template to generate 'dll.c' file used as 'distiller' backend.
  *
  * Generated C code is inserted here by Python frontend, the result is
  * compiled using make, and executed.
  *
  * \warning Please remember that template uses percent character '%%' to mark
  *          placeholders for embedded data. If you need percent in the final
  *          C code, write it here as double '%%%%' sign.
  *
  * Example:
  *    printf ("%%d %%le\n", i, x[i]);
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "type_marker.h"

/// Instead of comparing 'x = y' I compare '|x - y| < EPS*(|x| + |y|)', this
/// comparison is used to compare 'qDivM' parameters.
const double EPS = 1e-5;

static inline double
dot (double *a, double *b)
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// ---------------------------------------------------------------------------
/// Computes determinant of the 3x3 matrix.
// ---------------------------------------------------------------------------
static inline double
det (double a11, double a12, double a13,
     double a21, double a22, double a23,
     double a31, double a32, double a33)
{
   return a11*a22*a33 + a21*a32*a13 + a12*a23*a31 -
          a31*a22*a13 - a11*a23*a32 - a33*a21*a12;
}

// ---------------------------------------------------------------------------
/// Finds if particle is inside of box. To check, we write offset from the 'r0'
/// point as 'r - r0 = alpha*dr1 + beta*dr2 + gamma*dr3' and take dot products
/// with 'dr1', .., 'dr3' to derive equations to find 'alpha', 'beta', 'gamma'.
/// Point is inside of the box if all these coeffiecients are in [0,1].
// ---------------------------------------------------------------------------
static inline int
is_in_box (marker_t *p,
           double r0[3],
           double e1[3],
           double e2[3],
           double e3[3])
{
   double a11 = dot (e1, e1),
          a12 = dot (e1, e2),
          a13 = dot (e1, e3),
          a22 = dot (e2, e2),
          a23 = dot (e2, e3),
          a33 = dot (e3, e3);
   double e1dr = (p->x - r0[0])*e1[0] + (p->y - r0[1])*e1[1] + (p->z - r0[2])*e1[2],
          e2dr = (p->x - r0[0])*e2[0] + (p->y - r0[1])*e2[1] + (p->z - r0[2])*e2[2],
          e3dr = (p->x - r0[0])*e3[0] + (p->y - r0[1])*e3[1] + (p->z - r0[2])*e3[2];

   double sys = det (a11, a12, a13,
                     a12, a22, a23,
                     a13, a23, a33);
   double alpha_det = det (e1dr, a12, a13,
                           e2dr, a22, a23,
                           e3dr, a23, a33);
   if (alpha_det < 0 || alpha_det > sys)
      return 0;

   double beta_det = det (a11, e1dr, a13,
                          a12, e2dr, a23,
                          a13, e3dr, a33);
   if (beta_det < 0 || beta_det > sys)
        return 0;

   double gamma_det = det (a11, a12, e1dr,
                           a12, a22, e2dr,
                           a13, a23, e3dr);
   if (gamma_det < 0 || gamma_det > sys)
      return 0;

   return 1;
}

// ---------------------------------------------------------------------------
/// Tests if particle is inside of the sphere.
// ---------------------------------------------------------------------------
static inline int
is_in_sphere (marker_t *p, double r[3], double R)
{
   return (p->x - r[0])*(p->x - r[0]) +
          (p->y - r[1])*(p->y - r[1]) +
          (p->z - r[2])*(p->z - r[2]) <= R*R;
}

// ---------------------------------------------------------------------------
/// Tests if particle is inside of the cylinder. Cross product of displacement
/// vector with direction vector gives distance between point and axis. Squared
/// distances are compared to save square root computing.
// ---------------------------------------------------------------------------
static inline int
is_in_cylinder (marker_t *p, double r[3], double d[3], double R)
{
   double a[3] = {p->x - r[0], p->y - r[1], p->z - r[2]};
   double c[3] = {d[1]*a[2] - d[2]*a[1],
                  d[2]*a[0] - d[0]*a[2],
                  d[0]*a[1] - d[1]*a[0]};
   return dot (c, c) <= R*R*dot (d, d);
}

// ---------------------------------------------------------------------------
/// Tests if particle is under a plane specified by point and normal.
/// 'Under' means 'in half-space the normal vector doesn't point into'.
// ---------------------------------------------------------------------------
static inline int
is_under_plane (marker_t *p, double r[3], double n[3])
{
   double a[3] = {p->x - r[0], p->y - r[1], p->z - r[2]};
   return dot (a, n) <= 0;
}

// ---------------------------------------------------------------------------
/// Tests if particle is of given type.
// ---------------------------------------------------------------------------
static inline int
is_plasma (marker_t *p, double ref_qDivM)
{
   return fabs (p->qDivM - ref_qDivM) <
                                 EPS*0.5*(fabs (p->qDivM) + fabs (ref_qDivM));
}

// ---------------------------------------------------------------------------
/// Filters 'input' array to keep only particles which pass the test;
/// return value holds the number of particles left.
// ---------------------------------------------------------------------------
int
distiller_filter (marker_t *input, int N)
{
   // Loop-invariant constants (usually some geometry).
   %(invariants)s;

   int final = 0;
   for (int i = 0 ; i < N ; ++i) {
      marker_t *p = input + i;

      // Temporary variables (to shorten complex expressions).
      %(tmp_vars)s;

      if (%(final_test)s)
         input[final++] = input[i];
   }

   return final;
}
