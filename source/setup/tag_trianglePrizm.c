/** \file tag_trianglePrizm.c
  * This functions creates triangle prizm made of cold uniform plasma.
  */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <mpi.h>

#include "type_marker.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_partition.h"
#include "misc_parameters.h"

#include "main.h"

static double x1, y1, x2, y2, x3, y3, z1, z2;					///< Parameters of the prizm.

// ---------------------------------------------------------------------------
/// Tests if the (x,y,z) point is inside of the prizm.
// ---------------------------------------------------------------------------
// XXX static
int
pointIsInsidePrizm (double x, double y, double z)
{
  #define MF_VEC_PRODUCT(UX, UY, VX, VY) ((UX)*(VY) - (UY)*(VX))
  double chirality1 = MF_VEC_PRODUCT(x1 - x, y1 - y, x2 - x, y2 - y);
  double chirality2 = MF_VEC_PRODUCT(x2 - x, y2 - y, x3 - x, y3 - y);
  double chirality3 = MF_VEC_PRODUCT(x3 - x, y3 - y, x1 - x, y1 - y);
  #undef MF_VEC_PRODUCT

  return (chirality1*chirality2 >= 0) && (chirality1*chirality3 >= 0) && (chirality2*chirality3 >= 0) &&
         mc_have_z*(z - z1)*(z2 - z) >= 0;
}

// One-liners to keep from messing up with defines or from multiple evaluations of an arguments.
static inline int 	min (int a, int b) 	{  	return (a < b) ? a : b; 	}
static inline int 	max (int a, int b)	{	return (a > b) ? a : b;		}

// ---------------------------------------------------------------------------
/// Creates triangle prizm made of a cold uniform plasma.
// ---------------------------------------------------------------------------
double
tag_trianglePrizm (FILE *fp)
{
#if 0
  x1 = cfg_readDouble (fp)*units (mc_micron);					// First vertex (all coordinates are in microns).
  y1 = cfg_readDouble (fp)*units (mc_micron);
  x2 = cfg_readDouble (fp)*units (mc_micron);					// Second vertex.
  y2 = cfg_readDouble (fp)*units (mc_micron);
  x3 = cfg_readDouble (fp)*units (mc_micron);					// Third vertex.
  y3 = cfg_readDouble (fp)*units (mc_micron);
  z1 = cfg_readDouble (fp)*units (mc_micron);					// Z coordinate of the lower end of the prizm.
  z2 = cfg_readDouble (fp)*units (mc_micron);					// Z coordinate of the upper end of the prizm.

  double rho = cfg_readDouble (fp);						// Gets charge density.

  if (z2 < z1) 									// Makes sure top and bottom are ok.
  {
    double tmp = z1;
    z1 = z2, z2 = tmp;
  }

  rho *= (rho < 0) ? (- units (mc_ne_critical)) : 1;				// Expresses concentrations in [cm^-3].

  const double q = cfg_readDouble (fp);						// Gets charge in [e].
  const double M = cfg_readDouble (fp);						// Gets mass in [m].

  int nx = cfg_readInt (fp);							// Gets particles allocation density.
  int ny = cfg_readInt (fp);
  int nz = cfg_readInt (fp);
  nx = mc_have_x*(nx - 1) + 1;							// Prepares allocation density for 1D/2D.
  ny = mc_have_y*(ny - 1) + 1;
  nz = mc_have_z*(nz - 1) + 1;
  assert (nx > 0 && ny > 0 && nz > 0 &&
      "Bad marker allocation density (positive 'nx', 'ny', and 'nz' expected).");

  say ("tag_trianglePrizm: ");						// Prints report about final configuration.
  say ("  o q = %.2f [|e|], M = %.2f [m_e]", q, M);
  say ("  o n = %f [n_e critical]", rho/units (mc_ne_critical));
  say ("  o n = %.3e [cm^-3]", rho);
  say ("  o vertex #1: (%f, %f) [micron],", x1/units (mc_micron), y1/units (mc_micron));
  say ("  o vertex #2: (%f, %f) [micron],", x2/units (mc_micron), y2/units (mc_micron));
  say ("  o vertex #3: (%f, %f) [micron],", x3/units (mc_micron), y3/units (mc_micron));
  say ("  o z1 = %f [micron], z2 = %f [micron]", z1/units (mc_micron), z2/units (mc_micron));
  say ("  o sampling density %d %d %d particles per cell", nx, ny, nz);

  int I1 = min (x1/h1, min (x2/h1, x3/h1)) - 1;					// Computes integer prizm bounding box;
  int J1 = min (y1/h2, min (y2/h2, y3/h2)) - 1;					//   additional shifts are to account for
  int K1 = min (z1/h3, z2) - 1; 		 				//   negative x1, .., z2.
  int I2 = max (x1/h1, max (x2/h1, x3/h1)) + 2;
  int J2 = max (y1/h2, max (y2/h2, y3/h2)) + 2;
  int K2 = max (z1/h3, z2/h3) + 2;

  int i1 = mc_have_x*nx*max (cpu_min[0], I1);					// Converts all indices to markers steps.
  int j1 = mc_have_y*ny*max (cpu_min[1], J1);
  int k1 = mc_have_z*nz*max (cpu_min[2], K1);
  int i2 = mc_have_x*(nx*min (cpu_max[0], I2) - 1);
  int j2 = mc_have_y*(ny*min (cpu_max[1], J2) - 1);
  int k2 = mc_have_z*(nz*min (cpu_max[2], K2) - 1);

  rho *= q*mc_CGS_e/(units (mc_rho0)*nx*ny*nz);					// Turns concentration into marker charge density.

  if(!memEstimateOnly)
  {
    cache_startChapter ();							// Starts new chapter (block of marker pages).
    cache_setSpecie (q/M);
  }

  double part = 0, dx = h1/nx, dy = h2/ny, dz = h3/nz;				// Number of particles and steps between them.
  for (int i = i1 ; i <= i2 ; ++i)						// Allocates target.
  {
    double x = mc_have_x*(i + 0.5)*dx;
    for (int j = j1 ; j <= j2 ; ++j)
    {
      double y = mc_have_y*(j + 0.5)*dy;
      for (int k = k1 ; k <= k2 ; ++k)
      {
        double z = mc_have_z*(k + 0.5)*dz;
        if (pointIsInsidePrizm (x, y, z))
        {
          if (!memEstimateOnly)
          {
            marker_t *marker = cache_getParticle ();			// Gets new cached marker pointer.
            marker->x = x;
            marker->y = y;
            marker->z = z;
            marker->vx = 0.0;
            marker->vy = 0.0;
            marker->vz = 0.0;
            marker->rho = rho;
          }
          ++part;
        }
      }
    }
  }

  cache_flush ();							// Flushes particles to main storage.

  double sum;								// Gets total number of particles on master.
  MPI_Reduce (&part, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (cpu_here)
    return part*sizeof (marker_t);

  // Master checks that foil is well resolved.
  if (sum < nx*ny*nz*max (I2 - I1 + 1, max (J2 - J1 + 1, K2 - K1 + 1)))
    SAY_WARNING ("\nToo small number of particles - prizm seems to be undersampled.\n");
  say ("  o %.3f particles is on cpu %d.", part, cpu_here);

  return part*sizeof (marker_t);
#endif
  assert (0);
  return 0;
}
