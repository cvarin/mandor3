/** \file tag_foil.c
  * Alternative sampling of the foil (see tag_foil.h).
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "type_marker.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_partition.h"
#include "misc_parameters.h"

#include "main.h"
#include "setup/plasma.h"

/// Small margin to prevent overlap of two layers due to round-off errors.
#define MC_MULTILAYER_MARGIN	(1e-10)

#define MF_MIN(a, b) (((a) < (b)) ? (a) : (b))	///< Minimum of 2 numbers.
#define MF_MAX(a, b) (((a) > (b)) ? (a) : (b))	///< Maximum of 2 numbers.

enum {PROF_LINEAR, PROF_EXP};

// ---------------------------------------------------------------------------
/// Estimates number of particles in the foil ('N' is a number of particles
/// per cell).
// ---------------------------------------------------------------------------
static double
tag_foil_estimate (double x,  double y,  double z,
                   double lx, double ly, double lz,
                   double width, int N)
{
   double sum = 0;
   for (int i = dmn_min[0] ; i <= dmn_max[0] ; ++i)
   for (int j = dmn_min[1] ; j <= dmn_max[1] ; ++j)
   for (int k = dmn_min[2] ; k <= dmn_max[2] ; ++k) {
      double distance = (i*h1 - x)*lx + (j*h2 - y)*ly + (k*h3 - z)*lz;
      if (fabs (distance) < 0.5*width)  sum += N;
   }
   return sum;
}

// ---------------------------------------------------------------------------
/// Creates foil with given profile of the denisty. Returns memory used to
/// store particles [bytes].
// ---------------------------------------------------------------------------
double
tag_foil (FILE *fp)
{
   // Reads position of the foil's center.
   double x0 = cfg_readDouble (fp)*units (mc_micron)*mc_have_x,
          y0 = cfg_readDouble (fp)*units (mc_micron)*mc_have_y,
          z0 = cfg_readDouble (fp)*units (mc_micron)*mc_have_z;

   // Reads direction of face normal to the foil.
   double lx = cfg_readDouble (fp)*mc_have_x,
          ly = cfg_readDouble (fp)*mc_have_y,
          lz = cfg_readDouble (fp)*mc_have_z;
   double l = sqrt (lx*lx + ly*ly + lz*lz);
   ENSURE (l > 1.0e-2, "too small foil normal vector: |(%e, %e, %e)| = %e",
                       lx, ly, lz, l);
   lx /= l;
   ly /= l;
   lz /= l;

   // Reads width and optional shift of the foil.
   double width = cfg_readDouble (fp)*units (mc_micron);
   double shift = 0;
   if (cfg_isOption (fp)) {
      shift = cfg_readOptDouble (fp)*units (mc_micron);
   }

   // Gets type of the profile.
   const char *profileType = cfg_readWord (fp);
   int profile = cfg_identifyWord (profileType,
                                   "linear",      PROF_LINEAR,
                                   "exponential", PROF_EXP,
                                   mc_cfgTermGuesses);
   ENSURE (profile != -1, "bad profile '%s' (use 'linear' or 'exponential')",
                          profileType);

   // Makes copy of the profile type string.
   char name[30];
   strncpy (name, profileType, 30);
   name[29] = 0;

   // Gets start/end concentration and charge to mass ratio.
   double n1 = cfg_readDouble (fp),
          n2 = cfg_readDouble (fp),
          q  = cfg_readDouble (fp),
          M  = cfg_readDouble (fp);

   // Converts concentrations in [cm^-3].
   n1 *= (n1 < 0) ? (- units (mc_ne_critical)) : 1;
   n2 *= (n2 < 0) ? (- units (mc_ne_critical)) : 1;

   // Gets number of particles per cell per axis.
   int nx = cfg_readInt (fp),
       ny = cfg_readInt (fp),
       nz = cfg_readInt (fp);

   // Checks width.
   ENSURE (width >= (mc_have_x*h1 + mc_have_y*h2 + mc_have_z*h3)
                  / (mc_have_x + mc_have_y + mc_have_z),
           "unresolved foil (width = %.3e), steps are (%.3e, %.3e, %.3e)",
           width, mc_have_x*h1, mc_have_y*h2, mc_have_z*h3);

   ENSURE (nx > 0 && ny > 0 && nz > 0,
           "bad particles density; nx (%d), ny(%d), nz(%d)",
           nx, ny, nz);

   say ("%s: ", __func__);
   say ("  - %s profile,", name);
   say ("  - q = %.2f [|e|], M = %.2f [m_e]", q, M);
   say ("  - n_start = %.3f [n_e critical], n_end = %.3f [n_e critical],",
        n1/units (mc_ne_critical),
        n2/units (mc_ne_critical));
   say ("  - n_start = %.3e [cm^-3], n_end = %.3e [cm^-3],", n1, n2);

   // Turns concentration to charge density.
   n1 *= q*mc_CGS_e/units (mc_rho0);
   n2 *= q*mc_CGS_e/units (mc_rho0);

   // Accounts for deactivated axises.
   nx = mc_have_x*(nx - 1) + 1;
   ny = mc_have_y*(ny - 1) + 1;
   nz = mc_have_z*(nz - 1) + 1;

   // Prints the rest of parameters.
   say ("  - center @ (%f, %f, %f) [micron]", x0/units (mc_micron),
                                              y0/units (mc_micron),
                                              z0/units (mc_micron));
   say ("  - width = %f [micron], shift = %f [micron]", width/units(mc_micron),
                                                        shift/units(mc_micron));
   say ("  - particles placement: %d x %d x %d", nx, ny, nz);

   if (memEstimateOnly) {
      return sizeof (marker_t)*tag_foil_estimate (x0, y0, z0,
                                                  lx, ly, lz,
                                                  width, nx*ny*nz);
   }

   // Allocates new plasma object.
   plasma_newObject ();

   // Gets spacing between particles.
   const double dx = h1*mc_have_x/(double) nx,
                dy = h2*mc_have_y/(double) ny,
                dz = h3*mc_have_z/(double) nz;

   const double margin = MC_MULTILAYER_MARGIN*		\
         (mc_have_x*h1 + mc_have_y*h2 + mc_have_z*h3) /	\
         ((mc_have_x + mc_have_y + mc_have_z)*width);

   // Number of particles.
   double  part    = 0;

   // Margins between foil and domain boundary.
   const int xMargin = 5*mc_have_x*(dmn_bc_min[0] != BC_PERIODIC);
   const int yMargin = 5*mc_have_y*(dmn_bc_min[1] != BC_PERIODIC);
   const int zMargin = 5*mc_have_z*(dmn_bc_min[2] != BC_PERIODIC);

   // Cells to allocate on local node.
   const int I1 = mc_have_x*MF_MAX(dmn_min[0] + xMargin,
                                   cpu_min[0]),
             J1 = mc_have_y*MF_MAX(dmn_min[1] + yMargin,
                                   cpu_min[1]),
             K1 = mc_have_z*MF_MAX(dmn_min[2] + zMargin,
                                   cpu_min[2]),
             I2 = mc_have_x*MF_MIN(dmn_max[0] - xMargin - 1,
                                   cpu_max[0] - 1),
             J2 = mc_have_y*MF_MIN(dmn_max[1] - yMargin - 1,
                                   cpu_max[1] - 1),
             K2 = mc_have_z*MF_MIN(dmn_max[2] - zMargin - 1,
                                   cpu_max[2] - 1);

   // Allocates particles.
   for (int I = I1 ; I <= I2 ; ++I)
   for (int J = J1 ; J <= J2 ; ++J)
   for (int K = K1 ; K <= K2 ; ++K)
   for (int i = 0  ; i < nx  ; ++i)
   for (int j = 0  ; j < ny  ; ++j)
   for (int k = 0  ; k < nz  ; ++k) {
      const double x = I*h1 + (i + 0.5)*dx;
      const double y = J*h2 + (j + 0.5)*dy;
      const double z = K*h3 + (k + 0.5)*dz;

      double alpha = shift - ((x - x0)*lx + (y - y0)*ly + (z - z0)*lz);

      // '0.5' is a shift to place the center into the ref-point.
      alpha = alpha/width + 0.5;
      if (alpha < margin || alpha > 1.0 - margin) continue;

      marker_t *marker = plasma_marker ();
      marker->x = x;
      marker->y = y;
      marker->z = z;
      marker->vx = 0;
      marker->vy = 0;
      marker->vz = 0;
      marker->qDivM = q/M;

      // Saves interpolation coefficient instead of charge density.
      marker->rho = alpha;
      ++part;
   }

   say ("  - %.3e particles added @ cpu %d.", part, cpu_here);

   // 'Charge density -> marker density' factor.
   const double weight = 1.0/((double) nx*ny*nz);

   // Gets foil as linear array.
   long int  N;
   marker_t *p   = plasma_getObject (0, &N),
            *end = p + N;
   switch (profile) {
      case PROF_LINEAR:
         for ( ; p < end ; ++p) p->rho = weight*(n1 + (n2 - n1)*p->rho);
         break;

      case PROF_EXP:
         for ( ; p < end ; ++p) p->rho = weight*n1*pow (n2/n1, p->rho);
         break;
   }

   return part*sizeof (marker_t);
}
