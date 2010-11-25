/** \file tag_maxwell.c
  * Creates maxwellian distribution function.
  */

#include <math.h>
#include <stdlib.h>

#include "type_mesh.h"

#include "misc_PIC.h"
#include "misc_units.h"
#include "log.h"
#include "misc_cfgReader.h"

#include "setup/main.h"
#include "setup_denavit.h"
#include "setup/plasma.h"
#include "setup/tag_plasma.h"
#include "setup_distrMapper.h"

#define MC_VMAX			6.0	///< exp (-16) = 1.125e-7, exp (16) = 8.9e6.
#define MC_WING_WIDTH		1.6	///< Width of the sampled region for uniform spacing placement.
#define MC_MAPPER_SAMPLES_N	5000	///< Number of sample points in the DF mapper.

/// Parameters of maxwellian DF.
static double ux, uy, uz;
static int    tuneEnergy, uniformWeight;

static void maxwell_fillDF (const double rho, const double qDivM);

// ---------------------------------------------------------------------------
/// Entry point.
// ---------------------------------------------------------------------------
double
tag_maxwell (FILE *fp)
{
   double omegaPeSquared = cfg_readDouble (fp);
   double qDivM          = cfg_readDouble (fp);
   ux                    = cfg_readDouble (fp);
   uy                    = cfg_readDouble (fp);
   uz                    = cfg_readDouble (fp);
   tuneEnergy            = cfg_readInt (fp);
   uniformWeight         = 1 == cfg_readInt (fp);

   // Turns frequency into charge density.
   double rho = omegaPeSquared*4*mc_pi*mc_pi/qDivM/(double) plasma_PPC;
   ENSURE (qDivM*rho > 0, "Inconsistent q/M(%.3e) and charge density (%.3e).", qDivM, rho);

   const double speedToTemperature = mc_CGS_c*mc_CGS_c*mc_CGS_m/fabs (qDivM*mc_CGS_eV);
   say ("tag_maxwell:");
   say ("  - forcing distribution function to Maxwellian");
   say ("  - rho = %f (n/n_cr = %f), q/M = %e", rho*plasma_PPC, rho*qDivM*plasma_PPC/(4*mc_pi*mc_pi), qDivM);
   say ("  - ux = %f (Tx = %.3e)", ux, speedToTemperature*(sqrt (ux*ux + 1) - 1));
   say ("  - uy = %f (Ty = %.3e)", uy, speedToTemperature*(sqrt (uy*uy + 1) - 1));
   say ("  - uz = %f (Tz = %.3e)", uz, speedToTemperature*(sqrt (uz*uz + 1) - 1));

   if (tuneEnergy)
      say ("  - thermal energy is adjusted explicitly");

   if (!memEstimateOnly) {
      maxwell_fillDF (rho, qDivM);
      say ("  - %d particles in a pattern (uniform %s)", plasma_PPC, (uniformWeight) ? "weight" : "spacing");
   }

   return (((cpu_max[0] - cpu_min[0])*mc_have_x + 1.0)*
           ((cpu_max[1] - cpu_min[1])*mc_have_y + 1.0)*
           ((cpu_max[2] - cpu_min[2])*mc_have_z + 1.0) + 1.0)*
           plasma_PPC*sizeof (marker_t);
}

// ---------------------------------------------------------------------------
/// Maxwell DF for mapper generator (see setup_distrMapper.c).
// ---------------------------------------------------------------------------
static double
func_Maxwell (double x)
{
   return exp (- x*x);
}

// ---------------------------------------------------------------------------
/// Creates pattern, aplies symmetrizator and copies it into all cells.
// ---------------------------------------------------------------------------
static void
maxwell_fillDF (const double rho, const double qDivM)
{
   // Initializes distribution function mappers.
   mapper_t *mapperAll      = mapper_create (func_Maxwell, - MC_VMAX, + MC_VMAX, MC_MAPPER_SAMPLES_N);
   mapper_t *mapperPositive = mapper_create (func_Maxwell, 0, + MC_VMAX, MC_MAPPER_SAMPLES_N);

   ENSURE (plasma_layers % 2 == 0, "For mirror balancing we need odd number of layers.");

   // Initializes the reference distribution (velocity distribution for markers in the single cell).
   marker_t pattern[plasma_PPC];
   for (int i = 0, p = 0 ; i < plasma_nx ; ++i)					// Fills reference pattern.
      for (int j = 0 ; j < plasma_ny ; ++j)
         for (int k = 0 ; k < plasma_nz ; ++k)
            for (int l = 0 ; l < plasma_layers/2 ; ++l) {
               double vx, vy, vz, tmp;
               denavit_createQuartet (p/2, plasma_PPC/2, &tmp, &vx, &vy, &vz);

               if (uniformWeight) {
                  vx             = mapper_invoke (vx, mapperPositive);
                  vy             = mapper_invoke (vy, mapperAll);
                  vz             = mapper_invoke (vz, mapperAll);
                  pattern[p].rho = 1.0;
               } else {
                  vx            *= MC_WING_WIDTH;
                  vy             = MC_WING_WIDTH*2*(vy - 0.5);
                  vz             = MC_WING_WIDTH*2*(vz - 0.5);
                  pattern[p].rho = func_Maxwell (vx)*func_Maxwell (vy)*func_Maxwell (vz);
               }

               pattern[p].vx = vx*ux;						// Allocates marker in one velocity hemispere.
               pattern[p].vy = vy*uy;
               pattern[p].vz = vz*uz;
               ++p;

               pattern[p].vx  = - vx*ux;					// Allocates marker in opposite hemisphere.
               pattern[p].vy  = - vy*uy;
               pattern[p].vz  = - vz*uz;
               pattern[p].rho = pattern[p-1].rho;
               ++p;
            }

   mapper_destroy (mapperAll);							// Deallocates mappers.
   mapper_destroy (mapperPositive);

   double norm = 0;								// Normalizes the total charge density.
   for (int p = 0 ; p < plasma_PPC ; ++p)
      norm += pattern[p].rho;
   norm = plasma_PPC*rho/norm;

   for (int p = 0 ; p < plasma_PPC ; ++p)
      pattern[p].rho *= norm;

   if (tuneEnergy) {								// Normalizes total thermal energy.
      double Wx = 0, Wy = 0, Wz = 0;
      for (int p = 0 ; p < plasma_PPC ; ++p) {
         Wx += pattern[p].vx*pattern[p].vx*pattern[p].rho;
         Wy += pattern[p].vy*pattern[p].vy*pattern[p].rho;
         Wz += pattern[p].vz*pattern[p].vz*pattern[p].rho;
      }

      ENSURE (Wx != 0 && Wy != 0 && Wz != 0, "Bad Tx(%.3e), Ty(%.3e), or Tz(%.3e) temperature.", Wx, Wy, Wz);

      Wx = ux*sqrt (rho*plasma_PPC*0.5/Wx);
      Wy = uy*sqrt (rho*plasma_PPC*0.5/Wy);
      Wz = uz*sqrt (rho*plasma_PPC*0.5/Wz);
      for (int p = 0 ; p < plasma_PPC ; ++p) {
         pattern[p].vx *= Wx;
         pattern[p].vy *= Wy;
         pattern[p].vz *= Wz;
      }

      if (!uniformWeight)
         say ("  o weight in nonuniform - make sure wings are wide enough or energy adjusting will screw DF");
   }

   long int N;
   for (marker_t *p = plasma_getObject (0, &N), *end = p + N; p < end ; ++p) {	// Forces velocity distribution.
      int l = MF_DBL_TO_INT (p->vx + 0.1);
      p->vx = pattern[l].vx;
      p->vy = pattern[l].vy;
      p->vz = pattern[l].vz;
      p->rho = pattern[l].rho;
      p->qDivM = qDivM;
   }
}
