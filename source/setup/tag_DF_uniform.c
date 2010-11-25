/** \file tag_DF_uniform.c
  * Forces uniform distribution function is allocated plasma (see
  * 'tag_DF_uniform.h').
  */

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "setup/plasma.h"
#include "setup/tag_plasma.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_partition.h"
#include "misc_parameters.h"

#include "main.h"

// ---------------------------------------------------------------------------
/// Returns uniformly distributed random number from [0, 1) (note: drand48 is
/// told to be obsolete so I don't use it).
// ---------------------------------------------------------------------------
static inline double
frand (void)
{
   return (double) rand ()/(double) RAND_MAX;
}

// ---------------------------------------------------------------------------
/// Sets uniform velocity distribution function and uniform charge/ionization.
// ---------------------------------------------------------------------------
void
tag_DF_uniform (FILE *fp)
{
   double rho_crit = mc_CGS_e*units(mc_ne_critical)/units(mc_rho0);

   double rho   = cfg_readDouble (fp),
          qDivM = cfg_readDouble (fp),
          ux    = cfg_readDouble (fp),	// Mean velocity.
          uy    = cfg_readDouble (fp),
          uz    = cfg_readDouble (fp),
          dux   = cfg_readDouble (fp),	// Velocity spread.
          duy   = cfg_readDouble (fp),
          duz   = cfg_readDouble (fp);

   ENSURE (rho*qDivM > 0,
           "negative mass; sign of ρ(%.3e) ≠ sign of 'q/M'(%.3e)", rho, qDivM);

   // Normalizes rho.
   rho *= rho_crit/plasma_PPC;

   // Applies velocity distribution.
   long int N;
   for (marker_t *p = plasma_getObject (0, &N), *end = p + N; p < end ; ++p) {
      p->vx    = ux + (frand () - 0.5)*dux;
      p->vy    = uy + (frand () - 0.5)*duy;
      p->vz    = uz + (frand () - 0.5)*duz;
      p->rho   = rho;
      p->qDivM = qDivM;
   }

   say ("%s: uniform velocity distribution is set", __func__);
   say ("  - rho = %.3e[rho_crit], q/M = %.3e", rho*plasma_PPC/rho_crit,
                                                qDivM);
   say ("  - <v> = (%e, %e, %e)", ux,  uy,  uz);
   say ("  - dv  = (%e, %e, %e)", dux, duy, duz);
}
