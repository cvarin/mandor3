/** \file em.c
  * Maxwell solver:
  *   + Yee, central differences, staggered mesh
  *   + second order of approximation in space and time
  *   + conserves \f$ div \vec E \f$ and \f$ div \vec H \f$ with numerical
  *     accuracy.
  *
  * Mesh layout is documented in 'em.h'.
  *
  * Full time step is done as few substeps which are traced by 'emStage'
  * variable to force a proper order of calls, which is
  * -# \f$ H^{n-1/2} \to H^n \f$ by 'em_HHalfStep'
  * -# \f$ H^n \to H^{n+1/2} \f$ by 'em_HStep'
  * -# \f$ E^n \to E^{n+1}   \f$ by 'em_EStep_start' and 'em_EStep_finish'
  *
  * \warning Current density source is added after parallel syncronization of
  *          the current density.
  */

#include <stdlib.h>
#include <assert.h>

#include "log.h"
#include "profiler.h"

#include "em.h"
#include "em_caps.h"
#include "em_TFSF.h"
#include "em_gaussSpot.h"

#include "scpic/em_sources.h"

/// Flag to check order of funcion calls (see top of the page).
static int emStage = 0;

// ---------------------------------------------------------------------------
/// Initializes em-subsystem.
// ---------------------------------------------------------------------------
void
em_init (meshVec_p E, meshVec_p H)
{
   ENSURE (1.0/(tau*tau) > 1.0/(h1*h1) + 1.0/(h2*h2) + 1.0/(h3*h3),
           "Yee-solver is unstable");

   cap_init       (E, H);
   gaussSpot_init ();
   TFSF_init      ();
   init_sources   ();
}

// ---------------------------------------------------------------------------
/// Half-step in time for magnetic field in inner region.
// ---------------------------------------------------------------------------
static void
em_HSmallStep (meshVec_RO_p E, meshVec_p H)
{
   const int i1 = (cpu_min[0] + 3)*mc_have_x,
             j1 = (cpu_min[1] + 3)*mc_have_y,
             k1 = (cpu_min[2] + 3)*mc_have_z;
   const int i2 = (cpu_max[0] - 2)*mc_have_x,
             j2 = (cpu_max[1] - 2)*mc_have_y,
             k2 = (cpu_max[2] - 2)*mc_have_z;
   double    c1 = 0.5*tau/h1*mc_have_x,
             c2 = 0.5*tau/h2*mc_have_y,
             c3 = 0.5*tau/h3*mc_have_z;

   for (int i = i1 ; i <= i2 ; i++)
   for (int j = j1 ; j <= j2 ; j++)
   for (int k = k1 ; k <= k2 ; k++) {
      mv_fx(H, i, j, k) += c3*(mv_fy(E, i, j, k) - mv_fy(E, i,   j,   k-1))
                         - c2*(mv_fz(E, i, j, k) - mv_fz(E, i,   j-1, k));
      mv_fy(H, i, j, k) += c1*(mv_fz(E, i, j, k) - mv_fz(E, i-1, j,   k))
                         - c3*(mv_fx(E, i, j, k) - mv_fx(E, i,   j,   k-1));
      mv_fz(H, i, j, k) += c2*(mv_fx(E, i, j, k) - mv_fx(E, i,   j-1, k))
                         - c1*(mv_fy(E, i, j, k) - mv_fy(E, i-1, j,   k));
   }
}

// ---------------------------------------------------------------------------
/// Initial half step for magnetic field: caches new fields in the caps, starts
/// parallel exchange, makes step in inner region, and updates field covered by
/// caps.
// ---------------------------------------------------------------------------
void
em_HHalfStep (meshVec_RO_p E, meshVec_p H)
{
   // Checks if previous step is finished.
   assert (emStage == 0);

   // Records field on TFSF-interface.
   profiler_begin (mc_prof_HHalf_TFSF_preH);
   TFSF_preHStep (mcast_meshVecI(E), mcast_meshVecI(H));

   // Caches and advances boundary fields.
   profiler_endBegin (mc_prof_HHalf_cacheEH);
   cap_cacheEH (E, H);

   // Half-step for PIC solver.
   profiler_endBegin (mc_prof_HHalf);
   em_HSmallStep (E, H);

   // Updates boundary layers for PIC solver.
   profiler_endBegin (mc_prof_HHalf_flushH);
   cap_flushHalfH (H);

   // Updates flags.
   profiler_end ();
   ++emStage;
}

// ---------------------------------------------------------------------------
/// Final half step for magnetic field.
// ---------------------------------------------------------------------------
void
em_HStep (meshVec_RO_p E, meshVec_p H)
{
   assert (emStage == 1);	// Checks if previous step is finished.

   em_HSmallStep (E, H);	// Advances field in inner subdomain.
   cap_flushH    (H);		// Flushes new field from for boundary cap.

   // Maintains interface discontinuity.
   const reg_t innerRegH = {{(cpu_min[0] + 3)*mc_have_x, (cpu_min[1] + 3)*mc_have_y, (cpu_min[2] + 3)*mc_have_z},
                            {(cpu_max[0] - 2)*mc_have_x, (cpu_max[1] - 2)*mc_have_y, (cpu_max[2] - 2)*mc_have_z}};
   TFSF_postHStep (mcast_meshVecI(H),       &innerRegH);
   add_H_sources  (H,                 Time, &innerRegH);

   ++emStage;
}

// ---------------------------------------------------------------------------
/// Time step for electric field in inner region.
// ---------------------------------------------------------------------------
void
em_EStep_start (meshVec_p E, meshVec_RO_p H)
{
   assert (emStage == 2);	// Checks if previous step is finished.

   const int i1 = (cpu_min[0] + 2)*mc_have_x,
             j1 = (cpu_min[1] + 2)*mc_have_y,
             k1 = (cpu_min[2] + 2)*mc_have_z;
   const int i2 = (cpu_max[0] - 2)*mc_have_x,
             j2 = (cpu_max[1] - 2)*mc_have_y,
             k2 = (cpu_max[2] - 2)*mc_have_z;
   double    c1 = tau/h1,
             c2 = tau/h2,
             c3 = tau/h3;

   for (int i = i1 ; i <= i2 ; i++)
   for (int j = j1 ; j <= j2 ; j++)
   for (int k = k1 ; k <= k2 ; k++) {
      mv_fx(E, i, j, k) -= c3*(mv_fy(H, i,   j,   k+1) - mv_fy(H, i, j, k))
                         - c2*(mv_fz(H, i,   j+1, k  ) - mv_fz(H, i, j, k));
      mv_fy(E, i, j, k) -= c1*(mv_fz(H, i+1, j,   k  ) - mv_fz(H, i, j, k))
                         - c3*(mv_fx(H, i,   j,   k+1) - mv_fx(H, i, j, k));
      mv_fz(E, i, j, k) -= c2*(mv_fx(H, i,   j+1, k  ) - mv_fx(H, i, j, k))
                         - c1*(mv_fy(H, i+1, j,   k  ) - mv_fy(H, i, j, k));
   }

   ++emStage;			// Steps to next stage.
}

// ---------------------------------------------------------------------------
/// Time step for electric field - accounts current density term in separate pass.
// ---------------------------------------------------------------------------
static void
em_EStep_addCurrents (meshVec_p E, meshVec_RO_p J)
{
   const int i1 = (cpu_min[0] - 1)*mc_have_x,
             j1 = (cpu_min[1] - 1)*mc_have_y,
             k1 = (cpu_min[2] - 1)*mc_have_z;
   const int i2 = (cpu_max[0] + 1)*mc_have_x,
             j2 = (cpu_max[1] + 1)*mc_have_y,
             k2 = (cpu_max[2] + 1)*mc_have_z;

   for (int i = i1 ; i <= i2 ; i++)
   for (int j = j1 ; j <= j2 ; j++)
   for (int k = k1 ; k <= k2 ; k++) {
      mv_fx(E, i, j, k) -= tau*mv_fx(J, i, j, k);
      mv_fy(E, i, j, k) -= tau*mv_fy(J, i, j, k);
      mv_fz(E, i, j, k) -= tau*mv_fz(J, i, j, k);
   }
}

// ---------------------------------------------------------------------------
/// Time step for electric field - completes time step by updating fields under caps, adding current density contribution and so on.
// ---------------------------------------------------------------------------
void
em_EStep_finish (meshVec_p E, meshVec_p H, meshVec_RO_p J)
{
   assert (emStage == 3);

   // Updates boundary layer.
   profiler_begin (mc_prof_plasma_ebc_cap);
   cap_flushE (E);

   // Adds to E an interface input.
   reg_t innerRegE = {{(cpu_min[0] + 2)*mc_have_x, (cpu_min[1] + 2)*mc_have_y, (cpu_min[2] + 2)*mc_have_z},
                      {(cpu_max[0] - 2)*mc_have_x, (cpu_max[1] - 2)*mc_have_y, (cpu_max[2] - 2)*mc_have_z}};
   TFSF_postEStep (mcast_meshVecI(E),       &innerRegE);

   // Gets parallel messages.
   profiler_endBegin (mc_prof_plasma_ebc_recv);
   cap_catchParallelMsgs (E, H);

   reg_t domain = {{cpu_min[0] - 1, cpu_min[1] - 1, cpu_min[2] - 1},
                   {cpu_max[0] + 1, cpu_max[1] + 1, cpu_max[2] + 1}};
   add_E_sources (E, Time, &domain);

   // Adds contribution of the currents.
   profiler_endBegin (mc_prof_plasma_ebc_jPass);
   em_EStep_addCurrents (E, J);

   // Forces hard gauss-spot source, adds contribution of soft sources.
   profiler_endBegin (mc_prof_plasma_ebc_misc);
   gaussSpot_forceE (E, Time + tau);

   TFSF_completeFrame ();		// Marks completed frame.
   emStage = 0;				// Updates flag.
   profiler_end ();
}

// ---------------------------------------------------------------------------
/// Calculates total energy of EM field in local cpu domain.
// ---------------------------------------------------------------------------
void
em_energy (meshVec_RO_p E, meshVec_RO_p H, double *W_E, double *W_M)
{
   const int i1 = cpu_min[0],
             j1 = cpu_min[1],
             k1 = cpu_min[2];
   const int i2 = cpu_max[0] + (1 - mc_have_x),
             j2 = cpu_max[1] + (1 - mc_have_y),
             k2 = cpu_max[2] + (1 - mc_have_z);

   double we = 0,
          wb = 0;
   for (int i = i1 ; i < i2 ; ++i)
   for (int j = j1 ; j < j2 ; ++j)
   for (int k = k1 ; k < k2 ; ++k) {
      double ex = mv_fx(E, i+1, j, k),
             ey = mv_fy(E, i, j+1, k),
             ez = mv_fz(E, i, j, k+1),
             hx = mv_fx(H, i, j+1, k+1),
             hy = mv_fy(H, i+1, j, k+1),
             hz = mv_fz(H, i+1, j+1, k);
      we += ex*ex + ey*ey + ez*ez;
      wb += hx*hx + hy*hy + hz*hz;
   }
   *W_E = 0.5*we*h1*h2*h3;
   *W_M = 0.5*wb*h1*h2*h3;
}
