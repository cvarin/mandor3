/** \file plasma.c
  * Particles solvers (Boris solver, current density evaluation, charge density
  * evaluation, interfaces to read kinetic energy). See details about currents
  * at comments to "plasmaMove_currentBoundConditions" routine.
  *
  * NOTES: if Courant condition holds, than
  *            jx[-1][?][?]  = jy[?][-1][?]  = jz[?][?][-1]  = 0
  *            rho[-1][?][?] = rho[?][-1][?] = rho[?][?][-1] = 0
  *        explicitly.
  *
  * \TODO Consider removing plasma_multistep for 1D/2D - never will happen and
  * so this big piece of code may be removed safely at compile time (to make
  * optimizer happy).
  */

#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "type_marker.h"

#include "misc_PIC.h"
#include "log.h"
#include "misc_cfgReader.h"
#include "misc_markerPlacer.h"

#include "timer.h"
#include "plasma.h"
#include "profiler.h"
#include "parr_meshes.h"

#include "core/plasma_rho.h"
#include "core/plasma_VSP.h"
#include "core/plasma_parallel.h"

#define MC_BUFFER 	(10<<20)        ///< Buffer for incoming/outgoing particles (see plasma_init ()).
#define MC_CHECK_FREQ	(2000)		///< Number of markers we advance between testing parallel exchange in plasma_timestep ().

// See 'plasma.h' for detailed legend.
marker_t *plasma     = NULL;	///< Storage.
long int  countAll   = 0, 	///< Storage size (in markers).
          countCore  = 0, 	///< Position of the first particle of the core.
          countShell = 0; 	///< Number of particles in shell.

/// Boundaries of the core and of the shell.
static double shellX1, shellX2, shellY1, shellY2, shellZ1, shellZ2;
static double coreX1,  coreX2,  coreY1,  coreY2,  coreZ1,  coreZ2;

/// Kinetic energies of the particles.
static double       plasmaWx = NAN,
                    plasmaWy = NAN,
                    plasmaWz = NAN;

/// Connection for current density.
static connection_t connectionJ = mf_connection_init (TAG_PLASMA_J, "plasma:J");
static meshVec_p    plasma_J;

/// Simple 'enum' to remove magic numbers in calls to 'plasma_timestep'.
enum stage_e {DOING_CORE = 0, DOING_SHELL = 1};

#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
static FILE *tracer = NULL;
#endif

#include "plasma_walls.c"
#include "plasma_currentWalls.c"

/*
 *     oooo o                                                      o
 *   M"    "M  oo   ooo    oo  oooo  oo  oooo   ooooo   oo oooo   oMoooo    oooo o
 *  M           M     M     MM"       MM"     o"     "o  M"    M   M       M    "o
 *  M           M     M     M         M       M""""""""  M     M   M        "ooooo
 *   Mo     o"  Mo   oM     M         M       "oo    oo  M     M   M    oo Mo   oM
 *     """""     """" ""   """"      """"       """""   """" """"   """"   " """"
 */

/*
 * Here I import plasmaCurrents_kernel - routine to evaluate current density
 * using initial and final positions.
 */
#include "plasma_currentKernel.c"

/*                                       M
 *   ""M"""""o
 *     M    oM    o"""""o   ""Moo"""   ""M      o""""o"
 *     M"""""o   M       M    M"         M      "o    "
 *     M      M  M       M    M          M      o """"M
 *   ooMoooooM"   "ooooo"   ooMooo    oooMooo   M"oooM"
 */

// ---------------------------------------------------------------------------
/// Plasma time step (external iterator 'p' is set by client and used by routers ).
// ---------------------------------------------------------------------------
static void
plasma_timestep (enum stage_e doing_shell, long int N,
                 meshVec_RO_p E, meshVec_RO_p H, meshVec_p J)
{
   profiler_endBegin (mc_prof_plasma_mainLoop);
   marker_t *p = plasma + (doing_shell ? 0 : countAll - N);
   for ( ; N > 0 ; ) {
      marker_t *shell = plasma + countShell,
               *core  = plasma + countCore;
      // Performs timestep for 'page' particles between parallel tests.
      long int page = (N < MC_CHECK_FREQ || doing_shell) ? N : MC_CHECK_FREQ;
      N -= page;
      for ( ; page > 0 ; --page, ++p) {
#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
         // TODO: "buffer with given size + flushes"?
         //       proc: smaller number of calls to 'fwrite'
         //       proc: streamed write
         //       cons: additional logic (about 7 lines of code)
         if (p->id < 0)
	    fwrite (p, sizeof(marker_t), 1, tracer);
#endif
         // '_' corresponds to indices, shifted by '-1/2'.
         int i, j, k, i_, j_, k_;
         int i2, j2, k2, di, dj, dk;
         const double x1 = p->x, y1 = p->y, z1 = p->z;
         double sigmaX, sigmaY, sigmaZ, sigmaX_, sigmaY_, sigmaZ_;

         MF_PIC_SIGMAS (x1, h1, i, sigmaX, i_, sigmaX_);
         MF_PIC_SIGMAS (y1, h2, j, sigmaY, j_, sigmaY_);
         MF_PIC_SIGMAS (z1, h3, k, sigmaZ, k_, sigmaZ_);

         const double Ex = mf_PIC_readX (E, i_, sigmaX_, j,  sigmaY,  k,  sigmaZ),
                      Ey = mf_PIC_readY (E, i,  sigmaX,  j_, sigmaY_, k,  sigmaZ),
                      Ez = mf_PIC_readZ (E, i,  sigmaX,  j,  sigmaY,  k_, sigmaZ_),
                      Hx = mf_PIC_readX (H, i,  sigmaX,  j_, sigmaY_, k_, sigmaZ_),
                      Hy = mf_PIC_readY (H, i_, sigmaX_, j,  sigmaY,  k_, sigmaZ_),
                      Hz = mf_PIC_readZ (H, i_, sigmaX_, j_, sigmaY_, k,  sigmaZ);

         // Half acceleration by electric field.
         const double tauMass = 0.5*tau*p->qDivM,
                      px      = p->vx + tauMass*Ex,
                      py      = p->vy + tauMass*Ey,
                      pz      = p->vz + tauMass*Ez,
                      gamma   = sqrt (1 + px*px + py*py + pz*pz);

         // Updates energy at t = n*tau + guard against p^2 = 0 (nans in diag.dat).
         const double W = (gamma - 1.0)*p->rho/(px*px + py*py + pz*pz + 1e-19)*h1*h2*h3/p->qDivM;
         plasmaWx += W*px*px;
         plasmaWy += W*py*py;
         plasmaWz += W*pz*pz;

         const double tauGamma = tauMass/gamma;
         const double aHx = tauGamma*Hx,
                      aHy = tauGamma*Hy,
                      aHz = tauGamma*Hz;

         // Rotation around magnetic field vector: constants.
         double u = px + py*aHz - pz*aHy,
                v = py + pz*aHx - px*aHz,
                w = pz + px*aHy - py*aHx;

         const double b_xx = aHx*aHx,
                      b_yy = aHy*aHy,
                      b_zz = aHz*aHz,
                      b_xy = aHx*aHy,
                      b_yz = aHy*aHz,
                      b_xz = aHz*aHx;
         const double delta = 1.0/(1 + b_xx + b_yy + b_zz);

         // Final step (rotation and last half acceleration by E).
         p->vx = tauMass*Ex + ((1 + b_xx)*u   + (b_xy + aHz)*v + (b_xz - aHy)*w)*delta;
         p->vy = tauMass*Ey + ((b_xy - aHz)*u + (1 + b_yy)*v   + (b_yz + aHx)*w)*delta;
         p->vz = tauMass*Ez + ((b_xz + aHy)*u + (b_yz - aHx)*v +   (1 + b_zz)*w)*delta;

         const double tauToGamma = tau/sqrt (1 + p->vx*p->vx + p->vy*p->vy + p->vz*p->vz);
         const double x2 = p->x + p->vx*tauToGamma,
                      y2 = p->y + p->vy*tauToGamma,
                      z2 = p->z + p->vz*tauToGamma;

         // Evaluates current density using r1 and r2. The most possible
         // branches are here (to remove function calls) while unused
         // branches are kept out of the main loop code.
         i2 = MF_DBL_TO_INT (x2/h1);
         j2 = MF_DBL_TO_INT (y2/h2);
         k2 = MF_DBL_TO_INT (z2/h3);
         di = abs (i2 - i);
         dj = abs (j2 - j);
         dk = abs (k2 - k);

         // Due to very rare round-off error it is possible to get y1 = y2
         // and dj = 1 inconsistency. Chanses are negligible but it is
         // possible, so this protection is added against any noise-triggered
         // jumps of this type (see example at 'arsenal/341/save#06').
         // It may be a consequence of adding and subtracting shift in
         // MF_DBL_TO_INT to avoid the changing of rounding direction for
         // negative x/y/z.
         // Test it later when 'fesetround(3)' will be implemented in gcc.
         if (fabs (x2 - x1) < 1e-10*h1 && mc_have_x) {
            di = 0;
            i2 = i;
         }

         if (fabs (y2 - y1) < 1e-10*h2 && mc_have_y) {
            dj = 0;
            j2 = j;
         }

         if (fabs (z2 - z1) < 1e-10*h3 && mc_have_z) {
            dk = 0;
            k2 = k;
         }

         if (di + dj + dk == 0) {
            // Single cell motion.
            plasmaCurrents_kernel (i, j, k, x1, x2, y1, y2, z1, z2, p->rho);
         } else if (di + dj + dk == 1) {
            // 2 stage pass (one boundary crossing).
            // Construction in denominator is a cheap insurance against
            // divizion by zero.
            const double alpha = dk*((k + k2 + 1)*0.5*h3 - z1)/(z2 - z1 + ((1 - dk)<<16)) \
                               + dj*((j + j2 + 1)*0.5*h2 - y1)/(y2 - y1 + ((1 - dj)<<16)) \
                               + di*((i + i2 + 1)*0.5*h1 - x1)/(x2 - x1 + ((1 - di)<<16));

            const double xC = x1 + alpha*(x2 - x1),
                         yC = y1 + alpha*(y2 - y1),
                         zC = z1 + alpha*(z2 - z1);

            plasmaCurrents_kernel (i,   j,  k, x1, xC, y1, yC, z1, zC, p->rho);
            plasmaCurrents_kernel (i2, j2, k2, xC, x2, yC, y2, zC, z2, p->rho);
         } else {
            plasmaCurrents_multiCellStep (i, j, k, i2, j2, k2, x1, y1, z1, x2, y2, z2, p->rho);
         }

         // Updates coordinates and sends particle to router.
#if mc_have_x
         p->x = x2;
#endif

#if mc_have_y
         p->y = y2;
#endif

#if mc_have_z
         p->z = z2;
#endif
         // Checks if particle is in core region.
         if (mc_have_x*(x2 - coreX1)*(coreX2 - x2) >= 0
          && mc_have_y*(y2 - coreY1)*(coreY2 - y2) >= 0
          && mc_have_z*(z2 - coreZ1)*(coreZ2 - z2) >= 0) {
            if (doing_shell) {
               *(--core) = *p;	// Sends particle to core.
               *p = *(--shell);	// Replaces by unprocessed particle.
               --p;		// Decrements counter to process new marker.
            }
         } else {
            if (!doing_shell) {
               *(shell++) = *p; // Sends particle to shell.
               *p = *(core++);  // Fills emptied cell by particle from tail.
            }
         }
      }
      // Updates counters.
      countShell = shell - plasma;
      countCore  = core  - plasma;

      // Tests inboxes while we are doing (massive) core.
      if (!doing_shell)
         comm_plasma_test_inbox ();
   }
}

// ---------------------------------------------------------------------------
/// Performs time step for the particles and starts parallel exchange.
// ---------------------------------------------------------------------------
void
plasma_move (meshVec_RO_p E, meshVec_RO_p H, meshVec_p J)
{
   profiler_begin (mc_prof_plasma_cleanJ);
   plasma_J = J;
   mf_mesh_clean (plasma_J);
   plasmaWx = plasmaWy = plasmaWz = 0;

#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
   fwrite (&Time, sizeof(Time), 1, tracer);
#endif

   // Remembers initial core size and makes shell timestep.
   long int unprocessedCore = countCore;
   plasma_timestep (DOING_SHELL, countShell, E, H, J);

   // Applies local boundary conditions.
   profiler_endBegin (mc_prof_plasma_pbcPop);
   for (int b = 0 ; b < 6 ; ++b)
      BC[b][cpu_bc_min[b]] ();

   // Sends outgoing particles.
   comm_plasma_send ();

   // Parallel exchange (particles and currents).
   profiler_endBegin (mc_prof_plasma_jbc);
   jbc_start (mcast_meshVec (plasma_J));

   // Final timestep (unprocessed core particles).
   plasma_timestep (DOING_CORE, countAll - unprocessedCore, E, H, J);

#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
   marker_t p = {.id = 0};
   fwrite (&p, sizeof(p), 1, tracer);
#endif

   profiler_end ();
}

// ---------------------------------------------------------------------------
/// Reports thermal energy of plasma.
// ---------------------------------------------------------------------------
void
plasma_temperature (double *WTx, double *WTy, double *WTz)
{
   *WTx = plasmaWx;
   *WTy = plasmaWy;
   *WTz = plasmaWz;
}

// ---------------------------------------------------------------------------
/// Simple charge density computing, no parallel scheduling/optimizations.
// ---------------------------------------------------------------------------
void
plasma_rho (meshDouble_p rho)
{
   mf_mesh_clean (rho);
   plasmaRho_add (rho, plasma, countShell);
   plasmaRho_add (rho, plasma + countCore, countAll - countCore);
   plasmaRho_startExchange  (rho);
   plasmaRho_finishExchange (rho);
}

#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
static void
close_tracer_file (void)
{
   fclose (tracer);
}
#endif

// ---------------------------------------------------------------------------
/// Plasma module initialization.
// ---------------------------------------------------------------------------
void
plasma_init (void)
{
   ENSURE ((1 << 16) > fmax (h1, fmax (h2, h3)),
           "plasmaCurrents_multiCellStep is not 'division by zero' safe");

   ENSURE (1e12*MC_EPS > fmax (dmn_max[0], fmax (dmn_max[1], dmn_max[2])),
           "round-off errors may appear (big domain size + small MC_EPS)");

   // Checks that at least for x > -5*h1, y > .., MF_DBL_TO_INT works nicely.
   ENSURE (MF_DBL_TO_INT (dmn_min[0] - 5.1) < dmn_min[0] - 5.1 &&
           MF_DBL_TO_INT (dmn_min[1] - 5.1) < dmn_min[1] - 5.1 &&
           MF_DBL_TO_INT (dmn_min[2] - 5.1) < dmn_min[2] - 5.1,
           "broken MF_DBL_TO_INT");

#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
   SAY_WARNING ("tracer activated (additional field is added to each marker)");
   char *mode = (Time > 0.9*tau) ? "ab" : "wb";
   tracer = cfg_open (_("binData/tracer_%03d.bin", cpu_here), mode, __func__);
   ENSURE(!atexit(close_tracer_file), "cannot submit tracer file finalizator");
#endif

   // Sets all || exchange infrastructure.
   plasma_setupConnections    ();
   plasmaRho_setupConnections ();

   // XXX  placer_exchange ();	// Forces all particles to the hosts' cpus.
   MPI_Barrier (MPI_COMM_WORLD);

   // Sets size of the core region.
   coreX1 = (cpu_min[0] + 3)*h1*mc_have_x;
   coreY1 = (cpu_min[1] + 3)*h2*mc_have_y;
   coreZ1 = (cpu_min[2] + 3)*h3*mc_have_z;
   coreX2 = (cpu_max[0] - 3)*h1*mc_have_x;
   coreY2 = (cpu_max[1] - 3)*h2*mc_have_y;
   coreZ2 = (cpu_max[2] - 3)*h3*mc_have_z;

   shellX1 = cpu_min[0]*h1*mc_have_x;
   shellY1 = cpu_min[1]*h2*mc_have_y;
   shellZ1 = cpu_min[2]*h3*mc_have_z;
   shellX2 = cpu_max[0]*h1*mc_have_x;
   shellY2 = cpu_max[1]*h2*mc_have_y;
   shellZ2 = cpu_max[2]*h3*mc_have_z;

   // Makes sure we have enough free space.
   if (countCore - countShell < MC_BUFFER/sizeof (marker_t)) {
      long int increment = MC_BUFFER/sizeof (marker_t) + 1;
      plasma = (marker_t*) realloc (plasma,
                                    sizeof (marker_t)*(countAll + increment));
      ENSURE (plasma, "cannot add %ld particles to buffer", increment);
      memmove (plasma + countCore + increment,
               plasma + countCore, (countAll - countCore)*sizeof (marker_t));
      countAll  += increment;
      countCore += increment;
   }

   SAY_DEBUG ("Plasma buffer: [%p, %p), %ld particles, %ld free.",
              plasma, plasma + countAll, countShell + countAll - countCore,
              countCore - countShell);

   // XXX: make reg_t constructor and use 'reg_inside ()' tester to test
   // condition in one line.
   ENSURE (cpu_max[0] - cpu_min[0] >= 7*mc_have_x &&
           cpu_max[1] - cpu_min[1] >= 7*mc_have_y &&
           cpu_max[2] - cpu_min[2] >= 7*mc_have_z,
           "core region is too small.");

   // Turns off all unused boundaries and ensures boundary condition.
   for (int i = 0 ; i < 6 ; ++i) {
      for (int j = 0 ; j < BC_ENUM_LENGHT && !ACTIVATOR[i] ; ++j)
	 BC[i][j] = dummy;
      BC[i][cpu_bc_min[i]] ();
   }

   comm_plasma_configure ();	// Initializes parallel exchange frame.
   comm_plasma_dump      ();	// Logs configuration.
}
