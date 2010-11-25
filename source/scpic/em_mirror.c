/** \file em_mirror.c
  * Source of EM waves: laser beam focused by parabolic mirror.
  *
  * See header for documentation and credits.
  *
  * \todo Use new mesh container to wrap all arrays of 'wave_x_t' type to
  *       automate boundary checks/unrolling.
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include <mpi.h>

#include "em_mirror.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

/// Set this flag to '1' to debugs code with plane wave source.
#define JUST_PLANE_WAVE 0

/// Offset between domain boundary and source; make it bigger than cap size to
/// keep real domain boundary code unmodified --- periodic and other physical
/// boundary condition code doesn't account for soft sources and one can lose
/// or screw the corresponding contributions otherwise.
#define OFFSET          5

/// X index of the laser soft source (read the previous comment too).
#define SOURCE_X_OFFSET 5

/// Set this flags to '0' to disable using symmetry of Stratton-Chu integrals.
#define USE_Y_SYMM 	1
#define USE_Z_SYMM 	1

/// Local wrapper to return two packed values.
typedef struct {
   double sin;
   double cos;
} vec2d_t;

static void compute_Stratton_Chu_integrals (mirror_t *src);

// ----------------------------------------------------------------------------
/// Calculates 'Ex' field at (xp, yp, zp) point as Stratton-Chu integral;
/// incoming pulse is plane/gaussian wave, best focus is at (0, 0, 0) point.
///
/// \warning This function is NOT USED but kept for the future.
/// \warning It doesn't feature 'sin + cos' optimization (simultaneous
///          integration for 'τ' and 'τ + π/2' to reuse heavy funcs).
// ----------------------------------------------------------------------------
double
focused_ex (double xp, double yp, double zp, double tau, mirror_t p)
{
   // Going to the mirror's frame.
   xp += (p.f0 - p.rm*p.rm/(4*p.f0));

#if JUST_PLANE_WAVE
   return 0;
#endif

   double h    = p.rm/p.NN,
          dydz = 4*h*h;

   double sum = 0,
          rr  = 1;
   for (int j = -p.NN/2 ; j < p.NN/2 ; j++)
   for (int k = -p.NN/2 ; k < p.NN/2 ; k++) {
      double y = 2*(j + 0.5)*h,
             z = 2*(k + 0.5)*h,
             x = (y*y + z*z - p.rm*p.rm)/(4*p.f0);
      if (p.rratio > 0)
         rr = exp (-(y*y + z*z)/(p.rm*p.rm*p.rratio*p.rratio));
      if (x <= 0) {
         double u   = sqrt ((x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp)),
                phi = 2*mc_pi*(x - u) + tau,
                ur  = 1.0/u;
         sum += rr*y*(sin (phi)*2*mc_pi*(1 - (x - xp)*ur)*ur +
                      cos (phi)*(x - xp)*ur*ur*ur);
      }
   }

   return sum*dydz/(4*mc_pi*p.f0);
}

// ----------------------------------------------------------------------------
/// Calculates 'Ey' field at (xp, yp, zp) point as Stratton-Chu integral;
/// incoming pulse is plane/gaussian wave, best focus is at (0, 0, 0) point.
/// Result is computed for 'τ' (in '.cos') and for 'τ+π/2' (in '.sin') fields.
// ----------------------------------------------------------------------------
static vec2d_t
focused_ey (double xp, double yp, double zp,
            double tau, mirror_t p)
{
   // Going to the mirror's frame.
   xp += (p.f0 - p.rm*p.rm/(4*p.f0));

#if JUST_PLANE_WAVE == 1
   return (vec2d_t) { .cos = cos (2*mc_pi*xp + tau),
                      .sin = sin (2*mc_pi*xp + tau) };
#endif

   double h    = p.rm/p.NN,
          dydz = 4*h*h;

   double  rr  = 1;
   vec2d_t sum = { .cos = 0,
                   .sin = 0 };
   for (int j = -p.NN/2 ; j < p.NN/2 ; j++)
   for (int k = -p.NN/2 ; k < p.NN/2 ; k++) {
      double y = 2*(j + 0.5)*h,
             z = 2*(k + 0.5)*h,
             x = (y*y + z*z - p.rm*p.rm)/(4*p.f0);
      if (p.rratio > 0)
         rr = exp (-(y*y + z*z)/(p.rm*p.rm*p.rratio*p.rratio));
      if (x <= 0) {
         double u   = sqrt ((x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp)),
                phi = 2*mc_pi*(x - u) + tau,
                ur  = 1.0/u;
         sum.cos += rr*(sin (phi)*2*mc_pi*(2*p.f0 - y*(y - yp)*ur)*ur
                                             + cos (phi)*y*(y - yp)*ur*ur*ur);
         sum.sin += rr*(cos (phi)*2*mc_pi*(2*p.f0 - y*(y - yp)*ur)*ur
                                             - sin (phi)*y*(y - yp)*ur*ur*ur);
      }
   }
   sum.cos *= dydz/(4*mc_pi*p.f0);
   sum.sin *= dydz/(4*mc_pi*p.f0);

   return sum;
}

// ----------------------------------------------------------------------------
/// Calculates 'Ez' field at (xp, yp, zp) point as Stratton-Chu integral;
/// incoming pulse is plane/gaussian wave, best focus is at (0, 0, 0) point.
/// Result is computed for 'τ' (in '.cos') and for 'τ+π/2' (in '.sin') fields.
// ----------------------------------------------------------------------------
static vec2d_t
focused_ez (double xp, double yp, double zp,
            double tau, mirror_t p)
{
   // Going to the mirror's frame.
   xp += (p.f0 - p.rm*p.rm/(4*p.f0));

#if JUST_PLANE_WAVE == 1
   return (vec2d_t) { .cos = 0,
                      .sin = 0 };
#endif

   double h    = p.rm/p.NN,
          dydz = 4*h*h;

   double  rr  = 1;
   vec2d_t sum = { .cos = 0,
                   .sin = 0 };
   for (int j = -p.NN/2 ; j < p.NN/2 ; j++)
   for (int k = -p.NN/2 ; k < p.NN/2 ; k++) {
      double y = 2*(j + 0.5)*h,
             z = 2*(k + 0.5)*h,
             x = (y*y + z*z - p.rm*p.rm)/(4*p.f0);
      if (p.rratio > 0)
         rr = exp (-(y*y + z*z)/(p.rm*p.rm*p.rratio*p.rratio));
      if (x <= 0) {
         double u   = sqrt ((x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp)),
                phi = 2*mc_pi*(x - u) + tau,
                ur  = 1.0/u,
                amp = rr*y*(z - zp)*ur*ur;
         sum.cos -= amp*(sin (phi)*2*mc_pi - cos (phi)*ur);
         sum.sin -= amp*(cos (phi)*2*mc_pi + sin (phi)*ur);
      }
   }
   sum.cos *= dydz/(4*mc_pi*p.f0);
   sum.sin *= dydz/(4*mc_pi*p.f0);

   return sum;
}

// ----------------------------------------------------------------------------
/// Calculates 'Hx' field at (xp, yp, zp) point as Stratton-Chu integral;
/// incoming pulse is plane/gaussian wave, best focus is at (0, 0, 0) point.
/// Result is computed for 'τ' (in '.cos') and for 'τ+π/2' (in '.sin') fields.
// ----------------------------------------------------------------------------
vec2d_t
focused_hx (double xp, double yp, double zp,
            double tau, mirror_t p)
{
   // Going to the mirror's frame.
   xp += (p.f0 - p.rm*p.rm/(4*p.f0));

#if JUST_PLANE_WAVE == 1
   return (vec2d_t) { .cos = 0,
                      .sin = 0 };
#endif

   double h    = p.rm/p.NN,
          dydz = 4*h*h;

   double  rr  = 1;
   vec2d_t sum = { .cos = 0,
                   .sin = 0 };
   for (int j = -p.NN/2 ; j < p.NN/2 ; j++)
   for (int k = -p.NN/2 ; k < p.NN/2 ; k++) {
      double y = 2*(j + 0.5)*h,
             z = 2*(k + 0.5)*h,
             x = (y*y + z*z - p.rm*p.rm)/(4*p.f0);
      if (p.rratio > 0)
         rr = exp (-(y*y + z*z)/(p.rm*p.rm*p.rratio*p.rratio));
      if (x <= 0) {
         double u   = sqrt ((x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp)),
                phi = 2*mc_pi*(x - u) + tau,
                ur  = 1.0/u,
                amp = 2*rr*p.f0*(z - zp)*ur*ur;
         sum.cos += amp*(sin (phi)*2*mc_pi - cos (phi)*ur);
         sum.sin += amp*(cos (phi)*2*mc_pi + sin (phi)*ur);
      }
   }
   sum.cos *= dydz/(4*mc_pi*p.f0);
   sum.sin *= dydz/(4*mc_pi*p.f0);

   return sum;
}

// ----------------------------------------------------------------------------
/// Calculates 'Hy' field at (xp, yp, zp) point as Stratton-Chu integral;
/// incoming pulse is plane/gaussian wave, best focus is at (0, 0, 0) point.
/// Result is computed for 'τ' (in '.cos') and for 'τ+π/2' (in '.sin') fields.
// ----------------------------------------------------------------------------
static vec2d_t
focused_hy (double xp, double yp, double zp,
            double tau, mirror_t p)
{
   return focused_ez (xp, yp, zp, tau, p);
}

// ----------------------------------------------------------------------------
/// Calculates 'Hz' field at (xp, yp, zp) point as Stratton-Chu integral;
/// incoming pulse is plane/gaussian wave, best focus is at (0, 0, 0) point.
/// Result is computed for 'τ' (in '.cos') and for 'τ+π/2' (in '.sin') fields.
// ----------------------------------------------------------------------------
static vec2d_t
focused_hz (double xp, double yp, double zp,
            double tau, mirror_t p)
{
   // Going to the mirror's frame.
   xp += (p.f0 - p.rm*p.rm/(4*p.f0));

#if JUST_PLANE_WAVE == 1
   return (vec2d_t) { .cos = cos (2*mc_pi*xp + tau),
                      .sin = sin (2*mc_pi*xp + tau) };
#endif

   double h    = p.rm/p.NN,
          dydz = 4*h*h;

   double  rr  = 1;
   vec2d_t sum = { .cos = 0,
                   .sin = 0 };
   for (int j = -p.NN/2 ; j < p.NN/2 ; j++)
   for (int k = -p.NN/2 ; k < p.NN/2 ; k++) {
      double y = 2*(j + 0.5)*h,
             z = 2*(k + 0.5)*h,
             x = (y*y + z*z - p.rm*p.rm)/(4*p.f0);
      if (p.rratio > 0)
         rr = exp (-(y*y + z*z)/(p.rm*p.rm*p.rratio*p.rratio));
      if (x <= 0) {
         double u   = sqrt ((x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp)),
                phi = 2*mc_pi*(x - u) + tau,
                ur  = 1.0/u,
                amp = rr*(y*(y - yp) - 2*p.f0*(x - xp))*ur*ur;
         sum.cos += amp*(2*mc_pi*sin (phi) - ur*cos (phi));
         sum.sin += amp*(2*mc_pi*cos (phi) + ur*sin (phi));
      }
   }
   sum.cos *= dydz/(4*mc_pi*p.f0);
   sum.sin *= dydz/(4*mc_pi*p.f0);

   return sum;
}

// █████████████████████████████████████████████████████████████████████████ //
// █████████████████████████████████████████████████████████████████████████ //
// █████████████████████████████████████████████████████████████████████████ //

// ----------------------------------------------------------------------------
/// Reads body of the tag; used in both 'setup.out' and in 'core.out'.
// ----------------------------------------------------------------------------
mirror_t
read_mirror_parameters (FILE *fp)
{
   mirror_t src = { .E = NULL,
                    .H = NULL };

   // Reads 'A' and converts it to dimensionless units.
   src.A = cfg_readDouble (fp);
   if (src.A > 0) {
      src.A = units (mc_r0)
              *mc_CGS_e
              *sqrt (8*mc_pi*src.A*mc_CGS_wattPerSquareCantimeter/mc_CGS_c)
              /(mc_CGS_m*mc_CGS_c*mc_CGS_c);
   } else {
      src.A *= (-2*mc_pi);
   }

   src.f0  = cfg_readDouble (fp);
   src.rm  = 4.0*sqrt (Ly*Ly + Lz*Lz);
   src.f0 *= 2*src.rm;
   src.NN  = ((int) (0.5/sqrt (h1*h1 + h2*h2 + h3*h3)))*((int) src.rm);

   src.xf  = cfg_readDouble (fp)*1e-4/units (mc_r0);	// μm -> r₀ conversion.
   src.yf  = cfg_readDouble (fp)*1e-4/units (mc_r0);
   src.zf  = cfg_readDouble (fp)*1e-4/units (mc_r0);
   src.dt  = cfg_readDouble (fp)*1e-15/units (mc_t0);	// fs -> t₀ conversion.
   src.T0  = cfg_readDouble (fp)*1e-15/units (mc_t0);

   // 'sqrt(2/log(2))' emerges from conversion of intensity FWHM to field.
   src.dt *= sqrt (2.0/log (2));

   src.phi0      = cfg_readDouble (fp);
   src.rratio    = cfg_readDouble (fp);
   src.X_forward = cfg_readInt    (fp) > 0;
   src.delay     = src.X_forward ? src.xf     - SOURCE_X_OFFSET*h1
                                 : dmn_max[0] - SOURCE_X_OFFSET*h1 - src.xf;
   return src;
}

// ----------------------------------------------------------------------------
/// Prints parameters of the source to the screen.
// ----------------------------------------------------------------------------
void
print_mirror_parameters (mirror_t src)
{
#if JUST_PLANE_WAVE == 1
   SAY_WARNING ("Adding a proxy (plane wave) instead of laser beam to debug.");
   return;
#endif
   say ("Adding laser source:");
   say ("    peak intensity is %.2e W/cm²",
       pow (src.A*mc_CGS_m*mc_CGS_c*mc_CGS_c/(mc_CGS_e*units (mc_r0)), 2)
     * mc_CGS_c/(8*mc_pi*mc_CGS_wattPerSquareCantimeter));
   say ("    incident wave is %s", (src.rratio == 0) ? "plane wave"
                                                     : "gaussian beam");
   say ("    mirror focus is in (%.2f, %.2f, %.2f) microns",
                                                      src.xf*1e4*units (mc_r0),
                                                      src.yf*1e4*units (mc_r0),
                                                      src.zf*1e4*units (mc_r0));
   say ("    mirror f-number is %.2f", 0.5*src.f0/src.rm);
   say ("    pulse duration is %.2lf fs, pulse arrives into focus in %.2lf fs",
                              src.dt*units (mc_t0)/1e-15/sqrt (2.0/log (2)),
                              src.T0*units (mc_t0)/1e-15);

   double scale = exp (- pow (2*(src.T0 - src.delay)/src.dt, 2));
   say ("    source to focus travel time is %.2lf fs", src.delay*units (mc_t0)
                                                                /1e-15);
   say ("    field amplitude on leading front is %.2le of a maximum", scale);
   say ("    direction of focused pulse is '%s X'", src.X_forward ? "positive"
                                                                  : "negative");
   say ("    ... more details are saved as debug info, check logs.");

   SAY_DEBUG ("Parameters of the laser source:");
   SAY_DEBUG ("  o maximum intensity in the best focus: %.2e W/cm^2",
       pow (src.A*mc_CGS_m*mc_CGS_c*mc_CGS_c/
            (mc_CGS_e*units(mc_r0)), 2) * mc_CGS_c / (8*mc_pi*mc_CGS_wattPerSquareCantimeter));
   SAY_DEBUG ("  o a = e E_max(best focus)/(m omega c) = %.2f", src.A/(2*mc_pi));
   SAY_DEBUG ("  o mirror parameters:");
   SAY_DEBUG ("    - f-number = %.2f", 0.5*src.f0/src.rm);
   SAY_DEBUG ("    - actual focal length of the mirror: %.1f microns", 1e4 * src.f0 * units(mc_r0));
   SAY_DEBUG ("    - actual diameter of the mirror: %.1f microns", 2e4 * src.rm * units(mc_r0));
   SAY_DEBUG ("    - mirror surface integration over %d x %d points", src.NN, src.NN);
   SAY_DEBUG ("  o best focus position:");
   SAY_DEBUG ("    - xf = %.2f micron", src.xf*1e4*units (mc_r0));
   SAY_DEBUG ("    - yf = %.2f micron", src.yf*1e4*units (mc_r0));
   SAY_DEBUG ("    - zf = %.2f micron", src.zf*1e4*units (mc_r0));
   SAY_DEBUG ("  o pulse duration: %.2f fs", src.dt*units (mc_t0)/1e-15/sqrt (2./log(2)));
   SAY_DEBUG ("  o pulse offset: %.2f fs",   src.T0*units (mc_t0)/1e-15);
   SAY_DEBUG ("  o phase shift: %.2f radian", src.phi0);
   if(src.rratio == 0)
      SAY_DEBUG ("  o pure plane wave incident onto the mirror");
   else
      SAY_DEBUG ("  o pulse of Gaussian transverse profile is incident onto the mirror: r_pulse/r_mirror = %.2f",
           src.rratio);
   SAY_DEBUG ("  o direction of propagation: %cx", src.X_forward ? '+' : '-');
}

// ----------------------------------------------------------------------------
/// Adds new laser beam focused by parabolic mirror.
/// Used by setup to compute and save data on HDD (see 'load_mirrors').
// ----------------------------------------------------------------------------
void
save_mirror_source (int mirror_number, mirror_t mirror)
{
   compute_Stratton_Chu_integrals (&mirror);

   // Master node saves the file.
   if (!cpu_here) {
      double to_microns  = units (mc_r0)/1e-4,
             to_femtosec = units (mc_t0)/1e-15;

      const char *filename = _(".EM_sources/laser/mirror_%02d", mirror_number);
      FILE *fp = cfg_open (filename, "wb", __func__);
      fprintf (fp, "Mirror focused laser beam #%d\n",  mirror_number);
      fprintf (fp, "---------------------------------------\n");
      fprintf (fp, "@ %.15le   Poynting vector at the best focus: positive => I [W/cmÂ²] / negative => e*E/(m*omega*c).\n", mirror.A/(-2*mc_pi));
      fprintf (fp, "@ %.15le   Mirror f-number (defined as F_mirror/D_mirror; the smaller the number, the tighter the focusing). The value has to be greater then 0.25.\n", mirror.f0/(2*mirror.rm));
      fprintf (fp, "@ %.15le   Focus X coordinate [micron].\n", mirror.xf*to_microns);
      fprintf (fp, "@ %.15le   Focus Y coordinate [micron].\n", mirror.yf*to_microns);
      fprintf (fp, "@ %.15le   Focus Z coordinate [micron].\n", mirror.zf*to_microns);
      fprintf (fp, "@ %.15le   Pulse duration, [fs] (defined by half of the maximum intensity), the envelope is Gaussian.\n", mirror.dt*to_femtosec/sqrt (2.0/log (2)));
      fprintf (fp, "@ %.15le   Pulse to focus arrival time [fs].\n", mirror.T0*to_femtosec);
      fprintf (fp, "@ %.15le   Phase shift.\n", mirror.phi0);
      fprintf (fp, "@ %.15le   Ratio between the radius of the pulse incident onto mirror and radius of the mirror (-1 to incident plane wave).\n", mirror.rratio);
      fprintf (fp, "@ %d       Direction of propagation (positive or negative).\n", mirror.X_forward ? +1 : -1);
      fprintf (fp, "---------------------------------------\n");

      long int nodes = (dmn_max[1] + 1 - 2*OFFSET)*(dmn_max[2] + 1 - 2*OFFSET);
      fwrite  (mirror.E, sizeof (wave_x_t), nodes, fp);
      fwrite  (mirror.H, sizeof (wave_x_t), nodes, fp);
      fclose  (fp);
   }

   // Releases memory to unburden 'setup.out'.
   free (mirror.E);
   free (mirror.H);
   mirror.E = mirror.H = NULL;
}

// ----------------------------------------------------------------------------
/// Allocates and loads subregion of mirror_t::E/H array from monolithic file.
// ----------------------------------------------------------------------------
static wave_x_t*
load_mirror_subregion (FILE *fp)
{
   int jmin = fmax (cpu_min[1] - 1, OFFSET             ) + 0.1,
       jmax = fmin (cpu_max[1] + 1, dmn_max[1] - OFFSET) + 0.1,
       kmin = fmax (cpu_min[2] - 1, OFFSET             ) + 0.1,
       kmax = fmin (cpu_max[2] + 1, dmn_max[2] - OFFSET) + 0.1;

   long int  nodes = (jmax - jmin + 1)*(kmax - kmin + 1);
   wave_x_t *dest  = calloc (nodes, sizeof (wave_x_t));
   ENSURE (dest, "cannot allocate memory for the source");

   // Simple strides and offsets to address points inside saved array.
   #define FPOS(j, k) sizeof(wave_x_t)*\
     ( (j - OFFSET)*(dmn_max[2] + 1 - 2*OFFSET) + (k - OFFSET))
   long int origin = ftell (fp);
   int      pos    = 0,
            line   = kmax - kmin + 1;
   for (int j = jmin ; j <= jmax ; ++j) {
      fseek (fp, origin + FPOS(j, kmin), SEEK_SET);
      fread (dest + pos, sizeof (wave_x_t), line, fp);
      pos += line;
   }
   // Pretends we read the entire global array: skips to the next chunk.
   fseek (fp,
          origin + FPOS (dmn_max[1] + 1 - OFFSET, dmn_min[2] + OFFSET),
          SEEK_SET);
   return dest;
}

// ----------------------------------------------------------------------------
/// Loads prepared mirror source or indicates 'no source' error.
// ----------------------------------------------------------------------------
int
mirror_source_is_loaded (int mirror_number, mirror_t *dest)
{
   const char *filename = _(".EM_sources/laser/mirror_%02d", mirror_number);
   FILE *fp = fopen (filename, "rb");
   if (!fp)
      return 0;

   // Reads header.
   fscanf (fp, "Mirror focused laser beam #%*d\n");
   fscanf (fp, "---------------------------------------\n");
   *dest = read_mirror_parameters (fp);
   fscanf (fp, "---------------------------------------\n");

   // Allocates storage and loads precomputed source.
   dest->E = load_mirror_subregion (fp);
   dest->H = load_mirror_subregion (fp);
   fclose (fp);

   return 1;
}

// █████████████████████████████████████████████████████████████████████████ //
// █████████████████████████████████████████████████████████████████████████ //
// █████████████████████████████████████████████████████████████████████████ //

// ----------------------------------------------------------------------------
/// Computes Stratton-Chu integrals (keeps both amplitude and phase info).
/// Uses (anti)symmetry to reuse computed values if possible.
// ----------------------------------------------------------------------------
static void
compute_Stratton_Chu_integrals (mirror_t *src)
{
   struct timeval start, end;		// XXX Use my timer after refactoring.
   gettimeofday (&start, NULL);

   // Defines and allocates global array to store the entire source.
   int      JMIN  = OFFSET,                  jmin = JMIN,
            JMAX  = dmn_max[1] - OFFSET,     jmax = JMAX,
            KMIN  = OFFSET,                  kmin = KMIN,
            KMAX  = dmn_max[2] - OFFSET,     kmax = KMAX;
   long int nodes = (JMAX - JMIN + 1)*(KMAX - KMIN + 1);
   src->E = calloc (nodes, sizeof (wave_x_t));
   src->H = calloc (nodes, sizeof (wave_x_t));
   ENSURE (src->E && src->H,
           "cannot allocate %.2fMb of RAM", 2*nodes*sizeof (wave_x_t)/1.049e6);
   if (sizeof (wave_x_t)*2*nodes*cpu_total > 2e9)
      SAY_WARNING ("Algorithm uses too much RAM! Mail me and I'll fix that.");

   // Activates symmetry passes and defines limits of a biggest independent
   // quarter of YZ coordinate plane.
   const int do_y = fabs (fmod (src->yf, 0.5*h2)) < 1e-12*h2 && USE_Y_SYMM,
             do_z = fabs (fmod (src->zf, 0.5*h3)) < 1e-12*h3 && USE_Z_SYMM;
   if (do_y) {
      if (src->yf < h2*(JMIN + JMAX)/2) {
         jmin = fmax (src->yf/h2 - 1, JMIN) + 0.1;
      } else {
         jmax = fmin (src->yf/h2 + 1, JMAX) + 0.1;
      }
   }
   if (do_z) {
      if (src->zf < 0.5*h3*(KMIN + KMAX)) {
         kmin = fmax (src->zf/h3 - 1, KMIN) + 0.1;
      } else {
         kmax = fmin (src->zf/h3 + 1, KMAX) + 0.1;
      }
   }
   if (do_y || do_z) {
      say ("    ... using '%s%s' (anti)symmetry to speed-up", do_y ? "Y" : "",
                                                              do_z ? "Z" : "");
   }
   SAY_DEBUG ("full region: %3d %3d - %3d %3d", JMIN, JMAX, KMIN, KMAX);
   SAY_DEBUG ("sub region:  %3d %3d - %3d %3d", jmin, jmax, kmin, kmax);

   // Calibrates amplitude of incident wave to get given intensity in focal
   // point; also finds 'zero phase', which corresponds to this maximum of Ey.
   // Mathematical proof is on website.
   vec2d_t ey   = focused_ey (0, 0, 0, 0, *src);
   double  phi0 = atan2 (ey.sin, ey.cos),
           ampl = sqrt (ey.cos*ey.cos + ey.sin*ey.sin);
   ENSURE (ampl, "bad source field, 'max(E_y(t)) = 0' is impossible");
   double a = src->A/ampl;

   // Parallel Stratton-Chu integration for each node of the source.
   double    xf = src->xf,
             yf = src->yf,
             zf = src->zf,
             dx = 0.5*h1,
             dy = 0.5*h2,
             dz = 0.5*h3;
   wave_x_t *E  = src->E,
            *H  = src->H;
   for (int j = jmin ; j <= jmax ; j++)
   for (int k = kmin ; k <= kmax ; k++) {
#define POS(J, K) ((J - JMIN)*(KMAX - KMIN + 1) + K - KMIN)
      long int cnt = POS(j, k);
      if (cnt % cpu_total != cpu_here)
         continue;

#define PACK(dest, scale, func) {	\
   vec2d_t res = func;			\
   dest ## _cos = (scale)*res.cos;	\
   dest ## _sin = (scale)*res.sin;	\
}
      double y = j*h2 - yf,
             z = k*h3 - zf;
      if (src->X_forward) {
         double x = SOURCE_X_OFFSET*h1 - xf;
         PACK(E[cnt].y,  a, focused_ey (x-h1, y-dy, z,    phi0, *src));
         PACK(E[cnt].z, -a, focused_ez (x-h1, y,    z-dz, phi0, *src));
         PACK(H[cnt].y, -a, focused_hy (x-dx, y,    z-dz, phi0, *src));
         PACK(H[cnt].z,  a, focused_hz (x-dx, y-dy, z,    phi0, *src));
      } else {
         double x = xf - (Lx - SOURCE_X_OFFSET*h1);
         PACK(E[cnt].y, -a, focused_ey (x-h1, y-dy, z,    phi0, *src));
         PACK(E[cnt].z,  a, focused_ez (x-h1, y,    z-dz, phi0, *src));
         PACK(H[cnt].y, -a, focused_hy (x-dx, y,    z-dz, phi0, *src));
         PACK(H[cnt].z,  a, focused_hz (x-dx, y-dy, z,    phi0, *src));
      }
      double progress = 100.0*(cnt - POS(jmin, kmin))/POS(jmax, kmax);
      say_doing ("Stratton-Chu integration: %.2f%% complete", progress);
   }

   // Defragmentation of result (THE RESULT IS COLLECTED ON ROOT NODE ONLY).
   ENSURE (sizeof (wave_x_t) == 4*sizeof (double),
           "'wave_x_t' is not an array of doubles");
   if (!cpu_here) {
// Good variant, but is replaced with a kludge below.
//       MPI_Reduce (MPI_IN_PLACE, E,    nodes*4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//       MPI_Reduce (MPI_IN_PLACE, H,    nodes*4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      // KLUDGE: LAMMPI in FIAN doesn't support MPI_IN_PLACE => I reduce using 
      // tmp array. Anyway, this kludge will be removed.
      // \todo XXX remove this kludge after 01.01.2011.
      wave_x_t *E_send = calloc (nodes, sizeof (wave_x_t)),
               *H_send = calloc (nodes, sizeof (wave_x_t));
      ENSURE (E_send && H_send, "cannot allocate memory for MPI_Reduce kludge");
      memcpy (E_send, E, sizeof (wave_x_t)*nodes);
      memcpy (H_send, H, sizeof (wave_x_t)*nodes);
      MPI_Reduce (E_send, E, nodes*4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce (H_send, H, nodes*4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      free (E_send);
      free (H_send);
   } else {
      MPI_Reduce (E,   NULL, nodes*4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce (H,   NULL, nodes*4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   }

   // Using symmetry to fill the rest of array on the root node.
   for (int j = jmin ; j <= jmax && !cpu_here && (do_y || do_z) ; j++)
   for (int k = kmin ; k <= kmax                                ; k++) {
      long int cnt = POS(j, k);

// Writes data on array with "out of the boundary" tests.
#define WRITE(J, K, ARR, COMP, SCALE) 					\
   if (J >= JMIN && J <= JMAX						\
   &&  K >= KMIN && K <= KMAX) {					\
      ARR[POS(J, K)]. COMP ## _cos = (SCALE)*ARR[cnt]. COMP ## _cos;	\
      ARR[POS(J, K)]. COMP ## _sin = (SCALE)*ARR[cnt]. COMP ## _sin;	\
}
      if (do_y) {
         // See section 'INDEXING OF MIRROR REFLECTION' in 'em_mirror.h'; note
         // that here 'shift' is '-1' (thus sign is changed).
         int jf = 2.0*yf/h2;
         WRITE(jf - j + 1, k, E, y,  1);
         WRITE(jf - j    , k, E, z, -1);
         WRITE(jf - j    , k, H, y, -1);
         WRITE(jf - j + 1, k, H, z,  1);
      }
      if (do_z) {
         int kf = 2.0*zf/h3;
         WRITE(j, kf - k    , E, y,  1);
         WRITE(j, kf - k + 1, E, z, -1);
         WRITE(j, kf - k + 1, H, y, -1);
         WRITE(j, kf - k    , H, z,  1);
      }
      if (do_y && do_z) {
         int jf = 2.0*yf/h2,
             kf = 2.0*zf/h3;
         WRITE(jf - j + 1, kf - k    , E, y, 1);
         WRITE(jf - j,     kf - k + 1, E, z, 1);
         WRITE(jf - j,     kf - k + 1, H, y, 1);
         WRITE(jf - j + 1, kf - k    , H, z, 1);
      }
      double progress = 100.0*(cnt - POS(jmin, kmin))/POS(jmax, kmax);
      say_doing ("Stratton-Chu symmetry pass: %.2f%%", progress);
   }
#undef POS
#undef PACK
#undef WRITE

   gettimeofday (&end, NULL);
   say_doing ("Stratton-Chu integration: done in %.2f seconds.",
       end.tv_sec - start.tv_sec + 1e-6*(end.tv_usec - start.tv_usec));
}

// ----------------------------------------------------------------------------
/// Scales and adds one array on a mesh subregion (basic soft source step).
/// \note Time dependence of intensity enters here.
// ----------------------------------------------------------------------------
static void
add_source (meshVec_p m, wave_x_t *src, int pos,
            double ampl, double phase, reg_t *to_update)
{
   if (pos < to_update->min[0] || pos > to_update->max[0])
      return;

   int    JMIN = fmax (cpu_min[1] - 1, OFFSET             ) + 0.1,
          JMAX = fmin (cpu_max[1] + 1, dmn_max[1] - OFFSET) + 0.1,
          KMIN = fmax (cpu_min[2] - 1, OFFSET             ) + 0.1,
          KMAX = fmin (cpu_max[2] + 1, dmn_max[2] - OFFSET) + 0.1;

   int    jmin = fmax (to_update->min[1], JMIN) + 0.1,
          jmax = fmin (to_update->max[1], JMAX) + 0.1,
          kmin = fmax (to_update->min[2], KMIN) + 0.1,
          kmax = fmin (to_update->max[2], KMAX) + 0.1;
   double cs   = cos (phase),
          sn   = sin (phase),
          A    = ampl*tau/h1;
   for (int j = jmin ; j <= jmax ; j++)
   for (int k = kmin ; k <= kmax ; k++) {
      long int cnt = (j - JMIN)*(KMAX - KMIN + 1) + k - KMIN;
      mv_fy(m, pos, j, k) += A*(src[cnt].z_cos*cs + src[cnt].z_sin*sn);
      mv_fz(m, pos, j, k) += A*(src[cnt].y_cos*cs + src[cnt].y_sin*sn);
   }
}

// ----------------------------------------------------------------------------
/// Adds sources of magnetic fields corresponding to the single mirror.
// ----------------------------------------------------------------------------
void
add_mirror_H (mirror_t *src, meshVec_p H, double time, reg_t *to_update)
{
   double t    = time - src->T0 + src->delay,
          ampl = exp (-pow (2*t/src->dt, 2));
   int    I    = (src->X_forward) ?              SOURCE_X_OFFSET
                                  : dmn_max[0] - SOURCE_X_OFFSET + 1;
   add_source (H, src->E, I, ampl, 2*mc_pi*t + src->phi0, to_update);
}

// ----------------------------------------------------------------------------
/// Adds sources of electric fields corresponding to the single mirror.
// ----------------------------------------------------------------------------
void
add_mirror_E (mirror_t *src, meshVec_p E, double time, reg_t *to_update)
{
   double t    = time + 0.5*tau - src->T0 + src->delay,
          ampl = exp (-pow (2*t/src->dt, 2));
   int    I    = (src->X_forward) ?              SOURCE_X_OFFSET - 1
                                  : dmn_max[0] - SOURCE_X_OFFSET + 1;
   add_source (E, src->H, I, ampl, 2*mc_pi*t + src->phi0, to_update);
}
