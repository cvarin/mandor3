/** \file em_gaussSpot.c
  * Hard source of the radiation in form of X-plane with radiating gauss-spot.
  */

#include <stdio.h>
#include <math.h>

#include <mpi.h>

#include "type_mesh.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

static int gaussSetupCheckOnly = 0;

// ---------------------------------------------------------------------------
/// Parameters of the gauss spot (packed into exchange structure).
// ---------------------------------------------------------------------------
struct gaussPack_s
{
    double ampl, frequency;
    double x, y, z;
    double width, duration, frontOffset, t_plateau;
    double Ey, Ez;
    int i;
} gauss;

// ---------------------------------------------------------------------------
/// Setups all parameters to numerically add gauss-spot type source.
// ---------------------------------------------------------------------------
void
gaussSpot_init (void)
{
    int        emptyFile = 0;
    static int first     = 1;

    ENSURE (first, "double initialization");
    first = 0;

    // Master reads and distributes all parameters.
    if (!cpu_here) {
        FILE  *fp;
        if ((fp = fopen (".gaussSpot.cfg", "rt")))
        {
            gauss.ampl      = cfg_readDouble (fp);
            gauss.frequency = 2*mc_pi*cfg_readDouble (fp);
            gauss.x         = cfg_readDouble (fp);
            gauss.y         = cfg_readDouble (fp);
            gauss.z         = cfg_readDouble (fp);
            gauss.width     = cfg_readDouble (fp);
            gauss.duration    = cfg_readDouble (fp);
            gauss.frontOffset = cfg_readDouble (fp);
            gauss.t_plateau   = cfg_readDouble (fp);
            gauss.Ey        = cfg_readInt (fp);
            gauss.Ez        = cfg_readInt (fp);
            fclose (fp);
        }
        else
            gaussSetupCheckOnly = emptyFile = 1;
    }

    MPI_Bcast (&gaussSetupCheckOnly, 1, MPI_INT, 0, MPI_COMM_WORLD);						// Master broadcasts status of the run.
    MPI_Bcast (&emptyFile, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&gauss, sizeof (gauss), MPI_BYTE, 0, MPI_COMM_WORLD);					// Broadcasts parameters of the source.

    if (emptyFile) {
        SAY_DEBUG ("gaussSpot_init: file '.gaussSpot.cfg' doesn't exist.");
        return;
    }

    gauss.Ey = (gauss.Ey > 0) - (gauss.Ey < 0);									// Clamps component flags to ensure amplitudes.
    gauss.Ez = (gauss.Ez > 0) - (gauss.Ez < 0);

    gauss.i = gauss.x/h1 + 0.5;

    const double I = 0.5*(gauss.Ey + gauss.Ez)*
                     gauss.ampl*gauss.ampl*units (mc_E0)*units (mc_E0)*mc_CGS_c/
                     (4*mc_pi)/mc_CGS_wattPerSquareCantimeter;
    say ("  gaussSpot_init: ");
    say ("  - eE_ampl/mc\\omega = %e, I = %.4e [W/cm^2]",
                gauss.ampl*units (mc_E0)/units (mc_A0), I);
    say ("  - P = %e [W/cm] = %e [W]",
                I*gauss.width*sqrt (mc_pi/2.0)*units (mc_r0),
                I*gauss.width*gauss.width*mc_pi/2.0*units (mc_r0)*units (mc_r0));
    say ("  - polarisation ((Ey:%.1f, Ez:%.1f)), omega/omega_0 = %.4f,",
                gauss.Ey, gauss.Ez, gauss.frequency/(2*mc_pi));
    say ("  - position (%f, %f, %f) [micron],",
                gauss.i*h1/units (mc_micron),
                gauss.y   /units (mc_micron),
                gauss.z   /units (mc_micron));
    say ("  - width (intensity FWHM) %f [laser lambda] = %f [micron],",
                gauss.width*sqrt (2*log (2)),
                gauss.width*sqrt (2*log (2))/units (mc_micron));
    say ("  - duration (intensity FWHM) %f [laser periods] = %f [fs] = %f [microns].",
                (gauss.duration + gauss.t_plateau)*sqrt (2*log (2)),
                (gauss.duration + gauss.t_plateau)*sqrt (2*log (2))/units (mc_femtosecond),
                (gauss.duration + gauss.t_plateau)*sqrt (2*log (2))/units (mc_femtosecond)*mc_CGS_c*1e-15/1e-4);
    say ("  - plateau duration %f [laser periods] = %f [fs] <=> %f [microns].",
                gauss.t_plateau*sqrt (2*log (2)), gauss.t_plateau*sqrt (2*log (2))/units (mc_femtosecond),
                gauss.t_plateau*sqrt (2*log (2))/units (mc_femtosecond)*mc_CGS_c*1e-15/1e-4);
    say ("  - launch delay %f [laser periods] = %f [fs] <=> %f [microns].",
                fabs (gauss.frontOffset), fabs (gauss.frontOffset)/units (mc_femtosecond),
                fabs (gauss.frontOffset)/units (mc_femtosecond)*mc_CGS_c*1e-15/1e-4);
    say ("  - longitudinal profile type is '%s'.",
                (gauss.t_plateau > 0.01*tau) ? "gaussian + const"
                                                   : "gaussian");

    gauss.Ey *= gauss.ampl;
    gauss.Ez *= gauss.ampl;
}

// ---------------------------------------------------------------------------
/// Tool to excplicitly shutdown all gauss-spot activity (overrides everything).
// ---------------------------------------------------------------------------
void
gaussSpot_deactivateAll (void)
{
    gaussSetupCheckOnly = 1;
}

// ---------------------------------------------------------------------------
/// Fills field in a focal plane with prescribed values regardless to everything
/// else (hard source).
// ---------------------------------------------------------------------------
void
gaussSpot_forceE (meshVec_p E, const double t)
{
    if (mf_mesh_pointIsOutside (E, gauss.i, E->jmin, E->kmin) || !mc_have_x ||
        gaussSetupCheckOnly)
        return;

    double timeFactor = 1;
    if (t < gauss.frontOffset) {
        double l = (t - gauss.frontOffset)/gauss.duration;
        timeFactor = exp (- l*l);
    }

    if (t > gauss.frontOffset + gauss.t_plateau) {
        double l = (t - gauss.frontOffset - gauss.t_plateau)/gauss.duration;
        timeFactor = exp (- l*l);
    }

    for (int j = E->jmin ; j <= E->jmax ; j++) {
        const double dy_ = mc_have_y*((j - 0.5)*h2 - gauss.y)/gauss.width,
                     dy  = mc_have_y*( j*h2        - gauss.y)/gauss.width;
        for (int k = E->kmin ; k <= E->kmax ; k++) {
            const double dz  = mc_have_z*( k*h3        - gauss.z)/gauss.width,
                         dz_ = mc_have_z*((k - 0.5)*h3 - gauss.z)/gauss.width;
            mv_fy(E, gauss.i, j, k) = gauss.Ey*sin (gauss.frequency*t)*timeFactor*exp (- dy_*dy_ - dz*dz);
            mv_fz(E, gauss.i, j, k) = gauss.Ez*cos (gauss.frequency*t)*timeFactor*exp (- dy*dy   - dz_*dz_);
        }
    }
}
