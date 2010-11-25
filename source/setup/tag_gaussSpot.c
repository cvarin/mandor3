/** \file tag_gaussSpot.c
  * \brief Hard source (gaussian/half-gaussian pulses) module, see also tag_gaussSpot.h.
  */

#include <math.h>
#include <stdlib.h>

#include "type_mesh.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"

static int gaussSpotIsPlaced = 0;										///< Flag to check uniqueness of the source.

// ---------------------------------------------------------------------------
/// Parameters of the gauss spot (packed into exchange structure).
// ---------------------------------------------------------------------------
struct gaussPack_s
{
  double ampl, frequency;
  double x, y, z;
  double width, duration, frontOffset, plateauDuration;
  double Ey, Ez;
  int i;
} gauss;

// ---------------------------------------------------------------------------
/// Reads written config file the way its done in modelling back-end and reports pulse parameters.
// ---------------------------------------------------------------------------
static void
tag_gaussSpot_report (void)
{
  if (cpu_here)													// Master reports everything.
    return;

  FILE  *fp = cfg_open (".gaussSpot.cfg", "rt", __func__);
  gauss.ampl = cfg_readDouble (fp);
  gauss.frequency = 2*mc_pi*cfg_readDouble (fp);
  gauss.x = cfg_readDouble (fp);
  gauss.y = cfg_readDouble (fp);
  gauss.z = cfg_readDouble (fp);
  gauss.width = cfg_readDouble (fp);
  gauss.duration = cfg_readDouble (fp);
  gauss.frontOffset = cfg_readDouble (fp);
  gauss.plateauDuration = cfg_readDouble (fp);
  gauss.Ey = cfg_readInt (fp);
  gauss.Ez = cfg_readInt (fp);
  fclose (fp);

  gauss.Ey = (gauss.Ey > 0) - (gauss.Ey < 0);									// Clamps component flags to ensure amplitudes.
  gauss.Ez = (gauss.Ez > 0) - (gauss.Ez < 0);

  gauss.i = gauss.x/h1 + 0.5;

  const double I = 0.5*(gauss.Ey + gauss.Ez)*gauss.ampl*gauss.ampl*units (mc_E0)*units (mc_E0)*mc_CGS_c/(4*mc_pi)/mc_CGS_wattPerSquareCantimeter;
  say ("  gaussSpot_init(draft launch using separate copy from em_gaussSpot.c): ");
  say ("  - eE_ampl/mc\\omega = %e, I = %.4e [W/cm^2]", gauss.ampl*units (mc_E0)/units (mc_A0), I);
  say ("  - P = %e [W/cm] = %e [W]", I*gauss.width*sqrt (mc_pi/2.0)*units (mc_r0), I*gauss.width*gauss.width*mc_pi/2.0*units (mc_r0)*units (mc_r0));
  say ("  - polarisation (Ey:%.1f, Ez:%.1f), omega/omega_0 = %.4f", gauss.Ey, gauss.Ez, gauss.frequency/(2*mc_pi));
  say ("  - position (%f, %f, %f) [micron]", gauss.i*h1/units (mc_micron), gauss.y/units (mc_micron), gauss.z/units (mc_micron));
  say ("  - width (intensity FWHM) %f [laser lambda] = %f [micron]", gauss.width*sqrt (2*log (2)), gauss.width*sqrt (2*log (2))/units (mc_micron));
  say ("  - duration (intensity FWHM) %f [laser periods] = %f [fs] <=> %f [microns].",
            (gauss.duration + gauss.plateauDuration)*sqrt (2*log (2)), (gauss.duration + gauss.plateauDuration)*sqrt (2*log (2))/units (mc_femtosecond),
            (gauss.duration + gauss.plateauDuration)*sqrt (2*log (2))/units (mc_femtosecond)*mc_CGS_c*1e-15/1e-4);
  say ("  - plateau duration %f [laser periods] = %f [fs] <=> %f [microns].",
            gauss.plateauDuration*sqrt (2*log (2)), gauss.plateauDuration*sqrt (2*log (2))/units (mc_femtosecond),
            gauss.plateauDuration*sqrt (2*log (2))/units (mc_femtosecond)*mc_CGS_c*1e-15/1e-4);
  say ("  - launch delay %f [laser periods] = %f [fs] <=> %f [microns]", fabs (gauss.frontOffset), fabs (gauss.frontOffset)/units (mc_femtosecond),
            fabs (gauss.frontOffset)/units (mc_femtosecond)*mc_CGS_c*1e-15/1e-4);
  say ("  - longitudinal profile type is '%s'", (gauss.plateauDuration > 0.01*tau) ? "gaussian + const" : "gaussian");
}

// ---------------------------------------------------------------------------
/// Processes all parameters and prepares file '<b>.gaussSpot.cfg</b>' for main simulation module.
// ---------------------------------------------------------------------------
void
tag_gaussSpot (FILE *fp, const char *name)
{
  double I = cfg_readDouble (fp);										// Peak intensity I [Watt/cm^2] or A0.
  double Omega = cfg_readDouble (fp);										// Frequency in omega_0.
  double X   = cfg_readDouble (fp)*units (mc_micron);								// Position of focal spot [microns].
  double Y   = cfg_readDouble (fp)*units (mc_micron);
  double Z   = cfg_readDouble (fp)*units (mc_micron);
  double R   = cfg_readDouble (fp);										// Gauss spot width
  double Tau = cfg_readDouble (fp);
  double frontOffset  = cfg_readDouble (fp)*units (mc_femtosecond);						// Cut-off delay [fs].
  double plateauDuration = 0;											// Duration fo the const part of pulse body.
  if (cfg_isOption (fp))
    plateauDuration = cfg_readOptDouble (fp)*units (mc_femtosecond);
  int EY  = cfg_readInt (fp);											// Polarization coefficients.
  int EZ  = cfg_readInt (fp);

  double TFSFpos = - 1.0e10;											// Optional TFSF position.
  if (cfg_isOption (fp))
    TFSFpos = cfg_readOptDouble (fp)*units (mc_micron);

  ENSURE (!gaussSpotIsPlaced, "only one gaussian hard-source, please");
  gaussSpotIsPlaced = 1;

  R = (R > 0) ? R*units (mc_micron)/sqrt (2.0*log (2.0)) : - R;							// Converts R to [r0].
  Tau = (Tau > 0) ? Tau*units (mc_femtosecond)/sqrt (2.0*log (2.0)) : - Tau;					// Converts duration to [t0].

  if (((X - 5*h1)*(Lx - 5*h1 - X) < 0 || !mc_have_x) || 							// Checks that spot is inside of the domain.
       (Y - 5*h2)*(Ly - 5*h2 - Y)*mc_have_y < 0 ||
       (Z - 5*h3)*(Lz - 5*h3 - Z)*mc_have_z < 0)
  {
    say ("tag_gaussSpot:\n  - spot position is (%e, %e, %e) [r0]", X, Y, Z);
    say ("  - spot mesh position is (%d, %d, %d)", (int) (X/h1), (int) (Y/h2), (int) (Z/h3));
    say ("  - domain size is (%d, %d, %d)", dmn_max[0], dmn_max[1], dmn_max[2]);
    DIE ("spot is badly placed (or X axis is deactivated)");
  }

  ENSURE (plateauDuration > - 1.0e-10*tau,
          "bad const-piece duration (%e)", plateauDuration);

  EY = (EY > 0) - (EY < 0);											// Clamps: coefficient = sign (coefficient).
  EZ = (EZ > 0) - (EZ < 0);
  ENSURE (EY || EZ, "degenerated polarisation (%d, %d)", EY, EZ);

  double Ampl = (I > 0) ? sqrt (4*mc_pi*I*mc_CGS_wattPerSquareCantimeter/((abs (EY) + abs(EZ))*0.5*mc_CGS_c))/units (mc_E0) : - I*units (mc_A0)/units (mc_E0);

  if (cpu_here)													// Master will update all files.
    return;

  FILE *fpCfg = cfg_open (".gaussSpot.cfg", "wt", __func__);							// Creates source dimensionless parameters file.
  fprintf (fpCfg, "@ %e       Amplitude of the source [E0].\n", Ampl);
  fprintf (fpCfg, "@ %e       Cyclic frequency [omega_0].\n", Omega);
  fprintf (fpCfg, "@ %e       Focus plane coordinate [r0].\n", X);
  fprintf (fpCfg, "@ %e       Focus Y coordinate [r0].\n", Y);
  fprintf (fpCfg, "@ %e       Focus Z coordinate [r0].\n", Z);
  fprintf (fpCfg, "@ %e       Gauss-spot size parameter [r0].\n", R);
  fprintf (fpCfg, "@ %e       Pulse duration [t0].\n", Tau);
  fprintf (fpCfg, "@ %e       Pulse front offset [t0].\n", frontOffset);
  fprintf (fpCfg, "@ %e       Pulse const-part duration [t0].\n", plateauDuration);
  fprintf (fpCfg, "@ %d       Polarization (Ey).\n", EY);
  fprintf (fpCfg, "@ %d       Polarization (Ez).\n", EZ);
  fprintf (fpCfg, "\n\n   This file is automatically created by setup routine using \"%s\" input file.\n", name);
  fclose (fpCfg);

  say ("tag_gaussSpot:");
  tag_gaussSpot_report ();

  say ("  - to pass 2 pulse delay times (%.4f [t0] = %.4f [fs]) you need %d steps", 2*frontOffset,
            2*frontOffset/units (mc_femtosecond), (int) (2*frontOffset/tau));

  if (TFSFpos > - 1.0e8)											// ASSUMES c = 1 IN DIMENSIONLESS UNITS.
    say ("  - to record pulse with TFSF recorder at %.4f [r0] = %.4f [micron] (I = %d) you need %d steps", TFSFpos, TFSFpos/units (mc_micron),
              (int) (TFSFpos/h1), (int)((2*frontOffset + plateauDuration + fabs (TFSFpos - X))/tau));
}
