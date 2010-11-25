/** \file misc_units.c
  * Units manager (see misc_units.h for details).
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"

static int initialized = 0;

// All simulation scales (in CGS units). To convert we simply multiply
// dimensionless parameter on the scale or divide dimensional one on the scale.
// Take a look on units_update.
static double units_lambda          = -1,
              units_omega           = -1,
              units_E               = -1,
              units_n               = -1,
              units_v               = -1,
              units_t               = -1,
              units_rho             = -1,
              units_micron          = -1,
              units_femtosecond     = -1,
              units_criticalElectronsConcentration = -1,
              units_plasmaPeriod    = -1,
              units_plasmaFrequency = -1,
              units_debayScale      = -1,
              units_S               = -1,
              units_A               = -1;

// ---------------------------------------------------------------------------
/// Saves units to the disk (normally is done only by 'setup.out' and only in
/// 'create' mode).
// ---------------------------------------------------------------------------
void
units_save (void)
{
   ENSURE (initialized, "units are not prepared to save them");

   FILE *fp = cfg_open ("binData/units.bin", "wt", __func__);
   fwrite (&units_lambda,          sizeof (double), 1, fp);
   fwrite (&units_v,               sizeof (double), 1, fp);
   fwrite (&units_t,               sizeof (double), 1, fp);
   fwrite (&units_omega,           sizeof (double), 1, fp);
   fwrite (&units_E,               sizeof (double), 1, fp);
   fwrite (&units_S,               sizeof (double), 1, fp);
   fwrite (&units_A,               sizeof (double), 1, fp);
   fwrite (&units_n,               sizeof (double), 1, fp);
   fwrite (&units_rho,             sizeof (double), 1, fp);
   fwrite (&units_criticalElectronsConcentration, sizeof (double), 1, fp);
   fwrite (&units_debayScale,      sizeof (double), 1, fp);
   fwrite (&units_micron,          sizeof (double), 1, fp);
   fwrite (&units_femtosecond,     sizeof (double), 1, fp);
   fwrite (&units_plasmaFrequency, sizeof (double), 1, fp);
   fwrite (&units_plasmaPeriod,    sizeof (double), 1, fp);
   fclose (fp);
}

// ---------------------------------------------------------------------------
/// Gets all units from disk (recorded by units_save).
// ---------------------------------------------------------------------------
void
units_load (void)
{
   ENSURE (!initialized, "double initialization");

   FILE *fp = cfg_open ("binData/units.bin", "rt", __func__);
   fread (&units_lambda,          sizeof (double), 1, fp);
   fread (&units_v,               sizeof (double), 1, fp);
   fread (&units_t,               sizeof (double), 1, fp);
   fread (&units_omega,           sizeof (double), 1, fp);
   fread (&units_E,               sizeof (double), 1, fp);
   fread (&units_S,               sizeof (double), 1, fp);
   fread (&units_A,               sizeof (double), 1, fp);
   fread (&units_n,               sizeof (double), 1, fp);
   fread (&units_rho,             sizeof (double), 1, fp);
   fread (&units_criticalElectronsConcentration, sizeof (double), 1, fp);
   fread (&units_debayScale,      sizeof (double), 1, fp);
   fread (&units_micron,          sizeof (double), 1, fp);
   fread (&units_femtosecond,     sizeof (double), 1, fp);
   fread (&units_plasmaFrequency, sizeof (double), 1, fp);
   fread (&units_plasmaPeriod,    sizeof (double), 1, fp);
   fclose (fp);

   initialized = 1;
}

// ---------------------------------------------------------------------------
/// Takes lambda and defines all other units. Caller is responsible for the
/// 'normalness' of the lambda (>= 0).
// ---------------------------------------------------------------------------
static void
units_update (double lambda)
{
   units_lambda = lambda;
   units_v      = mc_CGS_c;					// v = c
   units_t      = lambda/units_v;				// t = r0/v
   units_omega  = 2*mc_pi/units_t;				// ω = 2π/t0
   units_A      = mc_CGS_m*mc_CGS_c*units_omega/mc_CGS_e;	// A0 = mcω/e
   units_E      = units_A/(2*mc_pi);				// E0 = mcω/(2πe)
   units_n      = mc_CGS_m*units_omega*units_omega/
                  (16*mc_pi*mc_pi*mc_pi*mc_CGS_e*mc_CGS_e);	// n = mω²/(16π³e²)
   units_rho    = mc_CGS_e*units_n;				// ρ = en

   units_micron          = 1.0e-4/lambda;			// one micron in dimensionless units
   units_femtosecond     = 1.0e-15/units_t;			// one fs in dimensionless units
   units_plasmaFrequency = 2*mc_CGS_e*sqrt (mc_pi*units_n/mc_CGS_m)*units_t;	// sqrt (4πne²/m)⋅t0.
   units_plasmaPeriod    = 2*mc_pi/units_plasmaFrequency;	// 2π/(plasmaFrequency⋅t0) = 2π/units_plasmaFrequency.
   // n_cr = ω²m/(4πe²), which means ω_pe(n_cr) = ω_0.
   units_criticalElectronsConcentration = units_omega*units_omega*mc_CGS_m/(4*mc_pi*mc_CGS_e*mc_CGS_e);
}

// ---------------------------------------------------------------------------
/// Returns unit requested by id defined in header.
// ---------------------------------------------------------------------------
double
units (int id)
{
   // Packed parameters to replace switch. Order of parameters is defined by
   // header defines.
   double values[mc_unitsN] = {
      units_n, units_lambda, units_t, units_E, units_S, units_rho,
      units_criticalElectronsConcentration, units_plasmaPeriod,
      units_debayScale, units_micron, units_femtosecond, units_A, units_v
   };

   ENSURE (initialized, "units are not set");
   ENSURE (id >= 0 && id < mc_unitsN, "bad unit id: %d", id);
   ENSURE (values[id] > 0, "unsupported unit (value = %le)", values[id]);

   return values[id];
}

// ---------------------------------------------------------------------------
/// All units are set for this wavelength (laser-plasma interaction class of
/// problems).
// ---------------------------------------------------------------------------
void
units_setLambda (double lambda)
{
   ENSURE (!initialized, "double initialization");
   ENSURE (lambda > 0,   "bad lambda (%e)", lambda);

   units_update (lambda);
   initialized = 1;
}

// ---------------------------------------------------------------------------
/// All units are set to force this critical concentration of electrons
/// (electrostatic turbulence class of problems).
// ---------------------------------------------------------------------------
void
units_setCriticalDensity (double n_cr)
{
   ENSURE (!initialized, "double initialization");
   ENSURE (n_cr > 0,     "bad ne_cr (%e)", n_cr);

   // ω0     = ω_pe (n_cr)
   // lambda = 2πc/ω0
   units_omega = sqrt (4*mc_pi*n_cr*mc_CGS_e*mc_CGS_e / mc_CGS_m);
   units_update (2*mc_pi*mc_CGS_c/units_omega);

   initialized = 1;
}
