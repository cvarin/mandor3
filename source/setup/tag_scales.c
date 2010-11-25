/** \file tag_scales.c
  * \brief Prints all key scales to hint possible limitations of PIC.
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "misc_units.h"

#include "log.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

enum {eSystems_laser, eSystems_selfHeating};				///< Enumerated types of the simulation cases.

// ---------------------------------------------------------------------------
/// Outputs dispersion properties of the Yee solver.
// ---------------------------------------------------------------------------
static void
tag_scalesLaser (FILE *fp)
{
  DIE ("not re-implemented yet (should show dispersion properties of the scheme)");
}

// ---------------------------------------------------------------------------
/// Outputs restrictions on the spatial and time steps due to solver limitations.
// ---------------------------------------------------------------------------
static void
tag_scalesSelfHeat (FILE *fp)
{
  double nMin = cfg_readOptDouble (fp); // Max (n_cr) to print data for.
  double nMax = cfg_readOptDouble (fp); // Max (n_cr) to print data for.
  int N = cfg_readOptDouble (fp);       // Number of levels.

  ENSURE (N > 0, "bad number of sublevels (N = %d)", N);

  const double mc2 = mc_CGS_m*mc_CGS_c*mc_CGS_c/mc_CGS_eV;	// mc2 [eV]
  say ("tag_scales:\n  self-heating scales/restrictions for e^- (mc^2 = %.2e [eV]) are", mc2);
  say ("   n/n_cr   skin(um)  r_D(1eV) r_D(10eV) r_D(100eV) r_D(1KeV) T_D(eV)");
  double h = h1;
  for (double n = nMin ; n <= nMax ; n += (nMax - nMin)/(double)(N))
  {
    const double r0 = units (mc_r0)/1.0e-4; 	// r0 in microns.
    const double rp = r0/sqrt (n);          	// r0*w_0/w_p.
    double skin = sqrt (2.0)*rp/(2*mc_pi);  	// Skin depth
    double r_D = rp/(2*mc_pi)*sqrt (1.0/mc2);	// r_D | T = 1 eV.

    // Computes self-heating stop temperatures.
    double T1 = mc2 * 4.0*mc_pi*mc_pi * n * mc_have_x*h1*h1;
    double T2 = mc2 * 4.0*mc_pi*mc_pi * n * mc_have_y*h2*h2;
    double T3 = mc2 * 4.0*mc_pi*mc_pi * n * mc_have_z*h3*h3;
    if (T1 < T2)
    {
      T1 = T2;
      h = h2;
    }
    if (T1 < T3)
    {
      T1 = T3;
      h = h3;
    }
    say ("  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e   %.2e",
         n, skin, r_D, r_D*sqrt (1e1), r_D*sqrt (1e2), r_D*sqrt (1e3), T1);
  }
  say ("    Debay temperature is computed using h = %.2e [r0] = %.2e [micron]",
       h, h/units (mc_micron));
}

// ---------------------------------------------------------------------------
/// Reads range of densities and outputs .
// ---------------------------------------------------------------------------
void
tag_scales (FILE *fp)
{
  const char *word = cfg_readWord (fp);
  const int ourCase = cfg_identifyWord (word, "laser", eSystems_laser, "self-heating", eSystems_selfHeating, mc_cfgTermGuesses);

  switch (ourCase)
  {
    case eSystems_laser:
      tag_scalesLaser (fp);
    break;

    case eSystems_selfHeating:
      tag_scalesSelfHeat (fp);
    break;

    default:
      DIE ("estimate mode '%s' is not installed", word);
    break;
  }
}
