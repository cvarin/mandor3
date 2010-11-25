/** \file tag_units.c
  * Sets reference scale to generate all other units.
  *
  * For now the basic scale is wavelength or critical density. All other units
  * are evaluated using speed of light and parameters of the electron (charge
  * and mass). All details are in misc_units.c. The main benefit is that
  * all dimensionalizing can be done straightforwardly in any directions.
  *
  * \warning All simulations in the main module are \b always carried on in
  *          dimensionless units.
  *
  * Example of the config file entry is
    <pre>
    [units]
    @ 1.0       Positive => works as wavelengh in microns [10^-4 cm].
                Negative => works as critical electrons' density [cm^-3].
    </pre>
  *
  * \note Setting wavelengh is convinient for EM simulations (laser-plasma
  * interaction, experimentalists usually use microns) while direct setting
  * of the critical charge denstity is really convinient for elestrostatic/EM
  * simulations of the initial value problems (instabilities) in plasma due to
  * automatic checks of the corresponding wavelength resolution.
  */

#include <stdlib.h>

#include "misc_units.h"
#include "log.h"
#include "misc_cfgReader.h"

// ---------------------------------------------------------------------------
/// Gets reference parameter and passes it to the scale evaluation routine.
// ---------------------------------------------------------------------------
void
tag_units (FILE *fp)
{
    double parameter = cfg_readDouble (fp);

    if (parameter > 0)
        units_setLambda (parameter*1e-4);
    else
        units_setCriticalDensity (- parameter);

    say ("tag_units:\n  o critical electrons concentration = %e cm^{-3}",
              units (mc_ne_critical));
    say ("  o lambda = %e [cm] = %e [micron].",
              units (mc_r0), units (mc_r0)/1.0e-4);
}
