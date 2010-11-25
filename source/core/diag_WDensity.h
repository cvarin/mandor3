/** \file diag_WDensity.h
  * Outputs averaged energies (thermal, electric field, magnetic field and so on) into text file (interfaces).
  */

#ifndef mc_diag_WDensity_header
#define mc_diag_WDensity_header					///< \internal Guard.

void wDensity_prepare (int startNew);
void wDensity_addPoint (double W_E, double W_M, double WTx, double WTy, double WTz);
void wDensity_flushData (void);

#endif
