/** \file misc_units.h
  * Units manager. Main purpose is to provide simple way to convert to/from
  * dimensionless units. It simplifies diagnostic output and intermidiate
  * setup computations which are much more transparent if done in CGS rather
  * than in dimensionless form.
  */

#ifndef mc_n0

/*
 * Dimensional units' ids for different plasma/experiment parameters
 * (parameters of the call to 'units (id)').
 *
 * IMPORTANT: if changed you MUST fix 'double units (int id)' function !
 */
#define mc_n0					0
#define mc_r0					1
#define mc_t0					2
#define mc_E0					3
#define mc_S0					4
#define mc_rho0					5
#define mc_ne_critical				6
#define mc_plasmaPeriod				7
#define mc_debayScale				8
#define mc_micron				9
#define mc_femtosecond				10
#define mc_A0					11
#define mc_v0					12
#define mc_unitsN				13

/*
 * Physical units are taken from "NRL Plasma Formulary".
 */
#define mc_CGS_c			2.9979e+10	///< [cm/sec], speed of light in vacuum.
#define mc_CGS_e			4.8032e-10	///< [statcoulomb (statcoul)], absolute value of the charge of electron.
#define mc_CGS_m			9.1094e-28	///< [gr], mass of the electron.
#define mc_CGS_eV			1.6022e-12	///< [erg], energy assiciated with one electron-volt.
#define mc_CGS_wattPerSquareCantimeter	1e7		///< [erg/(sec*cm^2)].
#define mc_pi		3.141592653589793238462643	///< [1], this value of \pi is used in MPI test routine as an exact result.

double units                    (int id);
void   units_setLambda          (double lambda);
void   units_setCriticalDensity (double n_cr);
void   units_save               (void);
void   units_load               (void);

#endif
