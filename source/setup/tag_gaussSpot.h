/** \file tag_gaussSpot.h
  * \brief Hard source (gaussian/half-gaussian pulses) module.
  *
  * <b>Hard source</b> sets electric field directly. This type of sources is not
  * suitable for general modelling because of the plane of this source reflects
  * all incoming radiation like a mirror (see core/em_TFSF.h on soft sources).
  * Nevertheless it may be used as a generator of Maxwell equations' solution to
  * be used with TF/SF interface.
  *
  * Example of the config file entry is
    <pre>
    [gaussSpot]
    @ 6.25e+22      Energy parameter: positive => I [W/cm<sup>2</sup>] / negative => e*E/(m*c*omega).
    @ 1.0           Cyclic frequency [omega_0]. Must be 1 (or you understand what do you do).
    @ 14.0          Focus X coordinate [micron].
    @ 5.0           Focus Y coordinate [micron].
    @ 5.0           Focus Z coordinate [micron].
    @ 3.0           Gauss-spot width (positive => FWHM for laser intensity [micron], negative => gauss envelope for field amplitude [r0]).
    @ 30.0          Pulse duration (positive => FWHM for laser intensity [fs], negative => gauss envelope for field amplitude [t0]).
    @ 40.0          Pulse leading front offset (leading edge to pulse max delay) [fs].
    > 10.0          Plato duration [fs].
    @ 0             EY != 0
    @ 1             EZ != 0
    > 20.0          TF/SF emitter face position [micron].
    </pre>
  *
  * Notes on parameters converstion:
  *
  * - 'Peak intensity' to 'electric field amplitude':
  *   - \f$ I = c/(4\pi)[E,H]\cdot<..> \f$,
  *   - averaging over time results in factor
  *     \f$ \displaystyle <..> = \left\{ \begin{array}{ll} \displaystyle 0.5\quad & \mbox{linear polarisation} \\
  *                                                        \displaystyle 1.0 & \mbox{circular polarisation} \end{array} \right.\f$,
  *   - \f$ E = \sqrt {4\pi I/c<..>} \f$.
  *
  * - 'FWHM of intensity' to 'exponential gaussian factor for field':
  *   - \b FWHM means Full Width at Half Magnitude, input parameters are
  *     \f$ T_{FWHM},\ \Delta_{FWHM} \f$;
  *   - electric field in focal plane is
  *     \f[ \displaystyle E (t, r_\perp) = E_0\cdot e^{\displaystyle - \frac{(t - t_0)^2}{DT^2} - \frac{(r_\perp - r_0)^2}{DR^2}}; \f]
  *   - by definition of FWHM we have
  *     \f[\displaystyle \frac {E^2(t_0 - 0.5\cdot T_{FWHM})}{E^2(t_0)} = e^{\displaystyle - 2\cdot\frac{T_{FWHM}^2}{DT^2}} \equiv \frac 12 \f]
  *   - \f$ DT = T_{FWHM}/\sqrt {2\ln(2)} \f$, \f$ DR = \Delta_{FWHM}/\sqrt {2\ln(2)} \f$.
  *
  * - Gaussian parameters to the power hint:
  *   - \f$ \int_R \exp (-2x^2/\Delta_x^2) dx = \Delta\cdot\sqrt {\pi/2} \f$
  *   - \f$ P = I\cdot \int_{R^2} \exp (-2x^2/\Delta_x^2) \exp (-2y^2/\Delta_y^2) dx dy = I \Delta_x \Delta_y\cdot \pi/2\f$ for 3D
  *   - \f$ P = I\cdot \int_{R} \exp (-2x^2/\Delta_x^2) dx = I \Delta_x \cdot \sqrt{\pi/2}\f$ for 2D
  */

#ifndef MC_TAG_GAUSSSPOT_HEADER
#define MC_TAG_GAUSSSPOT_HEADER

#include <stdio.h>

void tag_gaussSpot (FILE *fp, const char *name);

#endif
