/** \file tag_scales.h
  * \brief Prints all basic scales of the run and hints possible limitations of PIC.
  *
  * For now only laser-plasma interaction is under intensive development, so scales are printed for
  * numerical dispersion (laser) and plasma self-heating. Later other general types of simulations
  * might be considered as well.
  *
  * <h3>Self-heating (spatial) scales.</h3>
  *
  * In order to describe properly the evolution of the system we have to resolve all spatial scale
  * by our mesh. In decreasing order scales are:
  * - EM-wave wavelength: \f$ \lambda = 2\pi c/\omega_0 = r_0\f$ and \f$r_0\f$ is taken as a unit
  *   of size.
  * - Skin depth (\f$\delta\f$):\n
  *   According to Landau book skin depth \f$\delta\f$ is \f[\delta = c/\sqrt{2\pi\sigma\omega}\f]
  *   Dropping down temperature dispersion we have \f[\epsilon = 1 - \omega_p^2/\omega^2 = 1 + 4\pi i \sigma/\omega\f]
  *   It gives an estimate of the \f$\sigma\f$ \f[\sigma = i \omega_p^2/(4\pi\omega)\f]
  *   Finally we have \f[\delta \simeq \sqrt 2 c / \omega_p = \frac {\sqrt 2}{2\pi}\cdot \frac {\omega_0}{\omega_p}\cdot r_0\f]
  * - Debay scale (\f$r_D\f$):\n
  *   That is the cruelest restriction, because of PIC codes experience numerical instability if
  *   spatial step is bigger than \f$r_D\f$. From simple boundary case \f[r_D = \sqrt{\frac{T}{4\pi ne^2}} = h\f]
  *   we derive a low limit on the self-heating temperature (temperature of saturation):
  *   \f[T = 4\pi^2\cdot mc^2\cdot \left(\frac h{r_0}\cdot\frac{\omega_p}{\omega_0}\right)^2\f]
  */

#ifndef MC_TAG_SCALES_HEADER

/** \internal Guard. */
#define MC_TAG_SCALES_HEADER

void tag_scales (FILE *fp);

#endif
