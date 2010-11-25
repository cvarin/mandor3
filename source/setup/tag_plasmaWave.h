/** \file tag_plasmaWave.h
  * \brief Adds perturbation to the mean velocity of particles.
  *
  * Perturbation is in form of \f$ \displaystyle \vec v_j \to \vec v_j + A\cdot \frac qM \cdot \cos (\vec k\cdot \vec r_j), \f$
  * where \f$ \displaystyle \vec k = \left(\frac {2\pi m_x}{L_x}, \frac {2\pi m_y}{L_y}, \frac {2\pi m_z}{L_z}\right) \f$.
  *
  * Example of the config file entry is
    <pre>
    [plasmaWave]
    @ 1 	wave number along X.
    @ 0 	wave number along Y.
    @ 1 	wave number along Z.
    @ 0.1e-3	Amplitude of the perturbation [c].
    </pre>
  *
  */

#ifndef MC_TAG_PLASMA_WAVE_HEADER
#define MC_TAG_PLASMA_WAVE_HEADER

#include <stdio.h>

void tag_plasmaWave (FILE *fp);

#endif
