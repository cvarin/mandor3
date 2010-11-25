/** \file tag_photoelectrons.h
  * Prototypes of interfaces of 'tag_photoelectrons.c' module.
  *
  * Distribution function (DF) is
  *   \f$ \displaystyle f(\vec v) = \frac{3n_\alpha}{4\pi}\cdot\frac{v_z^2}{V_0^4}\cdot\delta (v - V_0) \f$.
  *
  * Sampling of the DF is done using quiet start (QS) method of Denavit and Wallsh.
  * Zero mean velocity and axial symmetry may be achieved with numerical accuracy by
  * sampling only subdomain of velocity space and forcing symmetry explicitly to
  * define the rest of the DF.
  *
  * Example of the config file entry is
    <pre>
    [photoelectrons]
    @ 1 	number of particles per mesh step along X.
    @ 1 	number of particles per mesh step along Y.
    @ 20 	number of particles per mesh step along Z.
    @ 20	number of particles in each spatial point.
    @ 1		mirror symmetrizator request (0 - no mirror).
    @ 1		axial symmetrizator request (nRotations, 0 - no rotational symmetry).
    @ 1.0e0	charge density parameter Q: positive => means (omega_pe/omega_0)^2, negative => density [cm^-3]
    @ -1.0	q/M [e/m]: charge to mass ratio.
    @ 0.1	thermal parameter V: positive => V0 [c], negative => energy [eV].
    @ 1		staggered mesh steps along X.
    @ 1		staggered mesh steps along Y.
    @ 1		staggered mesh steps along Z.
    @ 0		sampling: 1 to have uniform weight of markers, other to have uniform spacing of markers.
    </pre>
  *
  * Charge density parameter \b Q is interpreted as follow:
  * \f[ \left\{
  *   \begin{array}{ll}
  *   \displaystyle   Q > 0   &   \mbox{indirectly sets charge density so } (\omega_{pe}/\omega_0)^2 = n/n_{critical} = Q \\[6mm]
  *   \displaystyle   Q < 0   &   \mbox{Q works as } n_\alpha \mbox{ in } [cm^{-3}]: \rho = \mathrm{sign}(q/M)\cdot|Q|
  *   \end{array}
  * \right. \f]
  *
  * Photo-ionization energy (thermal) parameter \b V works like this:
  * \f[ \left\{
  *   \begin{array}{ll}
  *   \displaystyle   V > 0   &   \mbox{sets speed directly: }V_0 = V \\[6mm]
  *   \displaystyle   V < 0   &   \mbox{works like energy in eV: } \frac {m_eV_0^2}2 = |V|
  *   \end{array}
  * \right. \f]
  *
  * \note Uniform weight of markers seems more natural (distribution of markers mimics
  * distribution of real electrons), but uniform sampling of markers over surface of sphere
  * (weight changes to reproduce DF) helps to resolve regions with small concentration of
  * particles (possibly helping to uniformly resolve derivatives in the velocity space) and
  * gives MUCH better agreement in the growth rate.
  */

#ifndef MC_TAG_PHOTOELECTRONS_HEADER
#define MC_TAG_PHOTOELECTRONS_HEADER

#include <stdio.h>

#include "type_mesh.h"

/// Parameters ID for the tag_photoDF_parameters(int param, void *pntr).
#define mc_photoDF_V0		0	///< Signals to tag_photoDF_parameters(int, void *pntr) to write V<sub>0</sub> to the pointed location.
#define mc_photoDF_omega2_pe	1	///< Signals to tag_photoDF_parameters(int, void *pntr) to write \f$ \omega_{pe}^2 \f$ to the pointed location.
#define mc_photoDF_ID		3	///< Signals to tag_photoDF_parameters(int, void *pntr) to write chapter-ID to the pointed location.

double tag_photoelectrons (FILE *fp);
void tag_photoDF_parameters (int param, void *pntr);	// XXX - remove, use global parameters.

#endif
