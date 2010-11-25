/** \file tag_foil.h
  * Alternative sampling of the fully ionized cold foil.
  *
  * This functions creates plasma target (foil) as follow:
  * - plasma is fully preionized and cold
  * - profile can be linear (i.e. constant) or exponential
  * - orientation is set by center and normal
  * - only width of the foil is set (use filters to shape foil)
  * - 'shift' is used to position multiple layers using a common center
  *   as reference point (for user convinience)
  *
  * Example of the config file entry is
    <pre>
    [foil]
    @ 30.0              X coordinate of the center [micron].
    @ 20.0              Y coordinate of the center [micron].
    @ 20.0              Z coordinate of the center [micron].
    @ 0.96592583        x component of the foil normal.
    @ 0.25881905        y component of the foil normal.
    @ 0.0               z component of the foil normal.
    @ 2.0               Width of the foil [micron].
    > 2.0               Shift of the center along the normal, optional [micron].
    @ linear            Type of the profile (linear or exponential)
    @ -0.8e+0           Start concentration (positive => [cm^-3];  negative => [n_cr of electrons]).
    @ -2.0e+0           End   concentration (positive => [cm^-3];  negative => [n_cr of electrons]).
    @ -1.0              Charge of particle [e].
    @ +1.0              Mass of particle [m].
    @ 5                 Nx: number of particles per cell along X.
    @ 5                 Ny: number of particles per cell along Y.
    @ 1                 Nz: number of particles per cell along Z.
    </pre>
  *
  */

#ifndef MC_TAG_FOIL_HEADER
#define MC_TAG_FOIL_HEADER		///< \internal Guard.

#include <stdio.h>

double tag_foil (FILE *fp);

#endif
