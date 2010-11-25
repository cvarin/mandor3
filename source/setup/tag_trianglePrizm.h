/** \file tag_trianglePrizm.h
  * Plasma target allocation. Target is a triangle prizm made of cold uniform plasma.
  * Example of input file entry:
  * <pre>
[trianglePrizm]
@ 2.0		x - coordinate of the first vertex, micron.
@ 2.0		y - coordinate of the first vertex, micron.
@ 4.0		x - coordinate of the second vertex, micron.
@ 2.0		y - coordinate of the second vertex, micron.
@ 2.0		x - coordinate of the third vertex, micron.
@ 4.0		y - coordinate of the third vertex, micron.
@ 1.0		z - coordinate of the lower end of the prizm, micron.
@ 2.0		z - coordinate of the upper end of the prizm, micron.
@ -100		concentration (positive => [cm^-3];  negative => [n_cr for electrons]).
@ -1.0		charge of particle [e].
@ +1.0		mass of particle [m].
@ 4		nx: number of particles per cell along X.
@ 4		ny: number of particles per cell along Y.
@ 1		nz: number of particles per cell along Z.
  * </pre>
  */

#ifndef MC_TAG_TRIANGLEPRIZM_HEADER
#define MC_TAG_TRIANGLEPRIZM_HEADER	///< Guard

#include <stdio.h>

double tag_trianglePrizm (FILE *fp);

#endif
