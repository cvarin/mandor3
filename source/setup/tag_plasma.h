/** \file tag_plasma.h
  * Plasma storage and access subsystem.
  *
  * tag_plasma() samples piece of plasma in space and postpones the velocity
  * distribution to other tags, while doing all shaping/partitioning for them.
  * Shaping is done using filters and all objects like foils, spikes and so
  * on are shaped here in the [plasma] tag. During a creation all markers are
  * enumerated from zero to 'particles per cell' and this number is stored in
  * marker_t::vx field. All routines of quiet start may use this number to
  * assign the velocity properly (see tag_maxwell(), for example).
  *
  * Plasma shaping is done by removing the pieces. If particle fits one of
  * the blocks it is removed. Blocks are polyhedrons (limited by set of planes,
  * normals point to the plasma-free half-spaces), spherical or cylindrical
  * domains (positive radius means we want to keep internal regions and negative
  * means we want outside regions). Option/parameter interleaving helps to group
  * together all blocks.
  *
  * Example of the config file entry is
  * <pre>
  * [plasma]
  * @ 10		nx - number of samples per cell along X axis.
  * @ 5			ny - number of samples per cell along Y axis.
  * @ 5			nz - number of samples per cell along Z axis.
  * @ 4			Number of points per spatial point (velocity layers).
  * @ planes		All optional arguments below set cut-off planes.
  * > 0.2 0.2 0.2	r0		block #0 / plane #0
  * > 1 -1 0		n0
  * > 1			shift0
  * > 0.2 0.2 0.2	r0		block #0 / plane #1
  * > 1 -1 0		n0
  * > 1			shift0
  * @ spheres		Adds or removes particles in the sphere.
  * > -1 1 0 0.5	x0, y0, z0 and r ('r >= 0' - add particles, 'r < 0' - remove).
  * @ cylinders		Adds or removes particles in the cylinder.
  * > -1 1 0		r0
  * > 1 1 0		n
  * > 0.5		Radius ('r >= 0' - add particles, 'r < 0' - remove).
  *          ...
  * </pre>
  */

#ifndef MC_TAG_PLASMA_HEADER
#define MC_TAG_PLASMA_HEADER

#include <stdio.h>

// Parameters of the last chunk created by tag_plasma().
extern const int plasma_nx, plasma_ny, plasma_nz, plasma_layers, plasma_PPC;

double    tag_plasma (FILE *fp);

#endif
