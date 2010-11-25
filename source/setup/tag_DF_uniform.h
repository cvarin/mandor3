/** \file tag_DF_uniform.h
  * Sets uniform velocity distribution.
  *
  * This distribution function is useful to test boundary conditons and charge conservation
  * during development. Also it is very simple way to create cold background or cold beam.
  *
  * Example of the config file entry is
  * <pre>
  * [DF:uniform]
  * @ 5.0	rho	Charge density in ρ_critical for electrons.
  * @ -1.0	q/m	Charge to mass ratio.
  * @ 0.0	ux	Mean velocity.
  * @ 0.0	uy
  * @ 0.0	uz
  * @ 2.0	dUx	Width of uniform distribution.
  * @ 2.0	dUy
  * @ 2.0	dUz
  * </pre>
  */

#ifndef MC_TAG_DF_UNIFORM_HEADER
#define MC_TAG_DF_UNIFORM_HEADER

#include <stdio.h>

void tag_DF_uniform (FILE *fp);

#endif
