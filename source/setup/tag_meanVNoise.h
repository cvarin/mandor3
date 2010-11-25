/** \file tag_meanVNoise.h
  * \brief Creates noise perturbation of the mean velocity. Used as cheap way to introduce initial electrostatic (ES) perturbation when
  * exact expression for corresponding ES wave is complicated or unknown (typically it is used on the first stages of the new plasma
  * distribution function research).
  */

#ifndef tag_meanVNoise_header
#define tag_meanVNoise_header					///< \internal Guard against multiple inclusions.

#include <stdio.h>

void tag_meanVNoise (FILE *fp);

#endif
