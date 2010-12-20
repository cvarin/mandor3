/** \file tag_gradient.h
  * Applies density gradient to plasma
  * Example of config file:
  * <pre>
  * [gradient]
  * @ x1 y1 z1 starting point of the gradient
  * @ x2 y2 z2 end point
  * @ d1 d2 type d1 and d2 - start and end densities, type - 0 for linear, 1 for exponential
  * </pre>
  */


#ifndef MC_TAG_GRADIENT_HEADER
#define MC_TAG_GRADIENT_HEADER

#include <stdio.h>

void tag_gradient(FILE *fp);
#endif
