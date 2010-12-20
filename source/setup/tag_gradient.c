/** \file tag_gradient.c
  * Apply density gradient to plasma (see tag_gradient.h)
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "setup/plasma.h"
#include "setup/tag_plasma.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_partition.h"
#include "misc_parameters.h"

#include "main.h"

void tag_gradient(FILE *fp) {
  double x1, y1, z1, x2, y2, z2, d1, d2;
  int type;

  ENSURE(3 == fscanf(fp, "> %le %le %le %*[^\n] ", &x1, &y1, &z1), "Cannot read point 1");
  ENSURE(3 == fscanf(fp, "> %le %le %le %*[^\n] ", &x2, &y2, &z2), "Cannot read point 2");
  ENSURE(3 == fscanf(fp, "> %le %le %d %*[^\n] ", &d1, &d2, &type), "Cannot read point 3");

  double lx = x2 - x1,
         ly = y2 - y1,
         lz = z2 - z1;

  double width = lx*lx + ly*ly + lz*lz;

  // In fact we apply not density, but charge gradient.
  // It seems to me that it is just the same
  long int N;
  for (marker_t *p = plasma_getObject(0, &N), *end = p + N; p < end; ++p) { 
    double distance = ((p->x - x1)*lx + (p->y - y1)*ly + (p->z - z1)*lz) / width;
    if ((distance > 0) && (distance < 1)) {
      // we're inside the gradient region
      if (type == 0)
        p->rho *= d1 + (d2-d1) * distance; // linear
      else
        p->rho *= d1 * pow(d2 / d1, distance); // exponential 
    }
  }

  say("%s: density gradient is applied", __func__);
  say("  - points (%e, %e, %e)-(%e, %e, %e)", x1, y1, z1, x2, y2, z2);
  say("  - normal (%e, %e, %e)", lx, ly, lz);
  say("  - width = %e", width);
}

