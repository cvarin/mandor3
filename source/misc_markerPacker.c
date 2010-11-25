/** \file misc_markerPacker.c
  * \brief Set of routines to pack/unpack markers (may be of different types) into continuous buffer (mainly for || exchange).
  *
  * Main structure is just a wrapper around pointer to piece of memory to support continuous packing process. All stuff is
  * \b inline to help compiler to optimize (hopefully inline) all this routines.
  *
  * \todo Checks excessive calls to finalize/unit during all2all sessions - optimize order of calls to remove it.
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "log.h"
#include "misc_markerPacker.h"

