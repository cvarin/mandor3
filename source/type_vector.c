/** \file type_vector.c
  * Helpers to work with 3-component vectors (see 'type_vector.h' for details).
  */

#include <stdio.h>

#include "type_vector.h"

// ----------------------------------------------------------------------------
/// Prints double vector in form of (%.3e, %.3e, %.3e) into internal buffer and
/// returns the pointer (up to 16 slots can be used simultaneously).
// ----------------------------------------------------------------------------
const char *
vd_print (vec3d_t v)
{
   static char buffers[16][100];
   static int  num = -1;
   num = (num + 1) & 15;
   sprintf (buffers[num], "(%.3e, %.3e, %.3e)", v.v.x, v.v.y, v.v.z);
   return buffers[num];
}

