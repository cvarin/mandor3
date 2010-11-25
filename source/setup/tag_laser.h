/** \file tag_laser.h
  * Tag to add laser focused by parabolic mirror.
  *
  * Initial variant is written by Kostya Popov in his forked version of Mandor.
  * This is a soft source which is precomputed and cached by setup to be used
  * in main computational cycle directly from the cache.
  *
  * So here we just call interfaces. You should read 'em_mirror.h' for the
  * exact structure of '[focused laser]' tag.
  *
  * \warning The most important thing is that the source is computed using
  *          Stratton-Chu integrals and for now it is correct only for 3D.
  */

#ifndef MC_TAG_LASER
#define MC_TAG_LASER

#include <stdio.h>

double tag_laser (FILE *fp);

#endif
