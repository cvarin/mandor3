/** \file tag_scissors.h
  * Scissor tag for plasma created.
  *
  * Main goal is to create plasma in given subregion of the space without adding parameters to creator tag (scissor type of operation).
  *
  * Example of the config file entry is
    <pre>
    [scissors]
    @ -1e100        cache filter xMin [micron].
    @ 10            cache filter xMax [micron].
    @ -1e100        cache filter yMin [micron].
    @ 10            cache filter yMax [micron].
    @ -1e100        cache filter zMin [micron].
    @ 10            cache filter zMax [micron].

      Cuts off plasma outside the box [<domain left-bottom-rear corner> / (10, 10, 10)].
    </pre>
  */

#ifndef MC_TAG_CACHEFILTER_HEADER
#define MC_TAG_CACHEFILTER_HEADER	///< Multiple include guard.

#include <stdio.h>

void tag_scissors (FILE *fp);

#endif
