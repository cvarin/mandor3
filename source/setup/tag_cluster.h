/** \file tag_cluster.h
  * Cold uniform charge-neutral cluster creation interface (see tag_cluster.c).
  *
  * This functions creates plasma sphere (preionized and cold).
  *
  * Example of the config file entry is
    <pre>
    [cluster]
    @ 10            x-coordinate of the center [micron]
    @ 12            y-coordinate of the center [micron]
    @ 10            z-coordinate of the center [micron]
    @ 5             Radius [micron]
    @ -2.0          Density (positive => concentration in [cm^-3], negative => charge density in critical density for THIS specie).
    @ -1            Charge of particle [|e|]
    @ +1            Mass of particle [m_e]
    @ 5             Markers per cell per X direction
    @ 5             Markers per cell per Y direction
    @ 5             Markers per cell per Z direction

      Cold cluster with infinitely heavy ions.
    </pre>
  */

#ifndef mc_tag_cluster_header
#define mc_tag_cluster_header			///< \internal Guard.

#include <stdio.h>

double tag_cluster (FILE *fp);

#endif
