/** \file misc_markerPlacer.h
  * \brief This particle sorter is used to restore distribution of particles over subdomains (global, relatively expensive,
  * no limitations of partitioning).
  *
  * This library is designed to do \e initialization operation - sorts all particles to ensure that all particles are located on
  * the proper cpu nodes (according to the domain decomposition). It is not high performance routine but it provides functionality
  * necessary for the global domain rearrangement, particles loading and repartitioning.
  *
  * \sa misc_markerPlacer.c
  * \sa misc_all2all.c
  */

#ifndef misc_markerPlacer
#define misc_markerPlacer					///< Guard against multiple include.

void placer_exchange (void);

#endif
