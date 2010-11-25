/** \file probe.h
  * \brief Structures and predefined constants for \b probe diagnostic. For futher details see example of the usage (diag_probe.h, probe2tec, ..).
  */

#ifndef mc_probe_header
#define mc_probe_header				///< \internal Guard

#include <string.h>				// For strncmp.

/**
  * Sample.
  */
typedef struct
{
  double time;			///< Time.
  double Ex, Ey, Ez;		///< E-field.
  double Hx, Hy, Hz;		///< H-field.
} probeSample_t;

/**
  * Header.
  */
typedef struct
{
  char   magic[8];		///< Magic is just string \b "probes" used to verify the content of the file(s).
  int    i;			///< I-index of the probe node.
  int    j;			///< J-index of the probe node.
  int    k;			///< K-index of the probe node.
} probeHeader_t;

/**
  * Header for packed parameters of the simulation (bounding box, number of nodes and so on).
  */
typedef struct
{
  double tau;			///< Time step.
  double h1;			///< Spatial step along X (zero for degenerated axises).
  double h2;			///< Spatial step along Y (zero for degenerated axises).
  double h3;			///< Spatial step along Z (zero for degenerated axises).
  double min[3];		///< Min-corner of the domain.
  double max[3];		///< Max-corner of the domain.
} probeGlobals_t;

#define mf_probeMagic(headerPointer)	(! strncmp ((headerPointer)->magic, "probes", 8))

#endif
