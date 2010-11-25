/** \file type_marker.h
  * Plasma storage subsystem.
  */

#ifndef MC_TYPE_MARKER_HEADER
#define MC_TYPE_MARKER_HEADER

/**
  * Lagrangian marker (macro-particle) representing point on the characteristic line in 7D (r+v+t) space.
  */
typedef struct
{
  double x;			///< X coordinate.
  double y;			///< Y coordinate.
  double z;			///< Z coordinate.
  double vx;			///< \f$ \gamma v_x \f$ relativistic speed component.
  double vy;			///< \f$ \gamma v_y \f$ relativistic speed component.
  double vz;			///< \f$ \gamma v_z \f$ relativistic speed component.
  float  rho;			///< Weight of the marker (rho = q/(h1*h2*h3)).
  float  qDivM;			///< Charge to mass ratio of the marker.
} marker_t;

#endif
