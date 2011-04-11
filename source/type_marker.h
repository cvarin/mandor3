/** \file type_marker.h
  * Plasma storage subsystem.
  */

#ifndef MC_TYPE_MARKER_HEADER
#define MC_TYPE_MARKER_HEADER

/**
  * Lagrangian marker (macro-particle) representing the characteristic line in
  * 7D (r+v+t) phase space.
  */
typedef struct {
   double x;			///< X coordinate.
   double y;			///< Y coordinate.
   double z;			///< Z coordinate.
   double vx;			///< γ·v_x, γ is relativistic factor.
   double vy;			///< γ·v_y.
   double vz;			///< γ·v_z.
   float  rho;			///< Weight of the marker (rho = q/(h1*h2*h3)).
   float  qDivM;		///< Charge to mass ratio of the marker.
#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
   /// ID of the marker (FOR TRACER ONLY).
   int    id;
#endif
} marker_t;

#endif
