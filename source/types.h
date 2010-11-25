/** \file types.h
  * This is types which do not require additional support in form of functions to work with objects of 
  * these types. 
  */

#ifndef mc_types_header
#define mc_types_header				///< Guard against multiple incudes.

typedef struct {float  r[3]} vec3f_t;		///< 3D-vector with single accuracy (struct to permit assigment '=' operation).
typedef struct {double r[3]} vec3d_t;		///< 3D-vector with double accuracy (struct to permit assigment '=' operation).
typedef struct {float x, y, z} xyz3f_t;		///< 3D-vector with single accuracy and named components.
typedef struct {double x, y, z} xyz3d_t;	///< 3D-vector with double accuracy and named components.

#endif
