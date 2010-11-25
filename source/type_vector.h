/** \file type_vector.h
  * Type definition of integer 3-component vector and macroses to do common vector
  * manipulations regardless to the type.
  *
  * Main goal is to simplify writing of simple vector algebra operations and to simplify
  * access to the vector components by name (like 'vec.v.y') or by index (like 'vec.r[1]).
  * Access by name is very human readable and access by index is quite handy for loops.
  *
  * <h3>Special notes:</h3>
  * - vectors are supposed to be 3-dimentional (there are no type/limit checks in macroses);
  * - macroses operate on arrays with 3 elements (to access member of union transparently
  *   and to be used with any array of size 3 across the code);
  * - even in 1D/2D there are NO restriction on degenerated component of vector to be zero
  *   (this invariant is only for regions).
  *
  * <h2>Here I added few functions to do the job instead of macroses.</h2>
  *
  * Motivation is that I use this functions only in setup / startup, which is executed only
  * once and doesn't consume resources. That means that this part should have the best possible
  * readability but not performance, so macroses is probably an overkill. If functions will
  * work fine I'll keep them. I return vectors, not pointers, for simplicity of oneliners,
  * like
  *         c = vd_add (c, vd_scale (sin (alpha), vd (1, 0, 0)))
  */

#ifndef MF_VEC_DOT

/// Dot product: \f$ (\vec a, \vec b) = a_x\cdot b_x + a_y\cdot b_y + a_z\cdot b_z \f$.
#define MF_VEC_DOT(a, b)		((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

/// Cross product: \f$ \vec r = [\vec a, \vec b] \f$.
#define MF_VEC_CROSS(res, a, b)		((res)[0] = (a)[1]*(b)[2] - (a)[2]*(b)[1], (res)[1] = (a)[2]*(b)[0] - (a)[0]*(b)[2], (res)[2] = (a)[0]*(b)[1] - (a)[1]*(b)[0])

/// Vector subtraction: \f$ \vec r = \vec a - \vec b \f$.
#define MF_VEC_SUB(res, a, b)		((res)[0] = (a)[0] - (b)[0], (res)[1] = (a)[1] - (b)[1], (res)[2] = (a)[2] - (b)[2])

/// Vector addition: \f$ \vec r = \vec a + \vec b \f$.
#define MF_VEC_ADD(res, a, b)		((res)[0] = (a)[0] + (b)[0], (res)[1] = (a)[1] + (b)[1], (res)[2] = (a)[2] + (b)[2])

/// Vector assignment (copy): \f$ \vec r = \vec a \f$.
#define MF_VEC_COPY(res, a)		((res)[0] = (a)[0], (res)[1] = (a)[1], (res)[2] = (a)[2])

/// Vector clearing: \f$ \vec r = 0 \f$.
#define MF_VEC_CLEAR(res)		((res)[0] = (res)[1] = (res)[2] = 0)

/// Vector negation: \f$ \vec r = - \vec a \f$.
#define MF_VEC_NEG(res, a)		((res)[0] = - (a)[0], (res)[1] = - (a)[1], (res)[2] = - (a)[2])

/// Vector scaling: \f$ \vec r = \vec a * factor \f$.
#define MF_VEC_SCALE(res, factores, a)	((res)[0] = (a)[0]*(factores), (res)[1] = (a)[1]*(factores), (res)[2] = (a)[2]*(factores))

/// Vector packing: \f$ \vec r = (x, y, z) \f$.
#define MF_VEC_PACK(res, x, y, z)	((res)[0] = (x), (res)[1] = (y), (res)[2] = (z))

/// Vector unpacking: \f$ \vec r = (x, y, z) \f$.
#define MF_VEC_UNPACK(resX, resY, resZ, a)	((resX) = (a)[0], (resY) = (a)[1], (resZ) = (a)[2])

/// Checks if all components are the same: \f$ \vec a == \vec b \f$.
#define MF_VEC_EQ(a, b)			((a)[0] == (b)[0] && (a)[1] == (b)[1] && (a)[2] == (b)[2])

/**
  * Integer vector with 3 components stored as union (array/structure with named field).
  */
typedef union
{
   double r[3];		///< Vector as an array with 3 elements.

   struct {
      double x;		///< x component (mapped to r[0]).
      double y;		///< y component (mapped to r[1]).
      double z;		///< z component (mapped to r[2]).
   } v;			///< Vector as structure with 3 fields.
} vec3d_t;

/// Simple constructor from scalars.
static inline vec3d_t
vd (double x, double y, double z)
{
   vec3d_t v = { .r = {x, y, z} };
   return v;
}

/// Wraps 3 component into returned vec3d_t variable.
static inline vec3d_t
vd2 (const double *r)
{
   vec3d_t v = { .r = {r[0], r[1], r[2]} };
   return v;
}

/// Returns dot product '(a, b) = a_x*b_x + a_y*b_y + a_z*b_z'.
static inline double
vd_dot (vec3d_t a, vec3d_t b)
{
   return a.v.x*b.v.x + a.v.y*b.v.y + a.v.z*b.v.z;
}

/// Returns cross product '[a, b]', i.e. 'a×b'.
static inline vec3d_t
vd_cross (vec3d_t a, vec3d_t b)
{
   vec3d_t v;
   v.r[0] = a.r[1]*b.r[2] - a.r[2]*b.r[1];
   v.r[1] = a.r[2]*b.r[0] - a.r[0]*b.r[2];
   v.r[2] = a.r[0]*b.r[1] - a.r[1]*b.r[0];
   return v;
}

/// Returns 'a + b'.
static inline vec3d_t
vd_add (vec3d_t a, vec3d_t b)
{
   vec3d_t v;
   v.r[0] = a.r[0] + b.r[0];
   v.r[1] = a.r[1] + b.r[1];
   v.r[2] = a.r[2] + b.r[2];
   return v;
}

/// Returns 'a - b'.
static inline vec3d_t
vd_sub (vec3d_t a, vec3d_t b)
{
   vec3d_t v;
   v.r[0] = a.r[0] - b.r[0];
   v.r[1] = a.r[1] - b.r[1];
   v.r[2] = a.r[2] - b.r[2];
   return v;
}

/// Returns 'a*scale'.
static inline vec3d_t
vd_scale (vec3d_t a, double scale)
{
   vec3d_t v;
   v.r[0] = scale*a.r[0];
   v.r[1] = scale*a.r[1];
   v.r[2] = scale*a.r[2];
   return v;
}

const char* vd_print (vec3d_t v);

#endif
