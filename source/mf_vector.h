/** \file vector.h
  * \brief Macroses to do common vector manipulations regardless to the type.
  * \attention Vectors are supposed to be 3-dimentional and there are no type/limit checks.
  */

#ifndef mf_dotProduct

/// Dot product: \f$ (\vec a, \vec b) = a_x\cdot b_x + a_y\cdot b_y + a_z\cdot b_z \f$.
#define mf_dotProduct(a,b)		((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

/// Cross product: \f$ \vec c = [\vec a, \vec b] \f$.
#define mf_crossProduct(a,b,c)		((c)[0] = (a)[1]*(b)[2] - (a)[2]*(b)[1], (c)[1] = (a)[2]*(b)[0] - (a)[0]*(b)[2], (c)[2] = (a)[0]*(b)[1] - (a)[1]*(b)[0])

/// Vector subtraction: \f$ \vec c = \vec a - \vec b \f$.
#define mf_vectorSubtract(a, b, c)	((c)[0] = (a)[0] - (b)[0], (c)[1] = (a)[1] - (b)[1], (c)[2] = (a)[2] - (b)[2])

/// Vector addition: \f$ \vec c = \vec a + \vec b \f$.
#define mf_vectorAdd(a, b, c)		((c)[0] = (a)[0] + (b)[0], (c)[1] = (a)[1] + (b)[1], (c)[2] = (a)[2] + (b)[2])

/// Vector assignment: \f$ \vec b = \vec a \f$.
#define mf_vectorCopy(a, b)		((b)[0] = (a)[0], (b)[1] = (a)[1], (b)[2] = (a)[2])

/// Vector clearing: \f$ \vec a = 0 \f$.
#define mf_vectorClear(a)		((a)[0] = (a)[1] = (a)[2] = 0)

/// Vector negation: \f$ \vec b = - \vec a \f$.
#define mf_vectorNegate(a, b)		((b)[0] = - (a)[0], (b)[1] = - (a)[1], (b)[2] = - (a)[2])

/// Vector scaling: \f$ \vec b = \vec a * factor \f$.
#define mf_vectorScale(factor, a, b)	((b)[0] = (a)[0]*(factor), (b)[1] = (a)[1]*(factor), (b)[2] = (a)[2]*(factor))

/// Vector initalization: \f$ \vec r = (x, y, z) \f$.
#define mf_vectorSet(x, y, z, r)	((r)[0] = (x), (r)[1] = (y), (r)[2] = (z))

#endif
