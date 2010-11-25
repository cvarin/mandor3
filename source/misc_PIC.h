/** \file misc_PIC.h
  * Typical linear interpolation / distribution macroses for square mesh.
  *
  * \todo Make faster float->int conversion using SSE or other recipy from
  *       http://stackoverflow.com/questions/429632/how-to-speed-up-floating-point-to-integer-number-conversion
  */
#ifndef MC_MISC_PIC_HEADER
#define MC_MISC_PIC_HEADER

#include "type_mesh.h"

/**
  * Rounding downward (my version of gcc has bad implementation of 'fesetround' or
  * I do something wrong). Anyway I keep this workaround to have normal rounding for
  * negative indices. Another variant is in './source/real2int.c' but it is left
  * there just in case.
  */
#define MF_DBL_TO_INT(x) ((int) ((x) + 128.0) - 128)


/// Converts coordinate into index and diplacement (normalized on one).
#define MF_PIC_SIGMA(x, stepX, i, sigma)		\
{							\
  (sigma) = (x)/(stepX);				\
  (i) = MF_DBL_TO_INT (sigma);				\
  (sigma) -= (i);					\
}

/// Converts coordinate into index and diplacement for two meshes (saves one division).
#define MF_PIC_SIGMAS(x, stepX, i, sigma, i_, sigma_)	\
{							\
  (sigma) = (x)/(stepX);				\
  (sigma_) = (sigma) + 0.5;				\
  (i) = MF_DBL_TO_INT (sigma);				\
  (i_) = MF_DBL_TO_INT (sigma_);			\
  (sigma) -= i;						\
  (sigma_) -= i_;					\
}

/*
 *                                                                               oooo
 *  """M"""                M                                                        M                   M
 *     M      "Mo""""o    "M""""     o"""""o    ""Moo"""  "Mo""""o    o"""""o       M       o""""o     "M""""     o"""""o    ""Moo"""   o""""o"
 *     M       M     M     M        MoooooooM     M"       M      M  M       M      M       oooooM      M        M       M     M"       "o    "
 *     M       M     M     M        M             M        M      M  M       M      M      M"    M      M        M       M     M        o """"M
 *  oooMooo   oMoo ooMo    "oooo""   "Moooo""   ooMooo     M"oooo"    "ooooo"    oooMooo   "oooo"Moo    "oooo""   "ooooo"    ooMooo     M"oooM"
 *                                                         M
 *                                                        """""
 */

/* Lineyno (PIC) sobiraet znachenie s setki mesh using given memory mapper. */
#if 0
#define mf_PIC_general_read(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, func)\
(									\
  func(mesh, i, j, k)*(1 - (sigmaX))*(1 - (sigmaY))*(1 - (sigmaZ)) +	\
  func(mesh, i+1, j, k)*(sigmaX)*(1 - (sigmaY))*(1 - (sigmaZ)) +	\
  func(mesh, i, j+1, k)*(1 - (sigmaX))*(sigmaY)*(1 - (sigmaZ)) +	\
  func(mesh, i, j, k+1)*(1 - (sigmaX))*(1 - (sigmaY))*(sigmaZ) +	\
  func(mesh, i+1, j+1, k)*(sigmaX)*(sigmaY)*(1 - (sigmaZ)) +		\
  func(mesh, i, j+1, k+1)*(1 - (sigmaX))*(sigmaY)*(sigmaZ) +		\
  func(mesh, i+1, j, k+1)*(sigmaX)*(1 - (sigmaY))*(sigmaZ) +		\
  func(mesh, i+1, j+1, k+1)*(sigmaX)*(sigmaY)*(sigmaZ)			\
)
#endif

/*
 * That is multi-stage approach to PIC interpolation of a given value using memory
 * mapper function provided. Goal of it is to remove unnecessary computations
 * explicitly for 1D/2D (in this version we no longer need sigma{XYZ} for unreferenced
 * dimension(s)).
 */
#if mc_have_x
  #define mf_interpolate_X(mesh, func, i, sigmaX, j, k, expression)			\
    func(mesh, i, j, k)*(1 - (sigmaX))*(expression) + 					\
    func(mesh, (i)+1, j, k)*(sigmaX)*(expression)
#else
  #define mf_interpolate_X(mesh, func, i, sigmaX, j, k, expression)			\
    func(mesh, i, j, k)*(expression)
#endif

#if mc_have_y
  #define mf_interpolate_XY(mesh, func, i, sigmaX, j, sigmaY, k, expression)		\
    mf_interpolate_X(mesh, func, i, sigmaX, j, k, (1 - (sigmaY))*(expression)) +	\
    mf_interpolate_X(mesh, func, i, sigmaX, (j)+1, k, (sigmaY)*(expression))
#else
  #define mf_interpolate_XY(mesh, func, i, sigmaX, j, sigmaY, k, expression)		\
    mf_interpolate_X(mesh, func, i, sigmaX, j, k, (expression))
#endif

#if mc_have_z
  #define mf_PIC_general_read(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, func)		\
  (											\
    mf_interpolate_XY(mesh, func, i, sigmaX, j, sigmaY, k, (1 - (sigmaZ))) +		\
    mf_interpolate_XY(mesh, func, i, sigmaX, j, sigmaY, (k)+1, (sigmaZ))		\
  )
#else
  #define mf_PIC_general_read(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, func)		\
  (											\
    mf_interpolate_XY(mesh, func, i, sigmaX, j, sigmaY, k, 1)				\
  )
#endif


/* Lineyno (PIC) sobiraet scalar so scalar field mesh */
#define mf_PIC_read(mesh, i, sigmaX, j, sigmaY, k, sigmaZ)		\
(									\
  mf_PIC_general_read(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, mv_f)	\
)

/* Lineyno (PIC) sobiraet X-componentu s vector field mesh */
#define mf_PIC_readX(mesh, i, sigmaX, j, sigmaY, k, sigmaZ)		\
(									\
  mf_PIC_general_read(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, mv_fx)	\
)

/* Lineyno (PIC) sobiraet Y-componentu s vector field mesh */
#define mf_PIC_readY(mesh, i, sigmaX, j, sigmaY, k, sigmaZ)		\
(									\
  mf_PIC_general_read(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, mv_fy)	\
)

/* Lineyno (PIC) sobiraet Z-componentu s vector field mesh */
#define mf_PIC_readZ(mesh, i, sigmaX, j, sigmaY, k, sigmaZ)		\
(									\
  mf_PIC_general_read(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, mv_fz)	\
)


/*
 *                 M                                           M     ooo
 *  "M"""""o                           M                               M                      M
 *   M     "M    ""M       o""""o"    "M""""     ""Moo"""    ""M       Mo""""o   "M   ""M    "M""""     o"""""o    ""Moo"""   o""""o"
 *   M      M      M       "o    "     M           M"          M       M      M   M     M     M        M       M     M"       "o    "
 *   M     o"      M       o """"M     M           M           M       M      M   M     M     M        M       M     M        o """"M
 *  oMooooM"    oooMooo    M"oooM"     "oooo""   ooMooo     oooMooo   oM"oooo"    "Mooo"Mo    "oooo""   "ooooo"    ooMooo     M"oooM"
 *
 */

/*
 * Lineyno (PIC) raspredelyaet velichinu value na setku mesh m/y uzlami i, i+1, .., k, k+1 s vesami sigmaX, .., sigmaZ
 */
#if 0
#define mf_PIC_distribute(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, value)		  \
{										  \
  mv_f(mesh, i, j, k) += (value)*(1 - (sigmaX))*(1 - (sigmaY))*(1 - (sigmaZ));	  \
										  \
  mv_f(mesh, i+1, j, k) += (value)*(sigmaX)*(1 - (sigmaY))*(1 - (sigmaZ));	  \
  mv_f(mesh, i, j+1, k) += (value)*(1 - (sigmaX))*(sigmaY)*(1 - (sigmaZ));	  \
  mv_f(mesh, i, j, k+1) += (value)*(1 - (sigmaX))*(1 - (sigmaY))*(sigmaZ);	  \
										  \
  mv_f(mesh, i, j+1, k+1) += (value)*(1 - (sigmaX))*(sigmaY)*(sigmaZ);		  \
  mv_f(mesh, i+1, j, k+1) += (value)*(sigmaX)*(1 - (sigmaY))*(sigmaZ);		  \
  mv_f(mesh, i+1, j+1, k) += (value)*(sigmaX)*(sigmaY)*(1 - (sigmaZ));		  \
										  \
  mv_f(mesh, i+1, j+1, k+1) += (value)*(sigmaX)*(sigmaY)*(sigmaZ);		  \
}
#endif

/*
 * That is multi-stage approach to PIC distribution of a given value on given mesh.
 * Goal of it is to remove unnecessary computations explicitly for 1D/2D (in this
 * version we no longer need sigma? for unreferenced dimension(s)).
 */
#if mc_have_x
  #define mf_distribute_X(mesh, i, sigmaX, j, k, expression)			\
  {										\
    mv_f(mesh, i, j, k) += (1 - (sigmaX))*(expression);				\
    mv_f(mesh, (i)+1, j, k) += (sigmaX)*(expression);				\
  }
#else
  #define mf_distribute_X(mesh, i, sigmaX, j, k, expression)			\
  {										\
    mv_f(mesh, i, j, k) += (expression);					\
  }
#endif

#if mc_have_y
  #define mf_distribute_XY(mesh, i, sigmaX, j, sigmaY, k, expression)		\
  {										\
    mf_distribute_X(mesh, i, sigmaX, j, k, (1 - (sigmaY))*(expression));	\
    mf_distribute_X(mesh, i, sigmaX, (j)+1, k, (sigmaY)*(expression));		\
  }
#else
  #define mf_distribute_XY(mesh, i, sigmaX, j, sigmaY, k, expression)		\
  {										\
    mf_distribute_X(mesh, i, sigmaX, j, k, (expression));			\
  }
#endif

#if mc_have_z
  #define mf_PIC_distribute(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, value)	\
  {										\
    mf_distribute_XY(mesh, i, sigmaX, j, sigmaY, k, (1 - (sigmaZ))*(value));	\
    mf_distribute_XY(mesh, i, sigmaX, j, sigmaY, (k)+1, (sigmaZ)*(value));	\
  }
#else
  #define mf_PIC_distribute(mesh, i, sigmaX, j, sigmaY, k, sigmaZ, value)	\
  {										\
    mf_distribute_XY(mesh, i, sigmaX, j, sigmaY, k, (value));			\
  }
#endif

#endif
