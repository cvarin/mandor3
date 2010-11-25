/** \file type_mesh.h
  * \brief Subroutines to create and work with meshes which allows negative indices and 1D/2D/3D access (governed by predefined constants
  * mc_have_x, mc_have_y, mc_have_z and macroses to avoid any memory/performance penalty). See type_mesh.h for macroses and type definitions.
  */

#ifndef mesh_header
#define mesh_header				///< \internal Multiple include protection.

#include <stdio.h>
#include <string.h>

#include "misc_partition.h"


/* ======================================================================================= */
/*                            Definitions of the node types.                               */
/* ======================================================================================= */

/**
  * Types for node of the interlaced vector field.
  */
typedef struct
{
  double x;				///< \b X component of vector.
  double y;				///< \b Y component of vector.
  double z;				///< \b Z component of vector.
} vec3D_t;

/**
  * Types for node of the interlaced vector field. \todo Consider unions?
  */
typedef struct
{
  double r[3];				///< Packed vector to access component by index.
} vec3Di_t;

/* Indexes to address components in r[3] array. */
#define mc_x 0				///< Symbolic name for \b X - component index.
#define mc_y 1				///< Symbolic name for \b Y - component index.
#define mc_z 2				///< Symbolic name for \b Z - component index.

/* Possible types of the nodes. */
#define mc_double   0			///< Type of the mesh field \b scalar (double).
#define mc_vec3D_t  1			///< Type of the mesh field \b vector (double<sup>3</sup>).

/**
  * That is fast replacement of the call to sizeof/switch/whatever. I use binary 'and' to truncate index to array limits.
  */
static const int mesh_sizeofArray[2] = {sizeof (double), sizeof (vec3D_t)};
#define mf_mesh_sizeofNode(mesh) (mesh_sizeofArray[((mesh)->type)&1])		///< \b Sizeof - type macros.

/**
  * Size of the mesh name.
  */
#define mc_meshNameSize 30



/* ======================================================================================= */
/*                          Template of the mesh definition.                               */
/* ======================================================================================= */

/**
  * \brief This template can take qualifier (const) and type of the pointers 'storage', 'origin'. So instead of making typecast using mesh->origin
  * I do tapecast one time for entire mesh. It removes many typecasts and hopefully helps to keep const qualifiers as long as possible.
  */
#define mf_makeMeshPrototype(mod, node)															\
{																			\
  mod int    imin, imax;	/* Basic sizes for mesh stored in storage and accessed through origin. 						*/	\
  mod int    jmin, jmax;																\
  mod int    kmin, kmax;																\
  mod int    type;		/* Type of the content (double, vec3D_t). 									*/	\
  mod int    width_z;		/* Shifts used in unrolling into 1D array. 									*/	\
  mod int    width_yz;																	\
  mod int    size;		/* Total size of data (in bytes!) = N_nodes*sizeof (double). 							*/	\
  mod node * mod storage;	/* Malloc'd storage array. All nodes after shift are mapped within this array. 					*/	\
  mod node * mod origin;	/* Shifted pointer to linearly shift from (i, j, k) -> (i - iMin, j - jMin, k - kMin) -> unroll into 1D array.	*/	\
  mod char   name[mc_meshNameSize];															\
}

/**
  * General mesh container, used for all non-math operations (copy, clear, check sizes and coverage, IO); I meah any operations I need to know only
  * size of memory chunk or sizeof of node in bytes. Typically I typecast to this type before calling functions from type_mesh.c module.
  */
typedef struct mesh_s    mf_makeMeshPrototype( , void) mesh_t;
typedef struct mesh_RO_s mf_makeMeshPrototype(const, void) mesh_RO_t;

/*
  * \brief ALL THIS TYPES ARE USED TO SPECIFY TYPE OF THE POINTERS 'storage', 'origin' WITHOUT TYPECAST (TO PRESERVE POSSIBLE CONST ATTRIBUTE
  * IN THE MV_* MACROSES). ALL FIELDS MUST BE ABSOLUTELY THE SAME TO ENSURE SAFETY OF THE  TYPECAST BETWEEN ALL MESH TYPES IN THE PROGRAM.
  * NOTE: RO stands for read-only content.
  */
typedef struct mf_makeMeshPrototype( , double) meshDouble_t;		///< Scalar field addressed by mv_f and mv_unrl macroses.
typedef struct mf_makeMeshPrototype( , vec3D_t) meshVec_t;		///< Vector field addressed by mv_fx, mv_fy, mv_fz, mv_v and mv_unrl macroses.
typedef struct mf_makeMeshPrototype( , vec3Di_t) meshVecI_t;		///< Vector field addressed by mv_fi, mv_v and mv_unrl macroses.
typedef struct mf_makeMeshPrototype(const, double) meshDouble_RO_t;	///< Scalar field addressed by mv_f and mv_unrl macroses.
typedef struct mf_makeMeshPrototype(const, vec3D_t) meshVec_RO_t;	///< Vector field addressed by mv_fx, mv_fy, mv_fz and mv_unrl macroses.
typedef struct mf_makeMeshPrototype(const, vec3Di_t) meshVecI_RO_t;	///< Vector field addressed by mv_fi and mv_unrl macroses.

/**
  * That is pointers to mesh structures.
  */
typedef mesh_t       *mesh_p;
typedef meshDouble_t *meshDouble_p;
typedef meshVec_t    *meshVec_p;
typedef meshVecI_t   *meshVecI_p;
typedef const mesh_RO_t       * const mesh_RO_p;
typedef const meshDouble_RO_t * const meshDouble_RO_p;
typedef const meshVec_RO_t    * const meshVec_RO_p;
typedef const meshVecI_RO_t   * const meshVecI_RO_p;

/**
  * That is typecast macroses to cast to particular type.
  */
#define mcast_mesh(mesh)          ((mesh_p)      (void*)(mesh))
#define mcast_meshDouble(mesh)    ((meshDouble_p)(void*)(mesh))
#define mcast_meshVec(mesh)       ((meshVec_p)   (void*)(mesh))
#define mcast_meshVecI(mesh)      ((meshVecI_p)  (void*)(mesh))
#define mcast_mesh_RO(mesh)       ((mesh_RO_p)      (void*)(mesh))
#define mcast_meshDouble_RO(mesh) ((meshDouble_RO_p)(void*)(mesh))
#define mcast_meshVec_RO(mesh)    ((meshVec_RO_p)   (void*)(mesh))
#define mcast_meshVecI_RO(mesh)   ((meshVecI_RO_p)  (void*)(mesh))

#undef mf_makeMeshPrototype



/* ======================================================================================= */
/*                          Template of the mesh definition.                               */
/* ======================================================================================= */


/**
  * Macroses to address particular element of mesh.
  */
#define mv_f(mesh, i, j, k)          ((mesh)->origin[(i)*mc_have_x*(mesh)->width_yz + (j)*mc_have_y*(mesh)->width_z + (k)*mc_have_z])
#define mv_fx(mesh, i, j, k)         ((mesh)->origin[(i)*mc_have_x*(mesh)->width_yz + (j)*mc_have_y*(mesh)->width_z + (k)*mc_have_z].x)
#define mv_fy(mesh, i, j, k)         ((mesh)->origin[(i)*mc_have_x*(mesh)->width_yz + (j)*mc_have_y*(mesh)->width_z + (k)*mc_have_z].y)
#define mv_fz(mesh, i, j, k)         ((mesh)->origin[(i)*mc_have_x*(mesh)->width_yz + (j)*mc_have_y*(mesh)->width_z + (k)*mc_have_z].z)
#define mv_fi(mesh, i, j, k, compnt) ((mesh)->origin[(i)*mc_have_x*(mesh)->width_yz + (j)*mc_have_y*(mesh)->width_z + (k)*mc_have_z].r[compnt])
#define mv_v(mesh, i, j, k)          ((mesh)->origin[(i)*mc_have_x*(mesh)->width_yz + (j)*mc_have_y*(mesh)->width_z + (k)*mc_have_z])
#define mv_unrl(mesh, pos)           ((mesh)->origin[pos])
#define mf_offset(mesh, i, j, k)     ((i)*mc_have_x*(mesh)->width_yz + (j)*mc_have_y*(mesh)->width_z + (k)*mc_have_z)

#define mf_mesh_bytePointer(mesh, i, j, k) (((char*)((mesh)->origin)) + mf_mesh_sizeofNode(mesh)*					\
                                            ((i)*mc_have_x*(mesh)->width_yz + (j)*mc_have_y*(mesh)->width_z + (k)*mc_have_z))


/* ======================================================================================= */
/*                              Some basic validity checks                                 */
/* ======================================================================================= */

/**
  * Checks against non-defined components.
  */
#define mf_checkComponentID(id, name)				\
{                                                      	 	\
  ENSURE ((id) == mc_x || (id) == mc_y || (id) == mc_z,  	\
          "%s: bad component ID (%d).", name, (id));   		\
}

/**
  * Checks against non-existing type of the content.
  */
#define mf_checkTypeID(id, name)					      \
{                                                       		      \
  ENSURE ((id) == mc_double || (id) == mc_vec3D_t, "%s: bad type ID (%d)",    \
          (name), (id)); 						      \
}

/**
  * Checks that point is inside of the domain.
  */
#define mf_mesh_checkBounds(mesh, i, j, k, name)					\
{											\
  if (((i) - (mesh)->imin)*((mesh)->imax - (i))*mc_have_x < 0 ||			\
      ((j) - (mesh)->jmin)*((mesh)->jmax - (j))*mc_have_y < 0 ||			\
      ((k) - (mesh)->kmin)*((mesh)->kmax - (k))*mc_have_z < 0)				\
    DIE ("%s: boundary is violated.", (name));			\
}

/**
  * Checks if meshes have different size.
  */
#define mf_mesh_cmpSize(mesh1, mesh2)											\
(  abs ((mesh1)->imin - (mesh2)->imin) + abs ((mesh1)->jmin - (mesh2)->jmin) + abs ((mesh1)->kmin - (mesh2)->kmin) +	\
   abs ((mesh1)->imax - (mesh2)->imax) + abs ((mesh1)->jmax - (mesh2)->jmax) + abs ((mesh1)->kmax - (mesh2)->kmax) > 0)

/**
  * Checks if point is inside mesh and can be read with no page-fault.
  */
#define mf_mesh_pointIsOutside(mesh, i, j, k)										\
( ((i) - (mesh)->imin)*mc_have_x < 0 || ((mesh)->imax - (i))*mc_have_x < 0 || 						\
  ((j) - (mesh)->jmin)*mc_have_y < 0 || ((mesh)->jmax - (j))*mc_have_y < 0 ||						\
  ((k) - (mesh)->kmin)*mc_have_z < 0 || ((mesh)->kmax - (k))*mc_have_z < 0 )


/**
  * Prints all fields of the mesh-structure.
  */
#define mf_mesh_dumpInfo(mesh)																\
{																			\
  char names[2][10] = {"double", "vec3D_t"};														\
  SAY_DEBUG ("mesh_dumpInfo:");																\
  SAY_DEBUG ("  - mesh name is %s, type '%s'.", (mesh)->name, names[(mesh)->type == mc_vec3D_t]);							\
  SAY_DEBUG ("  - region can be written as");														\
  SAY_DEBUG ("    o (%d, %d, %d)-(%d, %d, %d)", (mesh)->imin, (mesh)->jmin, (mesh)->kmin, (mesh)->imax, (mesh)->jmax, (mesh)->kmax);			\
  SAY_DEBUG ("    o [%d, %d] x [%d, %d] x [%d, %d].", (mesh)->imin, (mesh)->imax, (mesh)->jmin, (mesh)->jmax, (mesh)->kmin, (mesh)->kmax);		\
}



/* ======================================================================================= */
/*                          Some basic mem-chunk manipulations.                            */
/* ======================================================================================= */

/**
  * Cleans content of the storage array.
  */
#define mf_mesh_clean(mesh)			\
{						\
  memset ((mesh)->storage, 0, (mesh)->size);	\
}

/**
  * Copies content of one mesh to identically shaped another one.
  */
#define mf_mesh_copy(src, dest)							\
{										\
  ENSURE (!mf_mesh_cmpSize (dest, src) && (src)->type == (dest)->type,		\
          "source and destination meshes are different");			\
  memcpy ((dest)->storage, (src)->storage, (src)->size);			\
}



/* ======================================================================================= */
/*                                  Function prototypes.                                   */
/* ======================================================================================= */

void    mesh_allocate   (mesh_p mesh, int imin, int jmin, int kmin, int imax, int jmax, int kmax, const char *name, int type);
void    mesh_free       (mesh_p mesh);
void    mesh_resize     (mesh_p mesh, int imin, int jmin, int kmin, int imax, int jmax, int kmax);
void    mesh_reallocate (mesh_p mesh, int new_imin, int new_jmin, int new_kmin, int new_imax, int new_jmax, int new_kmax);
double *mesh_pack       (mesh_p mesh, int i1, int j1, int k1, int i2, int j2, int k2, int *size);
void    mesh_unpack     (mesh_p mesh, int i1, int j1, int k1, int i2, int j2, int k2, double *data, int size);

void    mesh_save   (mesh_RO_p mesh, const char *name, int imin, int jmin, int kmin, int imax, int jmax, int kmax);
void    mesh_load   (mesh_p mesh, FILE *fp, int shrinkToValidPart);
void    mesh_upload (mesh_p mesh, const reg_t *REG, FILE *fp);

#endif
