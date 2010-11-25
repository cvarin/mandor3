/** \file type_mesh2.h
  * \brief Set of routines to work with unrolled arrays by using reg_t map to access
  * elements of the linear storage array.
  *
  * Main ideas:
  * - mesh is always an unrolled array plus something else on top to simplify access
  * - this something else is reg_t structure which already holds all sizes in usable
  *   form and have enough free fields to hold temporaries (3 strides and origin bias)
  * - separation of the storage may help to use \b restric keyword during the work and
  *   also to access a storage array as whole thing
  */

#ifndef MC_TYPE_MESH2_HEADER
#define MC_TYPE_MESH2_HEADER				///< \internal Multiple include protection.

#include <stdio.h>
#include <string.h>

#include "misc_parameters.h"

#include "log.h"
#include "type_reg.h"
#include "type_mesh.h"

/**
  * Updates fields \b wrap (strides) and \b cpu (shift of origin) in reg_t
  * structure used later by 'mv_map' macros as unrolling map.
  */
#define mf_mesh_initUnroll(reg)											\
{														\
  mf_reg_collapse(reg);												\
  (reg)->wrap[2] = mc_have_z;											\
  (reg)->wrap[1] = mc_have_y*((reg)->max[2] - (reg)->min[2] + 1);						\
  (reg)->wrap[0] = mc_have_x*((reg)->max[1] - (reg)->min[1] + 1)*((reg)->max[2] - (reg)->min[2] + 1);		\
  (reg)->cpu = - ((reg)->wrap[0]*(reg)->min[0] + (reg)->wrap[1]*(reg)->min[1] + (reg)->wrap[2]*(reg)->min[2]);	\
  ENSURE (reg_volume (reg) > 0, "bad unrolling mapper");			\
}

/**
  * \brief Maps (i, j, k)-node on the unrolled 1D coordinate. Array limits and unroll info are stored in
  * the reg_t-type \b map (see mf_mesh_initUnroll).
  */
#define mf_mesh2_offset(map, i, j, k)	((map)->cpu + (i)*mc_have_x*(map)->wrap[0] + (j)*mc_have_y*(map)->wrap[1] + (k)*mc_have_z)

/**
  * \brief Maps (i, j, k)-node on the location inside of the unrolled array. Array limits and unroll info
  * are stored in the reg_t-type \b map (see mf_mesh_initUnroll).
  */
#define mv_map(map, array, i, j, k)	((array)[mf_mesh2_offset(map, i, j, k)])

void  mesh2_save (const char *name, const reg_t *saveReg, const reg_t *map, const char *array, int elemSize);
void* mesh2_load (const char *name, reg_t *map, int *elemSize);
void  mesh2_mesh (reg_t *map, void *array, int elemSize, mesh_p mesh);

#endif
