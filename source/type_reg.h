/** \file type_reg.h
  * Fundamental type to work with mesh subdomains, bounding boxes, etc.
  *
  * <h2>Motivation.</h2>
  * All algorithms of decompositions, implementations of boundary condition, 1D/2D/3D general
  * array unrolling and so on deals with fundamental abstraction - 'sub-regions of mesh'. Type
  * reg_t is designed to simplify all work with such objects: reg_t may uniformly represent
  * line, rectangular or volume of mesh nodes where all faces and edges are aligned to the
  * axises and axis planes. Uniformly means that all geometric transformations (joining,
  * splitting, counting of the nodes number, etc) and processing are done without special care
  * about dimensionality.
  *
  * Set of routines developed here is used to do following operations with regions: allocate,
  * print/dump, compare, split, partition on given set of nodes and maintain a groups of regions
  * in form of lists. Parallel aspect of geometrical decompositions is handled in the parr_regLists.h.
  *
  * Special attention is paid to handle periodic wrapping of the regions around directions with
  * periodic boundary conditions and accounting of the 1D/2D degenerated cases.
  *
  * <h3>Notations used to describe the regions in the rest of the documentation.</h3>
  * Region in the mesh space may be written as a union of the intervals where indices are reside,
  * like '[imin, imax] x [jmin, jmax] x [kmin, kmax]' which means {i, j, k | imin <= i <= imax,
  * jmin <= j <= jmax, kmin <= k <= kmax"}. Sometimes instead of interval '[min, max]' a direct set
  * of all possible values of index is given: <b>{i1, i2, i3, ..}</b>. For the sake of clarity
  * minimal values of indexes sometimes are written as '0' and maximal value of index as 'I',
  * 'J', 'K' (depending on axis). Usig this way a domain of current cpu is written the most
  * shortly way as '[0, I] x [0, J] x [0, K]'.
  *
  * \attention Any function is permitted to reset min/max/wrap fields along degenerated
  * \axis(es) into 0 at any time at any place.
  *
  * \attention Coordinates of the corners in reg_t stucture must be sorted: '\b min' vector
  * coordinates should never be bigger that coordinates of '\b max'!
  */

#ifndef MC_TYPE_REG_HEADER
#define MC_TYPE_REG_HEADER		///< Guard against multiple incudes.

/**
  * Region is set by 'min' corner and 'max' corner. Additional information field 'barcode'
  * is used by creator. Two corners set the main diagonal of the object and only two points is
  * necessary to set either line, rectangular or volume if all faces/edges are aligned to the
  * coordinate axises and planes. This structure is also used to specify unrolling information
  * for the 1D/2D/3D meshes, in which case 'cpu' field stores 'shift' and 'wrap' field stores
  * 'strides' (see 'mv_map' and 'mf_mesh_updateMapper' macroses).
  */
typedef struct
{
  int min[3];		///< Corner of the region with minimal coordinates.
  int max[3];		///< Corner of the region with maximal coordinates.
  int cpu;		///< Cpu number of the object owner / shift of the origin used to map reg_t::min to the beginning of the unrolled array.
  int barcode;		///< Auxilary information tagged by creator to travel with the region.
  int wrap[3];		///< Wrapping of the region (used to trace the shift of region due to periodicity) / strides for the mesh.
} reg_t;

/**
  * \brief Group of regions. Used to hold regions for processing, passing and so on. Field \b list is a pointer to reallocatable memory chunk.
  */
typedef struct
{
  int    N;				///< Number of regions in the array.
  reg_t *list;				///< Array of regions.
} regList_t;

#define mc_reg_init     {{0, 0, 0}, {0, 0, 0}, 0, 0, {0, 0, 0}}			///< reg_t init constants (all zeroes).
#define mc_reg_initBad  {{-1, -1, -1}, {-5, -5, -5}, -1, -1, {-1, -1, -1}}	///< reg_t init (badly sorted to corrupt everything if used).
#define mc_regList_init {0, 0}							///< regList_t defaults (\b realloc is used => NULL pointer \b must be start value).

/**
  * Resets sizes along degenerated axises.
  */
#define mf_reg_collapse(reg)							\
{										\
  (reg)->min[0] *= mc_have_x;							\
  (reg)->min[1] *= mc_have_y;							\
  (reg)->min[2] *= mc_have_z;							\
  (reg)->max[0] *= mc_have_x;							\
  (reg)->max[1] *= mc_have_y;							\
  (reg)->max[2] *= mc_have_z;							\
  (reg)->wrap[0] *= mc_have_x;							\
  (reg)->wrap[1] *= mc_have_y;							\
  (reg)->wrap[2] *= mc_have_z;							\
}

/**
  * Defines 3 \b for loops to scan entire region (to save typing).
  */
#define mf_scanRegion(reg, i, j, k)										\
  for (int i = (reg)->min[0]*mc_have_x ; i <= (reg)->max[0]*mc_have_x ; ++i)					\
    for (int j = (reg)->min[1]*mc_have_y ; j <= (reg)->max[1]*mc_have_y ; ++j)					\
      for (int k = (reg)->min[2]*mc_have_z ; k <= (reg)->max[2]*mc_have_z ; ++k)

const char* reg_printRanges (const reg_t *reg);
const char* reg_printCorners (const reg_t *reg);

int  reg_isInside (const reg_t *region, const reg_t *domain);
int  reg_overlap (reg_t *region, const reg_t *domain, const reg_t *ghost);
long long reg_volume (const reg_t *reg);

reg_t reg_vv (const int *min, const int *max);

void reg_unwrap (reg_t *obj);
void reg_buildMap (const reg_t *dmn, regList_t *baseMap, regList_t *extMap, const reg_t *ghost, int wrap[3]);
void reg_distributeOnMap (const reg_t *reg, const regList_t *map, regList_t *list);

void regList_clean (regList_t *list);
void regList_add (regList_t *list, const reg_t *reg);
void regList_slice (regList_t *list, int axis, int coord);
int  regList_verify (regList_t *list);

#endif
