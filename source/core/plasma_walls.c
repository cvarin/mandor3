/** \file plasma_walls.c
  * Boundary conditions for particles.
  *
  * Plasma timestep is described in 'plasma.h'.
  *
  * Misc:
  *   o real boundary conditions (periodic, mirror, open) are applied first
  *
  *   o parallel exchange is implemented in plasma_parallel;
  *     main points:
  *       - incoming particles are received directly into shell array
  *       - outgoing particles are sent directly from plasma
  *       - communication is 'header+data' (two messages only)
  *
  *   o particles in core, or from the core, do not contribute to the boundary
  */

#ifndef MC_PLASMA_WALLS_HEADER
#define MC_PLASMA_WALLS_HEADER				///< \internal Guard.

/// Function to apply boundary condition to the shell.
typedef void (*pbc_func_t) (void);

// ---------------------------------------------------------------------------
/// Periodic boundary condition template.
// ---------------------------------------------------------------------------
#define MF_PERIODIC(POSTFIX, INDEX, TEST, BOUNDARY, SHIFT)		\
static void								\
pbc_periodic_ ## POSTFIX (void)						\
{									\
    for (marker_t *p = plasma, *end = p + countShell ; p < end ; ++p)	\
        if ((p-> INDEX) TEST (BOUNDARY))				\
            (p-> INDEX) += SHIFT;					\
}

MF_PERIODIC (xMin, x, <=, shellX1,  Lx)
MF_PERIODIC (xMax, x, >=, shellX2, -Lx)
MF_PERIODIC (yMin, y, <=, shellY1,  Ly)
MF_PERIODIC (yMax, y, >=, shellY2, -Ly)
MF_PERIODIC (zMin, z, <=, shellZ1,  Lz)
MF_PERIODIC (zMax, z, >=, shellZ2, -Lz)


// ---------------------------------------------------------------------------
/// Mirror boundary condition template.
// ---------------------------------------------------------------------------
#define MF_MIRROR(POSTFIX, INDEX, TEST, BOUNDARY)			\
static void								\
pbc_mirror_ ## POSTFIX (void)						\
{									\
    for (marker_t *p = plasma, *end = p + countShell ; p < end ; ++p)	\
        if ((p-> INDEX) TEST (BOUNDARY)) {				\
            (p-> INDEX) = 2*(BOUNDARY) - (p-> INDEX);			\
            (p->v ## INDEX) *= - 1;					\
        }								\
}

MF_MIRROR (xMin, x, <=, shellX1)
MF_MIRROR (xMax, x, >=, shellX2)
MF_MIRROR (yMin, y, <=, shellY1)
MF_MIRROR (yMax, y, >=, shellY2)
MF_MIRROR (zMin, z, <=, shellZ1)
MF_MIRROR (zMax, z, >=, shellZ2)


// ---------------------------------------------------------------------------
/// Sticky boundary condition template.
// ---------------------------------------------------------------------------
#define MF_DESTROY(POSTFIX, INDEX, TEST, BOUNDARY)			\
static void								\
pbc_destroy_ ## POSTFIX (void)						\
{									\
    marker_t *end = plasma + countShell;				\
    for (marker_t *p = plasma ; p < end ; ++p) {			\
        if ((p-> INDEX) TEST (BOUNDARY)) {				\
            *p = *(--end);						\
            --p;							\
        }								\
    }									\
    countShell = end - plasma;						\
}

MF_DESTROY (xMin, x, <=, shellX1 + h1)
MF_DESTROY (xMax, x, >=, shellX2 - h1)
MF_DESTROY (yMin, y, <=, shellY1 + h2)
MF_DESTROY (yMax, y, >=, shellY2 - h2)
MF_DESTROY (zMin, z, <=, shellZ1 + h3)
MF_DESTROY (zMax, z, >=, shellZ2 - h3)

// ---------------------------------------------------------------------------
/// Function does nothing - it replaces parallel exchange code which is called
/// independently and only after all local boundary conditions are applied.
// ---------------------------------------------------------------------------
static void dummy (void) {  return; }

// Map of boundary conditions: [face][condition]->func.
static pbc_func_t BC[6][BC_ENUM_LENGHT] = {
    {pbc_periodic_xMin, pbc_mirror_xMin, pbc_destroy_xMin, dummy},
    {pbc_periodic_yMin, pbc_mirror_yMin, pbc_destroy_yMin, dummy},
    {pbc_periodic_zMin, pbc_mirror_zMin, pbc_destroy_zMin, dummy},
    {pbc_periodic_xMax, pbc_mirror_xMax, pbc_destroy_xMax, dummy},
    {pbc_periodic_yMax, pbc_mirror_yMax, pbc_destroy_yMax, dummy},
    {pbc_periodic_zMax, pbc_mirror_zMax, pbc_destroy_zMax, dummy},
};

#endif
