/** \file em_caps.c
  * Unit of the separate time step in the boundary layer aligned region (to
  * advance field for parallel exchange ASAP).
  *
  * List of interface functions and full description of the decompisition
  * details is documented in the em_caps.h.
  *
  * Implementation notes:
  * - This source file includes files emCap_Mur.c, emCap_mirror.c,
  *   emCap_periodic.c, emCap_decomposition.c directly to provide compiler
  *   strong opportunity to figure out function dependencies.
  * - TF/SF correction is added before updating of the overlapping region
  *   because of regions cap_t::toFlush are shrunk in the tangential direction.
  * - TF/SF correction is added before using boundary conditions to stop
  *   propagating of the perturbation due to interaction with boundary.
  *
  * \todo Move search of NaN into 'type_mesh.c' or into some kind of
  *       'misc_<whatever>.c'.
  */

#include <math.h>

#include "core/em_TFSF.h"
#include "scpic/em_sources.h"

#include "log.h"

#include "parr_meshes.h"
#include "parr_regLists.h"
#include "parr_ghostCells.h"

// ---------------------------------------------------------------------------
/// \brief Buffer for the old time-layer values of the electric field (stored
/// in the order <b>Ep0, Eq0, Ep1, Eq1</b> and cached value of the
/// \f$ \displaystyle \frac {c\tau - h}{c\tau + h} \f$ constant (see
/// emCap_Mur.c for theory). The rest of parameters is placed here to avoid
/// passing of the big amount of parameters and remove some repetitions (old
/// call was looking like <b>const int bound = (ea<<1) | 1;
/// capBC_Mur_E (mcast_meshVecI(&capsE[(ep<<1) | 0].mesh), ea, 1,
///              capsE[(ep<<1) | 0].packs + bound);</b>.
// ---------------------------------------------------------------------------
typedef struct
{
  double    *buffer;			///< Buffer for old time layer data.
  double     alpha;			///< \f$ \displaystyle \frac {c\tau - h}{c\tau + h} \f$, spatial step in direction perpenicular to boundary.
  int        ea;			///< Axis perpendicular to the boundary (axis limited by boundary).
  int        top;			///< Boundary to apply Mur's ABC to (\b top or \b bottom).
  meshVecI_p E;				///< Pointer to the cap's mesh.
} capMur_t;

// ---------------------------------------------------------------------------
/// Structures to hold all stuff relevant to boundary cache.
// ---------------------------------------------------------------------------
typedef struct
{
  meshVec_t   mesh;			///< Mesh to store new H/E field data on the new time layer.
  reg_t       toUpdate;			///< Region of mesh to update explicitly.
  reg_t       toFlush;			///< Region of mesh to write back to the parent mesh.
  regList_t   toCopy;			///< List of regions of the data to get from other caps.
  capMur_t    packs[6];			///< Six packs (for all boundaries).
} cap_t;

static cap_t  capsE[6];			///< Caps to store new E.
static cap_t  capsH[6];			///< Caps to store new H.

static int capAxisFrameFinder[6] = {0, 1, 2, 0, 1, 2};				///< Precalc for i % 3, 0 <= i < 6.

static const int capClaimE[4][2]  = { {0, 1} /* periodic */, {0, 1} /* mirror */, {0, 1}  /* Mur */, {-1, 1} /* parallel */};	///< Claims for electric field.
static const int capClaimH[4][2]  = { {0, 1} /* periodic */, {0, 1} /* mirror */, {0, 1}  /* Mur */, { 0, 1} /* parallel */};	///< Claims for electric field.

// ---------------------------------------------------------------------------
/// Update ranges for electric field (BC/top/range).
// ---------------------------------------------------------------------------
static const int capUpdateE[4][2][2] =
{
  { {0, 1}, {-1, -1} } /* periodic */,
  { {1, 1}, {-1,  0} } /* mirror */,
  { {1, 1}, {-1,  0} } /* Mur */,
  { {0, 1}, {-1, -1} } /* parallel */
};

// ---------------------------------------------------------------------------
/// Flush ranges for electric field (BC/top/range).
// ---------------------------------------------------------------------------
static const int capFlushE[4][2][2] =
{
  { {0, 1}, {-1, +1} } /* periodic */,
  { {0, 1}, {-1, +1} } /* mirror */,
  { {0, 1}, {-1, +1} } /* Mur */,
  { {0, 1}, {-1, -1} } /* parallel */
};

// ---------------------------------------------------------------------------
/// Update ranges for magnetic field (BC/top/range).
// ---------------------------------------------------------------------------
static const int capUpdateH[4][2][2] =
{
  { {0, 2}, {-1, 0} } /* periodic */,
  { {1, 2}, {-1, 0} } /* mirror */,
  { {0, 2}, {-1, 0} } /* Mur */,
  { {0, 2}, {-1, 1} } /* parallel */
};

// ---------------------------------------------------------------------------
/// Flush ranges for magnetic field (BC/top/range).
// ---------------------------------------------------------------------------
static const int capFlushH[4][2][2] =
{
  { {0, 2}, {-1, +1} } /* periodic */,
  { {0, 2}, {-1, +1} } /* mirror */,
  { {0, 2}, {-1, +1} } /* Mur */,
  { {0, 2}, {-1, +1} } /* parallel */
};

// ---------------------------------------------------------------------------
/// Scans given mesh in the given range and reports \b NaNs founded.
// ---------------------------------------------------------------------------
void
cap_reportNaNs (const char *msg, meshVec_p H, const reg_t *reg)
{
  SAY_DEBUG ("%s: scanning mesh %s in region %s...", msg, H->name, reg_printRanges (reg));

  mf_scanRegion(reg, i, j, k)
  {
    if (isnan (mv_fx (H, i, j, k)) || isnan (mv_fy (H, i, j, k)) || isnan (mv_fz (H, i, j, k)))
      SAY_DEBUG ("%s @ (%d, %d, %d) = (%e, %e, %e).", H->name, i, j, k, mv_fx (H, i, j, k), mv_fy (H, i, j, k), mv_fz (H, i, j, k));
  }
}

// Boundary conditions.
#include "emCap_Mur.c"
#include "emCap_mirror.c"
#include "emCap_periodic.c"
#include "emCap_decomposition.c"

// ---------------------------------------------------------------------------
/// \brief Sets main parameters of the caps and allocates meshes.
/// \b offsets hold offsets for top/bottom (directed inside of the domain, measured from official boundary of region).
// ---------------------------------------------------------------------------
static void
cap_1_initCaps (const int capUpdate[4][2][2], const int capFlush[4][2][2], cap_t *caps, const char field)
{
  const int bc[2][3] = { {cpu_bc_min[0], cpu_bc_min[1], cpu_bc_min[2]}, 	/// Packed bc. \todo Replace by ext pointers.
                         {cpu_bc_max[0], cpu_bc_max[1], cpu_bc_max[2]}};
  const reg_t domain = {{cpu_min[0], cpu_min[1], cpu_min[2]}, {cpu_max[0], cpu_max[1], cpu_max[2]}};// Node domain \todo Replace by reg_t cpu.
  const int * const sides[2] = {domain.min, domain.max};			// Packed sides of the domain.

  static const char axises[3] = {'X', 'Y', 'Z'}, alignment[2][4] = {"min", "max"};// Precalcs.

  SAY_DEBUG ("Adding caps with domain %s.", reg_printRanges (&domain));
  for (int b = 0 ; b < 6 ; ++b)							// Sets electric field caps.
  {
    const int ea = b >> 1, ep = (ea + 1) % 3, eq = (ea + 2) % 3, top = b & 1;	// Cap orientation.

    if (!ACTIVATOR[ea])								// Do not cap degenerated sides.
      continue;

    // Defines sizes of cap regions.
    reg_t *flush  = &caps[b].toFlush,
          *update = &caps[b].toUpdate;
    flush->min[ea]  = sides[top][ea] + capFlush [ bc[top][ea] ] [top] [0];
    flush->max[ea]  = sides[top][ea] + capFlush [ bc[top][ea] ] [top] [1];
    flush->min[ep]  = sides[0]  [ep] + capFlush [ bc[0]  [ep] ] [ 0 ] [0];
    flush->max[ep]  = sides[1]  [ep] + capFlush [ bc[1]  [ep] ] [ 1 ] [1];
    flush->min[eq]  = sides[0]  [eq] + capFlush [ bc[0]  [eq] ] [ 0 ] [0];
    flush->max[eq]  = sides[1]  [eq] + capFlush [ bc[1]  [eq] ] [ 1 ] [1];
    update->min[ea] = sides[top][ea] + capUpdate[ bc[top][ea] ] [top] [0];
    update->max[ea] = sides[top][ea] + capUpdate[ bc[top][ea] ] [top] [1];
    update->min[ep] = sides[0]  [ep] + capUpdate[ bc[0]  [ep] ] [ 0 ] [0];
    update->max[ep] = sides[1]  [ep] + capUpdate[ bc[1]  [ep] ] [ 1 ] [1];
    update->min[eq] = sides[0]  [eq] + capUpdate[ bc[0]  [eq] ] [ 0 ] [0];
    update->max[eq] = sides[1]  [eq] + capUpdate[ bc[1]  [eq] ] [ 1 ] [1];

    mf_reg_collapse (flush);
    mf_reg_collapse (update);

    // Sets cap mesh region.
    const char *name = _("cap%c[%c/%s]", field, axises[ea], alignment[top]);
    mesh_allocate (mcast_mesh (&caps[b].mesh),
                   flush->min[0], flush->min[1], flush->min[2],
                   flush->max[0], flush->max[1], flush->max[2],
                   name, mc_vec3D_t);

    SAY_DEBUG ("  - %s: '%s' / mesh '%s'", name, reg_printRanges (update), reg_printRanges (flush));// Debug message.
  }
}

// ---------------------------------------------------------------------------
/// Macros to cut overlapped bottom of the region.
// ---------------------------------------------------------------------------
#define mf_cutBtm(regName, axis)											\
      if (cap->regName.min[axis] == cached->regName.min[axis] && cap->regName.max[axis] > cached->regName.max[axis])	\
        cap->regName.min[axis] = cached->regName.max[axis] + 1;

// ---------------------------------------------------------------------------
/// Macros to cut overlapped top of the region.
// ---------------------------------------------------------------------------
#define mf_cutTop(regName, axis)											\
      if (cap->regName.max[axis] == cached->regName.max[axis] && cap->regName.min[axis] < cached->regName.min[axis])	\
        cap->regName.max[axis] = cached->regName.min[axis] - 1;

// ---------------------------------------------------------------------------
/// \brief "Removes" overlapping of the caps by reducing cap_t::toUpdate and cap_t::toFlush and generating cap_t::toCopy list.
// ---------------------------------------------------------------------------
static void
cap_2_shrink (cap_t *caps)
{
  cap_t *cap = caps + 2;
  for (int b = 2 ; b < 6 ; ++b, ++cap)
  {
    if (!ACTIVATOR[b >> 1])							// Skips uncapped boundaries.
      continue;

    regList_clean (&cap->toCopy);						// Clears list.
    const int checkSum = reg_volume (&cap->toUpdate);				// Initial region size.
    int check = 0;

    SAY_DEBUG ("Original for cap %d: %s", b, reg_printRanges (&cap->toUpdate));

    const int ea = b >> 1, ep = capAxisFrameFinder[ea+1], eq = capAxisFrameFinder[ea+2];// Frame of the cap.

    cap_t *cached = caps;
    for (int ready = 0 ; ready < b ; ++ready, ++cached)				// Checks what is prepared already.
    {
      if (!ACTIVATOR[ready >> 1] || (ready >> 1) == (b >> 1))			// Skips uncapped boundaries.
        continue;

      reg_t cmmn = cap->toUpdate;						// Gets overlapping region.
      reg_overlap (&cmmn, &cached->toUpdate, 0);				// Gets overlapping region.

      cmmn.barcode = ready;							// Barcode holds offset to source cap.
      regList_add (&caps[b].toCopy, &cmmn);					// Adds common part to copy-list.

      mf_cutTop(toUpdate, ep);							// Removes overlap.
      mf_cutBtm(toUpdate, ep);
      mf_cutTop(toUpdate, eq);
      mf_cutBtm(toUpdate, eq);

      mf_cutTop(toFlush, ep);
      mf_cutBtm(toFlush, ep);
      mf_cutTop(toFlush, eq);
      mf_cutBtm(toFlush, eq);

      #undef mf_cutBtm
      #undef mf_cutTop

/*      SAY_DEBUG ("  reduced to        %s by %s / cap %d", reg_printRanges (&cap->toUpdate), reg_printRanges (&cmmn), ready);*/
      int tmp = reg_volume (&cmmn);
      ENSURE (tmp > 0, "bad overlap");
      check += tmp;
    }

    ENSURE (check + reg_volume (&cap->toUpdate) == checkSum,
            "bad algorithm, %d nodes became %d nodes",
            checkSum, check + reg_volume (&cap->toUpdate));

    SAY_DEBUG ("  reduced to        %s.", reg_printRanges (&cap->toUpdate));
  }
}

// ---------------------------------------------------------------------------
/// \brief Applies all boundary conditions for local ghost cell content without caching and using caps.
/// ghostSync_sync() is used to finish the update for overlapped regions. This function may be called
/// many times so absorbing boundary condition should be implemented accordingly.
// ---------------------------------------------------------------------------
void
cap_resetGhosts (meshVec_p E, meshVec_p H)
{
  for (int axis = 0 ; axis < 3 ; ++axis)					// Applies periodic bc to all caps.
  {
    if (!ACTIVATOR[axis])
      continue;

    switch (cpu_bc_min[axis])
    {
      case BC_PERIODIC:
        capBC_periodic_H (H, H, axis);
        capBC_periodic_E (E, E, axis);
      break;

      case BC_MIRROR:
        capBC_mirror_H (mcast_meshVecI (H), axis, 0);
        capBC_mirror_E (mcast_meshVecI (E), axis, 0);
      break;

      default:
        DIE ("unknown BC");
      case BC_SPLITTER:
      case BC_OPEN:
      break;
    }

    switch (cpu_bc_max[axis])
    {
      case BC_MIRROR:
        capBC_mirror_H (mcast_meshVecI (H), axis, 1);
        capBC_mirror_E (mcast_meshVecI (E), axis, 1);
      break;

      default:
        DIE ("unknown BC");
      case BC_SPLITTER:
      case BC_PERIODIC:
      case BC_OPEN:
      break;
    }
  }
}

// ---------------------------------------------------------------------------
/// Creates all data structures to support data exchange between CPUs.
// ---------------------------------------------------------------------------
static void
cap_dump (cap_t *cap)
{
   SAY_DEBUG ("Dumping caps: name | to flush | to update:");
   for (int b = 0 ; b < 6 ; ++b) {
      if (!ACTIVATOR[b >> 1])
         continue;

      SAY_DEBUG ("  - %10s | %15s | %15s", cap[b].mesh.name,
                                           reg_printRanges (&cap[b].toFlush),
                                           reg_printRanges (&cap[b].toUpdate));
      for (int i = 0 ; i < cap[b].toCopy.N ; ++i)
         SAY_DEBUG ("    o gets data for %s from caps %d",
                                      reg_printRanges (cap[b].toCopy.list + i),
                                      cap[b].toCopy.list[i].barcode);
   }
}

// ---------------------------------------------------------------------------
/// Creates all data structures required to support data exchange between CPUs (framework).
// ---------------------------------------------------------------------------
void
cap_init (meshVec_p E, meshVec_p H)
{
   // Syncronization (done here to be independent on loading and initialization
   // of the ghost cells).
   cap_resetGhosts (E, H);	// Resets local    ghost cells.
   ghostSync_sync  (E, H);	// Resets parallel ghost cells.

   /// Move to 'cap_clean ()'.
   memset (capsE, 0, 6*sizeof (cap_t));
   memset (capsH, 0, 6*sizeof (cap_t));

   // Initializes cap sizes to resolve overlaps later.
   cap_1_initCaps (capUpdateE, capFlushE, capsE, 'E');
   cap_1_initCaps (capUpdateH, capFlushH, capsH, 'H');

   // Initializes exchange BC.
   capBC_parr_init ();

   // Handles overlaps.
   cap_2_shrink (capsE);
   cap_2_shrink (capsH);

   // Initializes Mur's absorbing boundary condition storages.
   const double h[3] = {h1, h2, h3};
   for (int ea = 0 ; ea < 3 ; ++ea) {
      if (!ACTIVATOR[ea])
         continue;

      const int ep = capAxisFrameFinder[ea+1],
                eq = capAxisFrameFinder[ea+2];

      if (cpu_bc_min[ea] == BC_OPEN) {
         const int bound = (ea << 1) | 0;
         capBC_Mur_init (mcast_meshVecI (E), &capsE[bound].mesh, ea, 0, capsE[bound].packs + bound, tau, h[ea]);
         if (ACTIVATOR[ep]) {
            capBC_Mur_init (mcast_meshVecI (E), &capsE[(ep<<1) | 0].mesh, ea, 0, capsE[(ep<<1) | 0].packs + bound, tau, h[ea]);
            capBC_Mur_init (mcast_meshVecI (E), &capsE[(ep<<1) | 1].mesh, ea, 0, capsE[(ep<<1) | 1].packs + bound, tau, h[ea]);
         }
         if (ACTIVATOR[eq]) {
            capBC_Mur_init (mcast_meshVecI (E), &capsE[(eq<<1) | 0].mesh, ea, 0, capsE[(eq<<1) | 0].packs + bound, tau, h[ea]);
            capBC_Mur_init (mcast_meshVecI (E), &capsE[(eq<<1) | 1].mesh, ea, 0, capsE[(eq<<1) | 1].packs + bound, tau, h[ea]);
         }
      }

      if (cpu_bc_max[ea] == BC_OPEN) {
         const int bound = (ea << 1) | 1;
         capBC_Mur_init (mcast_meshVecI (E), &capsE[bound].mesh, ea, 1, capsE[bound].packs + bound, tau, h[ea]);
         if (ACTIVATOR[ep]) {
            capBC_Mur_init (mcast_meshVecI (E), &capsE[(ep<<1) | 0].mesh, ea, 1, capsE[(ep<<1) | 0].packs + bound, tau, h[ea]);
            capBC_Mur_init (mcast_meshVecI (E), &capsE[(ep<<1) | 1].mesh, ea, 1, capsE[(ep<<1) | 1].packs + bound, tau, h[ea]);
         }
         if (ACTIVATOR[eq]) {
            capBC_Mur_init (mcast_meshVecI (E), &capsE[(eq<<1) | 0].mesh, ea, 1, capsE[(eq<<1) | 0].packs + bound, tau, h[ea]);
            capBC_Mur_init (mcast_meshVecI (E), &capsE[(eq<<1) | 1].mesh, ea, 1, capsE[(eq<<1) | 1].packs + bound, tau, h[ea]);
         }
      }
   }

   cap_dump (capsE);
   cap_dump (capsH);

   SAY_DEBUG ("cap_init: all done.");
}

// ---------------------------------------------------------------------------
/// Caches new magnetic field: calculates new field in the region \b toUpdate,
/// adds TF/SF correction in this region, expands domain of definition of the
/// capped field using \b toCopy list and data from cap-owners, and applies
/// boundary conditions to fill ghost cells with new data.
// ---------------------------------------------------------------------------
static void
cap_cacheH (meshVec_RO_p E, meshVec_p H)
{
   const double c1 = tau/h1*mc_have_x,
                c2 = tau/h2*mc_have_y,
                c3 = tau/h3*mc_have_z;

   // Caches new H field in toUpdate regions.
   for (int b = 0 ; b < 6 ; ++b) {
      if (!ACTIVATOR[b>>1])
         continue;
/*
      SAY_DEBUG ("Caching H: block %s to caps %d...", reg_printRanges (& capsH[b].toUpdate), b);*/
      meshVec_t *cH = &(capsH[b].mesh);
      mf_scanRegion(&capsH[b].toUpdate, i, j, k) {
         mv_fx(cH, i, j, k) =     mv_fx(H, i, j, k)
                            + c3*(mv_fy(E, i, j, k) - mv_fy(E, i, j, k-1))
                            - c2*(mv_fz(E, i, j, k) - mv_fz(E, i, j-1, k));
         mv_fy(cH, i, j, k) =     mv_fy(H, i, j, k)
                            + c1*(mv_fz(E, i, j, k) - mv_fz(E, i-1, j, k))
                            - c3*(mv_fx(E, i, j, k) - mv_fx(E, i, j, k-1));
         mv_fz(cH, i, j, k) =     mv_fz(H, i, j, k)
                            + c2*(mv_fx(E, i, j, k) - mv_fx(E, i, j-1, k))
                            - c1*(mv_fy(E, i, j, k) - mv_fy(E, i-1, j, k));
//          if (isnan (mv_fx (cH, i, j, k))
//          ||  isnan (mv_fy (cH, i, j, k))
//          ||  isnan (mv_fz (cH, i, j, k))) {
//             SAY_DEBUG ("Caching H = (%e, %e, %e) at (%d, %d, %d) / caps %d.", mv_fx (cH, i, j, k), mv_fy (cH, i, j, k), mv_fz (cH, i, j, k), i, j, k, b);
//             SAY_DEBUG ("H = (%e, %e, %e)", mv_fx (H, i, j, k), mv_fy (H, i, j, k), mv_fz (H, i, j, k));
//             SAY_DEBUG ("E = (%e, %e, %e)", mv_fx (E, i, j, k), mv_fy (E, i, j, k), mv_fz (E, i, j, k));
//             SAY_DEBUG ("E_x = (%e, %e, %e)", mv_fx (E, i-1, j, k), mv_fy (E, i-1, j, k), mv_fz (E, i-1, j, k));
//             SAY_DEBUG ("E_y = (%e, %e, %e)", mv_fx (E, i, j-1, k), mv_fy (E, i, j-1, k), mv_fz (E, i, j-1, k));
//             SAY_DEBUG ("E_z = (%e, %e, %e)", mv_fx (E, i, j, k-1), mv_fy (E, i, j, k-1), mv_fz (E, i, j, k-1));
//          }
      }
   }

   // Adds TF/SF contribution to the cached H.
   for (int b = 0 ; b < 6 ; ++b) {
      if (ACTIVATOR[b>>1]) {
         TFSF_postHStep (mcast_meshVecI(&capsH[b].mesh),       &capsH[b].toFlush);
         add_H_sources  (mcast_meshVecI(&capsH[b].mesh), Time, &capsH[b].toFlush);
      }
   }

   // Copies new H field to fill all internal nodes.
   for (int b = 0 ; b < 6 ; ++b) {
      if (!ACTIVATOR[b>>1])
         continue;

      // Copies previously cached data.
      meshVec_t   *cH   = &(capsH[b].mesh);
      const reg_t *copy =   capsH[b].toCopy.list;
      for (const reg_t * const end = copy + capsH[b].toCopy.N ; copy < end ; ++copy) {
         meshVec_t *src = & capsH[copy->barcode].mesh;
/*       SAY_DEBUG ("Updating region %s @ cap %d using %d cap.", reg_printRanges (copy), b, copy->barcode);*/
         /// \todo Try memcpy.
         mf_scanRegion(copy, i, j, k) {
            mv_v(cH, i, j, k) = mv_v (src, i, j, k);
/*          if (isnan(mv_fx (cH, i, j, k)) || isnan(mv_fy (cH, i, j, k)) || isnan(mv_fz (cH, i, j, k)))
               SAY_DEBUG ("Copiing H-NaN at (%d, %d, %d) / caps %d.", i, j, k, copy->barcode);*/
         }
      }
   }

   // Applies local BC to all caps.
   for (int ea = 0 ; ea < 3 ; ++ea) {
      if (!ACTIVATOR[ea])
         continue;

      const int ep = capAxisFrameFinder[ea+1],
                eq = capAxisFrameFinder[ea+2];

      switch (cpu_bc_min[ea]) {
         case BC_PERIODIC:
            capBC_periodic_H (   &capsH[(ea<<1) | 0].mesh, &capsH[(ea<<1) | 1].mesh, ea);
            if (ACTIVATOR[ep]) {
               capBC_periodic_H (&capsH[(ep<<1) | 0].mesh, &capsH[(ep<<1) | 0].mesh, ea);
               capBC_periodic_H (&capsH[(ep<<1) | 1].mesh, &capsH[(ep<<1) | 1].mesh, ea);
            }
            if (ACTIVATOR[eq]) {
               capBC_periodic_H (&capsH[(eq<<1) | 0].mesh, &capsH[(eq<<1) | 0].mesh, ea);
               capBC_periodic_H (&capsH[(eq<<1) | 1].mesh, &capsH[(eq<<1) | 1].mesh, ea);
            }
            break;

         case BC_MIRROR:
            capBC_mirror_H (   mcast_meshVecI(&capsH[(ea<<1) | 0].mesh), ea, 0);
            if (ACTIVATOR[ep]) {
               capBC_mirror_H (mcast_meshVecI(&capsH[(ep<<1) | 0].mesh), ea, 0);
               capBC_mirror_H (mcast_meshVecI(&capsH[(ep<<1) | 1].mesh), ea, 0);
            }
            if (ACTIVATOR[eq]) {
               capBC_mirror_H (mcast_meshVecI(&capsH[(eq<<1) | 0].mesh), ea, 0);
               capBC_mirror_H (mcast_meshVecI(&capsH[(eq<<1) | 1].mesh), ea, 0);
            }
            break;

         case BC_OPEN:
            capBC_Mur_H (mcast_meshVecI(&capsH[(ea<<1) | 0].mesh), ea, 0);
            if (ACTIVATOR[ep]) {
               capBC_Mur_H (mcast_meshVecI(&capsH[(ep<<1) | 0].mesh), ea, 0);
               capBC_Mur_H (mcast_meshVecI(&capsH[(ep<<1) | 1].mesh), ea, 0);
            }
            if (ACTIVATOR[eq]) {
               capBC_Mur_H (mcast_meshVecI(&capsH[(eq<<1) | 0].mesh), ea, 0);
               capBC_Mur_H (mcast_meshVecI(&capsH[(eq<<1) | 1].mesh), ea, 0);
            }
            break;

         default:
            DIE ("unknown BC");
         case BC_SPLITTER:
            break;
      }

      switch (cpu_bc_max[ea]) {
         case BC_MIRROR:
            capBC_mirror_H (mcast_meshVecI(&capsH[(ea<<1) | 1].mesh), ea, 1);
            if (ACTIVATOR[ep]) {
               capBC_mirror_H (mcast_meshVecI(&capsH[(ep<<1) | 0].mesh), ea, 1);
               capBC_mirror_H (mcast_meshVecI(&capsH[(ep<<1) | 1].mesh), ea, 1);
            }
            if (ACTIVATOR[eq]) {
               capBC_mirror_H (mcast_meshVecI(&capsH[(eq<<1) | 0].mesh), ea, 1);
               capBC_mirror_H (mcast_meshVecI(&capsH[(eq<<1) | 1].mesh), ea, 1);
            }
            break;

         case BC_OPEN:
            capBC_Mur_H (mcast_meshVecI(&capsH[(ea<<1) | 1].mesh), ea, 1);
            if (ACTIVATOR[ep]) {
               capBC_Mur_H (mcast_meshVecI(&capsH[(ep<<1) | 0].mesh), ea, 1);
               capBC_Mur_H (mcast_meshVecI(&capsH[(ep<<1) | 1].mesh), ea, 1);
            }
            if (ACTIVATOR[eq]) {
               capBC_Mur_H (mcast_meshVecI(&capsH[(eq<<1) | 0].mesh), ea, 1);
               capBC_Mur_H (mcast_meshVecI(&capsH[(eq<<1) | 1].mesh), ea, 1);
            }
            break;

         default:
            DIE ("unknown BC");
         case BC_SPLITTER:
         case BC_PERIODIC:
            break;
      }
   }
}

// ---------------------------------------------------------------------------
/// Calculates new magnetic field to use it as a RHS source, calculates new
/// field in the region \b toUpdate, adds TF/SF correction in this region,
/// expands domain of definition of the capped field using toCopy list and
/// data from cap-owners, and applies boundary conditions to fill ghost cells
/// with new data.
// ---------------------------------------------------------------------------
void
cap_cacheEH (meshVec_RO_p E, meshVec_p H)
{
   double c1 = tau/h1*mc_have_x,
          c2 = tau/h2*mc_have_y,
          c3 = tau/h3*mc_have_z;

   // Prepares new magnetic field at 't = time + tau*(n + 1/2)'.
   cap_cacheH (E, H);

   // Caches new E field in 'toUpdate' regions.
   for (int b = 0 ; b < 6 ; ++b) {
      if (!ACTIVATOR[b>>1])
         continue;
/*
    SAY_DEBUG ("Caching E: block %s to caps %d...", reg_printRanges (& capsE[b].toUpdate), b);*/
      meshVec_t *cH = &(capsH[b].mesh),
                *cE = &(capsE[b].mesh);
      mf_scanRegion(&capsE[b].toUpdate, i, j, k) {
         mv_fx(cE, i, j, k) =     mv_fx(E,  i, j,   k)
                            - c3*(mv_fy(cH, i, j,   k+1) - mv_fy(cH, i, j, k))
                            + c2*(mv_fz(cH, i, j+1, k  ) - mv_fz(cH, i, j, k));
         mv_fy(cE, i, j, k) =     mv_fy(E,  i,   j,   k)
                            - c1*(mv_fz(cH, i+1, j,   k) - mv_fz(cH, i, j, k))
                            + c3*(mv_fx(cH, i,   j, k+1) - mv_fx(cH, i, j, k));
         mv_fz(cE, i, j, k) =     mv_fz(E,  i,   j,   k)
                            - c2*(mv_fx(cH, i,   j+1, k) - mv_fx(cH, i, j, k))
                            + c1*(mv_fy(cH, i+1, j,   k) - mv_fy(cH, i, j, k));
//          if (isnan (mv_fx (cE, i, j, k) + mv_fy (cE, i, j, k)
//                                         + mv_fz (cE, i, j, k))) {
//             SAY_DEBUG ("Caching E = (%e, %e, %e) at (%d, %d, %d) / caps %d.", mv_fx (cE, i, j, k), mv_fy (cE, i, j, k), mv_fz (cE, i, j, k), i, j, k, b);
//             SAY_DEBUG ("E = (%e, %e, %e)", mv_fx (E, i, j, k), mv_fy (E, i, j, k), mv_fz (E, i, j, k));
//             SAY_DEBUG ("cH = (%e, %e, %e)", mv_fx (cH, i, j, k), mv_fy (cH, i, j, k), mv_fz (cH, i, j, k));
//             SAY_DEBUG ("cH_x = (%e, %e, %e)", mv_fx (cH, i+1, j, k), mv_fy (cH, i+1, j, k), mv_fz (cH, i+1, j, k));
//             SAY_DEBUG ("cH_y = (%e, %e, %e)", mv_fx (cH, i, j+1, k), mv_fy (cH, i, j+1, k), mv_fz (cH, i, j+1, k));
//             SAY_DEBUG ("cH_z = (%e, %e, %e)", mv_fx (cH, i, j, k+1), mv_fy (cH, i, j, k+1), mv_fz (cH, i, j, k+1));
//          }
      }
   }

   // Adds TF/SF contribution to the cached E.
   for (int b = 0 ; b < 6 ; ++b) {
      if (!ACTIVATOR[b>>1])
         continue;

      TFSF_postEStep (mcast_meshVecI(&capsE[b].mesh),       &capsE[b].toFlush);
//       add_E_sources  (mcast_meshVecI(&capsE[b].mesh), Time, &capsE[b].toFlush);
   }

   // Copies new E field to fill all internal nodes.
   for (int b = 0 ; b < 6 ; ++b) {
      if (!ACTIVATOR[b>>1])
         continue;

      // Copies previously cached data.
      meshVec_t   *cE   = &(capsE[b].mesh);
      const reg_t *copy = capsE[b].toCopy.list;
      for (const reg_t * const end = copy + capsE[b].toCopy.N ; copy < end ; ++copy) {
         meshVec_t *src = & capsE[copy->barcode].mesh;
/*      SAY_DEBUG ("Updating region %s @ cap %d using %d cap.", reg_printRanges (copy), b, copy->barcode);*/
         /// \todo Try memcpy.
         mf_scanRegion(copy, i, j, k) {
            mv_v(cE, i, j, k) = mv_v (src, i, j, k);
/*        if (isnan (mv_fx (cE, i, j, k)) || isnan (mv_fy (cE, i, j, k)) || isnan (mv_fz (cE, i, j, k)))
          SAY_DEBUG ("Copying E-NaN at (%d, %d, %d) / caps %d.", i, j, k, copy->barcode);*/
         }
      }
   }

   // Applies local BC to all caps.
   for (int ea = 0 ; ea < 3 ; ++ea) {
      if (!ACTIVATOR[ea])
         continue;

      const int ep = capAxisFrameFinder[ea+1],
                eq = capAxisFrameFinder[ea+2];

      switch (cpu_bc_min[ea]) {
      case BC_PERIODIC:
         capBC_periodic_E (&capsE[(ea<<1) | 0].mesh, &capsE[(ea<<1) | 1].mesh, ea);

         if (ACTIVATOR[ep]) {
            capBC_periodic_E (&capsE[(ep<<1) | 0].mesh, &capsE[(ep<<1) | 0].mesh, ea);
            capBC_periodic_E (&capsE[(ep<<1) | 1].mesh, &capsE[(ep<<1) | 1].mesh, ea);
         }

         if (ACTIVATOR[eq]) {
            capBC_periodic_E (&capsE[(eq<<1) | 0].mesh, &capsE[(eq<<1) | 0].mesh, ea);
            capBC_periodic_E (&capsE[(eq<<1) | 1].mesh, &capsE[(eq<<1) | 1].mesh, ea);
         }
         break;

      case BC_MIRROR:
         capBC_mirror_E (mcast_meshVecI(&capsE[(ea<<1) | 0].mesh), ea, 0);
         if (ACTIVATOR[ep]) {
            capBC_mirror_E (mcast_meshVecI(&capsE[(ep<<1) | 0].mesh), ea, 0);
            capBC_mirror_E (mcast_meshVecI(&capsE[(ep<<1) | 1].mesh), ea, 0);
         }
         if (ACTIVATOR[eq]) {
            capBC_mirror_E (mcast_meshVecI(&capsE[(eq<<1) | 0].mesh), ea, 0);
            capBC_mirror_E (mcast_meshVecI(&capsE[(eq<<1) | 1].mesh), ea, 0);
         }
         break;

      case BC_OPEN:
         {
            const int bound = (ea<<1) | 0;
            capBC_Mur_E (capsE[(ea<<1) | 0].packs + bound);
            if (ACTIVATOR[ep]) {
               capBC_Mur_E (capsE[(ep<<1) | 0].packs + bound);
               capBC_Mur_E (capsE[(ep<<1) | 1].packs + bound);
            }
            if (ACTIVATOR[eq]) {
               capBC_Mur_E (capsE[(eq<<1) | 0].packs + bound);
               capBC_Mur_E (capsE[(eq<<1) | 1].packs + bound);
            }
         }
         break;

      default:
         DIE ("unknown BC");
      case BC_SPLITTER:
         break;
      }

      switch (cpu_bc_max[ea]) {
      case BC_MIRROR:
         capBC_mirror_E (mcast_meshVecI(&capsE[(ea<<1) | 1].mesh), ea, 1);
         if (ACTIVATOR[ep]) {
            capBC_mirror_E (mcast_meshVecI(&capsE[(ep<<1) | 0].mesh), ea, 1);
            capBC_mirror_E (mcast_meshVecI(&capsE[(ep<<1) | 1].mesh), ea, 1);
         }
         if (ACTIVATOR[eq]) {
            capBC_mirror_E (mcast_meshVecI(&capsE[(eq<<1) | 0].mesh), ea, 1);
            capBC_mirror_E (mcast_meshVecI(&capsE[(eq<<1) | 1].mesh), ea, 1);
         }
         break;

      case BC_OPEN:
         {
            const int bound = (ea<<1) | 1;
            capBC_Mur_E (capsE[(ea<<1) | 1].packs + bound);
            if (ACTIVATOR[ep]) {
               capBC_Mur_E (capsE[(ep<<1) | 0].packs + bound);
               capBC_Mur_E (capsE[(ep<<1) | 1].packs + bound);
            }
            if (ACTIVATOR[eq]) {
               capBC_Mur_E (capsE[(eq<<1) | 0].packs + bound);
               capBC_Mur_E (capsE[(eq<<1) | 1].packs + bound);
            }
         }
         break;

      default:
         DIE ("unknown BC");
      case BC_SPLITTER:
      case BC_PERIODIC:
         break;
      }
   }

   // Sends data.
   capBC_parr_throw (E, H);
}

// ---------------------------------------------------------------------------
/// Flushes cached values of the new field in the boundary layers.
// ---------------------------------------------------------------------------
void
cap_flushHalfH (meshVec_p H)
{
/*  SAY_DEBUG ("Flushing 1/2 H:");*/
  for (int b = 0 ; b < 6 ; ++b)
  {
    if (!ACTIVATOR[b>>1])
      continue;

    meshVec_t *cH = &(capsH[b].mesh);
/*    SAY_DEBUG ("  - cap %d / reg %s / mesh [%d,%d]x[%d,%d]x[%d,%d]", b, reg_printRanges(&capsH[b].toFlush),
               cH->imin, cH->imax, cH->jmin, cH->jmax, cH->kmin, cH->kmax);*/
    mf_scanRegion(&capsH[b].toFlush, i, j, k)
    {
      mv_fx(H, i, j, k) = 0.5*(mv_fx(H, i, j, k) + mv_fx(cH, i, j, k));
      mv_fy(H, i, j, k) = 0.5*(mv_fy(H, i, j, k) + mv_fy(cH, i, j, k));
      mv_fz(H, i, j, k) = 0.5*(mv_fz(H, i, j, k) + mv_fz(cH, i, j, k));
    }
  }
}

// ---------------------------------------------------------------------------
/// Flushes cached values of the new field in the boundary layers.
// ---------------------------------------------------------------------------
void
cap_flushH (meshVec_p H)
{
/*  SAY_DEBUG ("Flushing 1/2 H:");*/
  for (int b = 0 ; b < 6 ; ++b)
  {
    if (!ACTIVATOR[b>>1])
      continue;

    meshVec_t *cH = &capsH[b].mesh;
/*    cap_reportNaNs ("Before flushing E", cH, &capsH[b].toFlush);
    SAY_DEBUG ("Flushing H field from cap %d to %s.", b, reg_printRanges (&capsH[b].toFlush));*/

    mf_scanRegion(&capsH[b].toFlush, i, j, k)
      mv_f(H, i, j, k) = mv_f(cH, i, j, k);
  }
}

// ---------------------------------------------------------------------------
/// Flushes cached values of the new field in the boundary layers.
// ---------------------------------------------------------------------------
void
cap_flushE (meshVec_p E)
{
  for (int b = 0 ; b < 6 ; ++b)
  {
    if (!ACTIVATOR[b>>1])
      continue;

    meshVec_t *cE = &(capsE[b].mesh);
/*    cap_reportNaNs ("Before flushing E", cE, &capsE[b].toFlush);
    SAY_DEBUG ("Flushing E field from cap %d to %s.", b, reg_printRanges (&capsE[b].toFlush));*/

    mf_scanRegion(&capsE[b].toFlush, i, j, k)
      mv_f(E, i, j, k) = mv_f(cE, i, j, k);
  }
/*
  reg_t reg = {{cpu_min[0], cpu_min[1], cpu_min[2]}, {cpu_max[0], cpu_max[1], cpu_max[2]}};
  cap_reportNaNs ("After cap_flushE", E, &reg);
  SAY_DEBUG ("Round is completed.\n\n");*/
}

// ---------------------------------------------------------------------------
/// Receives all parallel messages sended by cap_cacheEH().
// ---------------------------------------------------------------------------
void
cap_catchParallelMsgs (meshVec_p E, meshVec_p H)
{
   capBC_parr_catchH (H);
   capBC_parr_catchE (E);
   capBC_parr_finishRound ();
/*
   reg_t reg = {{cpu_min[0], cpu_min[1], cpu_min[2]}, {cpu_max[0], cpu_max[1], cpu_max[2]}};
   cap_reportNaNs ("After cap_flushE", E, &reg);
   SAY_DEBUG ("Round is completed.\n\n");	*/
}
