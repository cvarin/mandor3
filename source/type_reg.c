/** \file type_reg.c
  * \brief Functions to perform simple manipulations with regions (see reg_t,
  * regList_t, type_reg.h).
  *
  * The main goal of this library is to provide simple language to create and
  * transform complex groups of subdomains. This task is very common for
  * parallel decompositions, adaptive mesh refinement patching, etc.
  *
  * \warning reg_printRanges(), reg_printWraps(): there are only \b 8 internal
  *          buffers to store the result so in case of bigger number of calls
  *          in one print statement the information will be corrupted
  *          (over-written by later calls).
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "type_reg.h"

#include "log.h"
#include "misc_parameters.h"		///< \todo Replace by \b dimensions.h.

// ---------------------------------------------------------------------------
/// Returns region's domain as string one may use with printf (internal buffer
/// stores up to 8 independent strings) - no free necessary.
// ---------------------------------------------------------------------------
const char *
reg_printRanges (const reg_t *reg)
{
   static char buffers[8][200];
   static int  num = -1;
   num = (num + 1) & 7;
   snprintf (buffers[num], 200, "[%d %d] [%d %d] [%d %d]",
             reg->min[0], reg->max[0], reg->min[1], reg->max[1], reg->min[2], reg->max[2]);
   buffers[num][199] = 0;			// Ensures NULL terminator.
   return buffers[num];
}

// ---------------------------------------------------------------------------
/// Returns region's domain as string one may use with printf (internal buffer
/// stores up to 8 independent strings) - no free necessary.
// ---------------------------------------------------------------------------
const char *
reg_printCorners (const reg_t *reg)
{
   static char buffers[8][200];
   static int  num = -1;
   num = (num + 1) & 7;
   snprintf (buffers[num], 200, "[%d %d %d] [%d %d %d]",
             reg->min[0], reg->min[1], reg->min[2], reg->max[0], reg->max[1], reg->max[2]);
   buffers[num][199] = 0;			// Ensures NULL terminator.
   return buffers[num];
}

// ---------------------------------------------------------------------------
/// Returns \b true if \b region is inside of the \b domain.
// ---------------------------------------------------------------------------
int
reg_isInside (const reg_t *region, const reg_t *domain)
{
  if ( (!mc_have_x || (domain->min[0] <= region->min[0] && domain->max[0] >= region->max[0])) &&
       (!mc_have_y || (domain->min[1] <= region->min[1] && domain->max[1] >= region->max[1])) &&
       (!mc_have_z || (domain->min[2] <= region->min[2] && domain->max[2] >= region->max[2])) )
    return 1;

  return 0;
}

// ---------------------------------------------------------------------------
/// Calculates common area of the \b region and \b domain; boundaries of
/// the domain are adjusted by \b ghost offsets (\b ghost = NULL means zero
/// offsets). Overlapped part is returned in the \b region itself, ret-value
/// \b 1 means no overlap.
// ---------------------------------------------------------------------------
int
reg_overlap (reg_t *region, const reg_t *domain, const reg_t *ghost)
{
   reg_t zero = {{0, 0, 0}, {0, 0, 0}};
   if (!ghost) {
      ghost = &zero;
   }

   // XXX use MIN, MAX.
   for (int ea = 0 ; ea < 3 ; ++ea) {
      region->min[ea] = (region->min[ea] < (domain->min[ea] + ghost->min[ea])) ? (domain->min[ea] + ghost->min[ea]) : region->min[ea];
      region->max[ea] = (region->max[ea] > (domain->max[ea] + ghost->max[ea])) ? (domain->max[ea] + ghost->max[ea]) : region->max[ea];
      region->min[ea] *= ACTIVATOR[ea];
      region->max[ea] *= ACTIVATOR[ea];
   }

   return reg_volume (region) <= 0;
}

// ---------------------------------------------------------------------------
/// Shifts region back to the normal position (inside of the actual
/// computational domain). Wrapping information stored in the reg_t::wrap field
/// is used to define the shift with no additional analysis.
// ---------------------------------------------------------------------------
void
reg_unwrap (reg_t *obj)
{
   for (int axis = 0 ; axis < 3 ; ++axis) {
      obj->min[axis] -= obj->wrap[axis];
      obj->max[axis] -= obj->wrap[axis];
      obj->wrap[axis] = 0;
   }
}

// ---------------------------------------------------------------------------
/// Takes partitioning of the main domain (baseMap), ghost cell offsets (to
/// resolve questions about who claims common boundary), request on wrapping
/// or extending (to find ghost cells' host cpu) and returns extended map.
///
/// Extended map is a map with additional copies (nodes may be wrapped around
/// domain in directions with non-zero wrap steps to handle periodic BC coupled
/// with parallel decomposition effects) and with domain extended to cover
/// local ghost cells if they have only one owner. This map can help to
/// identify unique owner of any mesh node for any cpu domain/ghost cell.
///
/// 'ghost == NULL' means no offsets.
// ---------------------------------------------------------------------------
void
reg_buildMap (const reg_t *dmn, regList_t *baseMap, regList_t *extMap, const reg_t *ghost, int wrap[3])
{
   reg_t zeroGhost = {{0, 0, 0}, {0, 0, 0}};
   if (!ghost) {
      ghost = &zeroGhost;
   }

   // Gets wrapping steps.
   for (int axis = 0 ; axis < 3 ; ++axis) {
      wrap[axis] = (wrap[axis]*ACTIVATOR[axis]) ? dmn->max[axis] - dmn->min[axis] : 0;
   }

   regList_clean (extMap);

   // Makes copy of the source partitioning ('27' to host all images possible).
   extMap->N = baseMap->N;
   extMap->list = (reg_t*) malloc (27*baseMap->N*sizeof (reg_t));
   memcpy (extMap->list, baseMap->list, baseMap->N*sizeof (reg_t));

   // Applies ghosts / extends domains.
   reg_t *reg = extMap->list;
   for (const reg_t * const end = reg + extMap->N ; reg < end ; ++reg) {
      reg->wrap[0] = reg->wrap[1] = reg->wrap[2] = 0;									// Cleans wrap.

      // Extends to grab ghost cells of the domain or shifts one one node
      // bounds to resolve claim to the common boundary.
      for (int axis = 0 ; axis < 3 ; ++axis) {
         // Uses ghost to attribute boundary to node.
         reg->min[axis] += (reg->min[axis] == dmn->min[axis] && ! wrap[axis]) ? -1 : ghost->min[axis];
         reg->max[axis] += (reg->max[axis] == dmn->max[axis] && ! wrap[axis]) ? +1 : ghost->max[axis];

         // Normalizes region for 1D/2D.
         reg->min[axis] *= ACTIVATOR[axis];
         reg->max[axis] *= ACTIVATOR[axis];
      }
   }

   // Produces wrapped copies.
   for (int axis = 0 ; axis < 3 ; ++axis) {
      if (!wrap[axis])
         continue;

      // Wraps all unprocessed regions.
      for (int sign = -1 ; sign < 2         ; sign += 2)
      for (int cpu  =  0 ; cpu  < extMap->N ; ++cpu) {
         // Drops region already wrapped along this axis.
         if (extMap->list[cpu].wrap[axis])
            continue;

         // Makes copy of region, wraps it, remembers inverse wrapping.
         extMap->list[extMap->N]            = extMap->list[cpu];
         extMap->list[extMap->N].min[axis] += sign*wrap[axis];
         extMap->list[extMap->N].max[axis] += sign*wrap[axis];
         extMap->list[extMap->N].wrap[axis] = sign*wrap[axis];
         ++(extMap->N);
      }
   }

   // Returns unused memory to system.
   extMap->list = (reg_t*) realloc (extMap->list, extMap->N*sizeof (reg_t));
}

// ---------------------------------------------------------------------------
/// Takes given list of regions and treats it as a partitioning map: given
/// region is tested against all map members and list of overlapped parts is
/// returned in \b list.
// ---------------------------------------------------------------------------
void
reg_distributeOnMap (const reg_t *reg, const regList_t *map, regList_t *list)
{
   int          resN = 0;			// Number of assigned pieces.
   reg_t       *res  = NULL;			// List of assigned pieces.
   const reg_t *node = map->list;

   for (const reg_t *end = node + map->N ; node < end ; ++node) {
      int   nonEmpty = 1;			// Volume of the node (in nodes).
      reg_t overlap;
      for (int axis = 0 ; axis < 3 ; ++axis) {	///< \todo Use reg_overlap().
         overlap.min[axis] = (node->min[axis] > reg->min[axis]) ? node->min[axis] : reg->min[axis];
         overlap.max[axis] = (node->max[axis] < reg->max[axis]) ? node->max[axis] : reg->max[axis];
         overlap.min[axis] *= ACTIVATOR[axis];
         overlap.max[axis] *= ACTIVATOR[axis];
         nonEmpty = nonEmpty && (overlap.max[axis] >= overlap.min[axis]);
      }

      // Adds new piece to the line.
      if (nonEmpty) {
         overlap.cpu     = node->cpu;
         overlap.wrap[0] = node->wrap[0];
         overlap.wrap[1] = node->wrap[1];
         overlap.wrap[2] = node->wrap[2];
         overlap.barcode = reg->barcode;
         res = (reg_t*) realloc (res, (++resN)*sizeof (reg_t));
         res[resN-1] = overlap;
      }
   }

   // Saves result in the envelope provided.
   list->N    = resN;
   list->list = res;
}

// ---------------------------------------------------------------------------
/// Calculates volume (number of nodes) of the region. Returns \b -1 in case of
/// badly sorted vertices of the region.
// ---------------------------------------------------------------------------
long long
reg_volume (const reg_t *reg)
{
   if ((mc_have_x && reg->max[0] < reg->min[0])
   ||  (mc_have_y && reg->max[1] < reg->min[1])
   ||  (mc_have_z && reg->max[2] < reg->min[2])) {
      return -1;
   }

  long long int res = 1;

#if mc_have_x
   res *= (reg->max[0] - reg->min[0] + 1);
#endif

#if mc_have_y
   res *= (reg->max[1] - reg->min[1] + 1);
#endif

#if mc_have_z
   res *= (reg->max[2] - reg->min[2] + 1);
#endif

   return res;
}

// ---------------------------------------------------------------------------
/// Cleans list of the regions.
// ---------------------------------------------------------------------------
void
regList_clean (regList_t *list)
{
   if (!list->N)
      return;

   free (list->list);
   list->N = 0;
   list->list = 0;
}

// ---------------------------------------------------------------------------
/// Adds new element to the list of the regions.
// ---------------------------------------------------------------------------
void
regList_add (regList_t *list, const reg_t *reg)
{
   list->list = (reg_t *) realloc (list->list, (++ list->N)*sizeof (reg_t));
   list->list[list->N-1] = *reg;
}

// ---------------------------------------------------------------------------
/// Slices all regions in the list by plane and stores all pieces in the same list.
// ---------------------------------------------------------------------------
void
regList_slice (regList_t *list, int axis, int coord)
{
   assert (axis >= 0 && axis < 3 && list && list->list);
   // Access by index only; list->list is changed in the cycle!
   for (int r = 0 ; r < list->N ; ++r) {															//   pointer is random (subject to realloc)!
      reg_t *reg = list->list + r;
      if (reg->min[axis] < coord && reg->max[axis] > coord) {
         reg_t r2 = *reg;													// Copy of sliced region.
         reg->max[axis] = r2.min[axis] = coord;										// Cuts top/bottoms out of two copies.
         regList_add (list, &r2);												// Note: list->list is changed!
      }
   }
}

// ---------------------------------------------------------------------------
/// \brief Checks that all elements have properly sorted vertices and all regions do not overlap. Returns \b 1 if there are any region with bad vertices,
/// and returns \b 2 if all regions are fine but come nodes are covered more than once.
// ---------------------------------------------------------------------------
int
regList_verify (regList_t *list)
{
  if (list && list->list)
  {
    reg_t *reg = list->list, *end = list->list + list->N;
    for ( ; reg < end ; ++reg)
      if (reg_volume (reg) <= 0)											// Checks that all regions are sorted normally.
        return 1;

    for ( ; reg < end ; ++reg)
      for (reg_t *test = reg + 1 ; test < end ; ++test)									// Checks that there are no overlapping.
        if ( (mc_have_x && reg->min[0] <= test->max[0] && reg->max[0] >= test->min[0]) ||
             (mc_have_y && reg->min[1] <= test->max[1] && reg->max[1] >= test->min[1]) ||
             (mc_have_z && reg->min[2] <= test->max[2] && reg->max[2] >= test->min[2]) )
          return 2;
  }
  return 0;
}

// ---------------------------------------------------------------------------
/// Simple constructor to generate region on the fly (uninited fields are
/// screwed on purpose to cause the error ASAP).
// ---------------------------------------------------------------------------
reg_t
reg_vv (const int *min, const int *max)
{
  reg_t r = {.min = {min[0], min[1], min[2]},
             .max = {max[0], max[1], max[2]},
             .cpu = -1, .barcode = 0,
             .wrap = {0, 0, 0}};
  return r;
}
