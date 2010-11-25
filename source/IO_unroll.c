/** \file IO_unroll.c
  * Analyzes memory allocation of the sub-region in memory and splits it into the biggest possible continuous chunks of memory.
  *
  * Routines presented in this file are used to optimize input/output of the subregion on mesh. For example, one may write:
  <pre>
  for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
    for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
      for (int k = reg->min[2] ; k <= reg->max[2] ; ++k)
        fwrite (&mv_f (mesh, i, j, k), sizeof (double), 1, fp);
  </pre>
  * and that will mean great many calls to \b fwrite. Than it may be optimized taking into account that nodes along k are always
  * consecutive and inner loop may be removed painlessly:
  <pre>
  for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
    for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
      fwrite (&mv_f (mesh, i, j, reg->min[2]), sizeof (double), reg->max[2] - reg->min[2] + 1, fp);
  </pre>
  * Routines presented here develops this idea further for case of 1D/2D simulation. Lets suppose we study 2D problem in XY plane.
  * It is 2D, so simulation domain is very big. Z axis is degenerated and that means that we still call \b fwrite for each node.
  *
  * Algorithm: data is stored in three-coordinate unrolled array. Module analyzes sizes of memory lines and strides and checks if
  * data spans are consecutive in memory. If they are, then set of spans is merged.
  *
  * \b Span is continuous piece of memory. Initially the entire subregion is a set of spans with first node corresponding to
  * address of mv_f(mesh, i, j, reg->min[2]) and size equal to reg->max[2] - reg->min[2] + 1. That program checks if spans cover
  * memory with no holes in between and joins spans if possible. Then all spans are enumerated and simple iterator allows to retrieve
  * them one by one to store/save/send/update and so on.
  *
  * Use it anywhere to get stored unrolled arrays by pieces which fit buffer.
  */

#ifndef MC_IO_UNROLL_MODULE
#define MC_IO_UNROLL_MODULE		///< Multiple include guard.

#include "log.h"
#include "type_mesh.h"

#include "misc_partition.h"

// ---------------------------------------------------------------------------
/// Structure to hold parameters of the span.
// ---------------------------------------------------------------------------
typedef struct
{
  int p;						///< Position of the cursor (\b p and \b q point span).
  int q;						///< Position of the cursor.
  int r;						///< Position of the cursor (\b r points position inside span).
  int dp;						///< Sizes of the spans (not necessary of the domain - spans may be merged).
  int dq;						///< Sizes of the spans (not necessary of the domain - spans may be merged).
  int dr;						///< Sizes of the spans (not necessary of the domain - spans may be merged).
  long int origin;					///< Start position of the first span (offset from the beginning of the storage array).
  long int strideP;					///< Stride along P direction.
  long int strideQ;					///< Stride along Q direction.
  long int strideR;					///< Stride along R direction.
  int done;						///< This flag is set to signal that unrolling is complete and all nodes were retrieved.
} spanSet_t;

// ---------------------------------------------------------------------------
/// \brief Decomposes memory occupied by sub-region \b reg in \b domain on continuous chunks which are enumerated in order of address' grow.
/// \b typeSize is the size of memory occupied by mesh node (may be \b double, \b vec3D_t, etc).
// ---------------------------------------------------------------------------
static void
span_init (const reg_t *domain, const reg_t *reg, spanSet_t *spans, int typeSize)
{
  ENSURE (reg_isInside (reg, domain), "region isn't inside of the domain");

  int sizeX = (domain->max[0] - domain->min[0])*mc_have_x + 1;								// Precalc of the sizes for strides.
  int sizeY = (domain->max[1] - domain->min[1])*mc_have_y + 1;
  int sizeZ = (domain->max[2] - domain->min[2])*mc_have_z + 1;

  long int sizes[6], strides[6];											// Stack of sizes.
  sizes[0] = (reg->max[2] - reg->min[2])*mc_have_z + 1;
  sizes[1] = (reg->max[1] - reg->min[1])*mc_have_y + 1;
  sizes[2] = (reg->max[0] - reg->min[0])*mc_have_x + 1;
  sizes[3] = 1;
  sizes[4] = sizes[5] = 0;

  strides[0] = mc_have_z*typeSize;											// Stack of strides.
  strides[1] = mc_have_y*typeSize*sizeZ;
  strides[2] = mc_have_x*typeSize*sizeY*sizeZ;
  strides[3] = typeSize*sizeX*sizeY*sizeZ;
  strides[4] = strides[5] = 0;

  spans->origin = (reg->min[0] - domain->min[0])*strides[2] + (reg->min[1] - domain->min[1])*strides[1] + 		// The start position of the first span.
                 (reg->min[2] - domain->min[2])*strides[0];

  // Checks if due to degeneracy (or coincidence of the region boundaries with mesh boundaries) spans may be merged.
  int gapDir = 0;
  for (int axis = 0 ; axis < 3 ; ++axis)
    if (strides[gapDir]*sizes[gapDir] == strides[gapDir+1] || !(strides[gapDir]))
      gapDir++;

  spans->dp = sizes[gapDir+2];												// 2D parametrization of all spans.
  spans->dq = sizes[gapDir+1];
  spans->dr = sizes[gapDir];
  spans->strideP = strides[gapDir+2];
  spans->strideQ = strides[gapDir+1];
  spans->strideR = strides[gapDir];

  spans->dr *= spans->strideR/typeSize;											// Granulalirity of r is always data-size.
  spans->strideR = typeSize;

  spans->p = spans->q = spans->r = spans->done = 0;									// Places cursor to the origin.
}

// ---------------------------------------------------------------------------
/// \brief Takes precalculated \b span and buffer \b capacity to return continuous chunk to process. Cursor is moved to new position (automatic update).
/// Returns size of the piece to process (guaranteed to be continuous of the biggest size to fit buffer) and offset of the data segment in \b shift.
// ---------------------------------------------------------------------------
static long int
span_iterate (spanSet_t *spans, int capacity, long *shift)
{
  *shift = spans->origin + spans->p*spans->strideP + spans->q*spans->strideQ + spans->r*spans->strideR;			// Gets old cursor position.

  long int retSize;
  if (capacity >= (spans->dr - spans->r)*spans->strideR)								// Checks if residual is small and fits buffer.
  {
    retSize = (spans->dr - spans->r)*spans->strideR;
    spans->r = 0;
    spans->q++;
    if (spans->q >= spans->dq)
    {
      spans->q = 0;
      spans->p++;
      if (spans->p >= spans->dp)
        spans->done = 1;
    }
  }
  else
  {
    retSize = capacity/spans->strideR;											// Returns number of data entries caller can
    spans->r += retSize;												// process using buffer with this capacity.
    retSize *= spans->strideR;
  }

  return retSize;
}

// ---------------------------------------------------------------------------
/// Checks if all spans are passed through unrolling sequence.
// ---------------------------------------------------------------------------
static int
span_allTouched (spanSet_t *span)
{
  return span->done;
}

#if 0
static void
dumpGeom (const char *msg, const reg_t *reg)
{
  say ("%s(%d, %d, %d) - (%d, %d, %d)", msg, reg->pos1[0], reg->pos1[1], reg->pos1[2], reg->pos2[0], reg->pos2[1], reg->pos2[2]);
}

static void
span_test (void)
{
  meshDouble_t mesh1, mesh2;
  spanSet_t spans1, spans2;
  reg_t reg = {{-1, -6, 11}, {3, 12, 12}, 0}, domain;

  reg.pos1[0] *= mc_have_x;											// Shrinks reg along degenerated axises.
  reg.pos1[1] *= mc_have_y;
  reg.pos1[2] *= mc_have_z;
  reg.pos2[0] *= mc_have_x;
  reg.pos2[1] *= mc_have_y;
  reg.pos2[2] *= mc_have_z;

  mesh_allocate (mcast_mesh (&mesh1), reg.pos1[0], reg.pos1[1], reg.pos1[2], reg.pos2[0], reg.pos2[1], reg.pos2[2], "testSpan:mesh1", mc_double);
  mesh_allocate (mcast_mesh (&mesh2), reg.pos1[0] - 1, reg.pos1[1] - 3, reg.pos1[2], reg.pos2[0], reg.pos2[1], reg.pos2[2] + 3, "testSpan:mesh2", mc_double);

  say ("Testing span coverage generator...");
  mf_mesh_dumpInfo (&mesh1);
  mf_mesh_dumpInfo (&mesh2);

  for (int i = reg.pos1[0] ; i <= reg.pos2[0] ; ++i)								// Randomizes content of the meshes.
    for (int j = reg.pos1[1] ; j <= reg.pos2[1] ; ++j)
      for (int k = reg.pos1[2] ; k <= reg.pos2[2] ; ++k)
      {
        mv_f(&mesh1, i, j, k) = rand ();
        mv_f(&mesh2, i, j, k) = rand ();
      }

  dumpGeom ("Testing copy of mesh using region ", &reg);
  domain.pos1[0] = mesh1.imin;
  domain.pos1[1] = mesh1.jmin;
  domain.pos1[2] = mesh1.kmin;
  domain.pos2[0] = mesh1.imax;
  domain.pos2[1] = mesh1.jmax;
  domain.pos2[2] = mesh1.kmax;
  span_init (&domain, &reg, &spans1, mf_mesh_sizeofNode (&mesh1));

  domain.pos1[0] = mesh2.imin;
  domain.pos1[1] = mesh2.jmin;
  domain.pos1[2] = mesh2.kmin;
  domain.pos2[0] = mesh2.imax;
  domain.pos2[1] = mesh2.jmax;
  domain.pos2[2] = mesh2.kmax;
  span_init (&domain, &reg, &spans2, mf_mesh_sizeofNode (&mesh2));

  while (! span_allTouched (&spans1))
  {
    long int offset;
    char buffer[1024];
    long int size = span_iterate (&spans1, 1024, &offset);							// Gets span and packs it.
    say ("pack: %ld(%ld bytes)", offset, size);
    memcpy (buffer, ((char*)mesh1.storage) + offset, size);

    int pos = 0;												// Puts unrolled data back to place.
    do {
      int s = span_iterate (&spans2, size - pos, &offset);
      memcpy (((char*)mesh2.storage) + offset, buffer + pos, s);
      say ("unpack: %ld(%ld bytes) from %d", offset, s, pos);
      pos += s;
    } while (size > pos);
  }

  for (int i = reg.pos1[0] ; i <= reg.pos2[0] ; ++i)								// Randomizes content of the meshes.
    for (int j = reg.pos1[1] ; j <= reg.pos2[1] ; ++j)
      for (int k = reg.pos1[2] ; k <= reg.pos2[2] ; ++k)
      {
/*        say ("(%02d, %02d, %02d): %ld -> %ld", i, j, k,
                  (long int)(mf_mesh_bytePointer(&mesh1, i, j, k) - (char*)mesh1.origin),
                  (long int)(mf_mesh_bytePointer(&mesh2, i, j, k) - (char*)mesh2.origin));*/
         say ("(%02d, %02d, %02d): %+e -> %+e (%+e)", i, j, k, mv_f(&mesh1, i, j, k), mv_f(&mesh2, i, j, k), mv_f(&mesh1, i, j, k) - mv_f(&mesh2, i, j, k));
      }

  mesh_free (mcast_mesh (&mesh1));
  mesh_free (mcast_mesh (&mesh2));
}
#endif

#endif
