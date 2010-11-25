/** \file tag_plasma_convex.c
  * Helpers to read and apply filters composed from many polyhedrons. File is
  * included directly into tag_plasma.c and is separated for readability only.
  */

#include <assert.h>

/// Plane in space 'A*x + B*y + C*z + D = 0'.
typedef struct {
   double  A, B, C, D;
} plane_t;

// ---------------------------------------------------------------------------
/// Region in space limited by infinite planes.
// ---------------------------------------------------------------------------
typedef struct {
  int      N;				///< Number of planes.
  plane_t *planes;			///< Array of planes.
} convex_t;

static int       polysN = 0;		///< Number of polyhedrons.
static convex_t *polys = NULL;		///< Polyhedrons.

// ---------------------------------------------------------------------------
/// Reads parameters of planes in uniformly indented block of config file.
// ---------------------------------------------------------------------------
static void
convex_read (FILE *fp)
{
   polys = (convex_t*) realloc (polys, (++polysN)*sizeof (convex_t));
   convex_t *b = polys + polysN - 1;
   b->N      = 0;
   b->planes = NULL;

   // Scans all '>' indented data blocks.
   while (cfg_isOption (fp)) {
      double rx, ry, rz, nx, ny, nz, shift;
      ENSURE (3 == fscanf (fp, "> %le %le %le %*[^\n] ", &rx, &ry, &rz),
              "cannot read reference point");
      ENSURE (3 == fscanf (fp, "> %le %le %le %*[^\n] ", &nx, &ny, &nz),
              "cannot read normal direction");
      ENSURE (1 == fscanf (fp, "> %le %*[^\n] ", &shift),
              "cannot read shift along the normal");

      // Computes factor to reset length of plane normal.
      double n = sqrt (nx*nx + ny*ny + nz*nz);
      ENSURE (n > 0.1,
              "too small normal length(%.3e), n = [%.3e %.3e %.3e]",
              n, nx, ny, nz);

      // Adds new plane to the storage.
      b->planes = (plane_t *) realloc (b->planes, (++b->N)*sizeof (plane_t));
      assert (b->planes);
      plane_t plane = { .A = nx/n,
                        .B = ny/n,
                        .C = nz/n,
                        .D = - (nx*rx + ny*ry + nz*rz)/n - shift };
      b->planes[b->N-1] = plane;
   }
}

// ---------------------------------------------------------------------------
/// Takes an array with markers and keeps only markers inside of the bounded region.
// ---------------------------------------------------------------------------
static int
convex_filter (marker_t *p, int N, convex_t *reg)
{
   for (marker_t *end = p + N ; p < end ; ++p) {
      // Tests against all planes.
      int outside = 0;
      for (plane_t *s = reg->planes, *s2 = s + reg->N ; (s < s2) && !outside ; ++s) {
         outside = (s->A*p->x + s->B*p->y + s->C*p->z + s->D > 0);
      }

      if (outside) {
         *p = *(--end);		// Replaces deleted particles by untested one.
         --p;			// Unrolls counter to test the new marker.
         --N;
      }
   }
   return N;
}

// ---------------------------------------------------------------------------
/// Removes loaded filters to prepare the tag for a new use.
// ---------------------------------------------------------------------------
static void
convex_free (void)
{
   for (convex_t *b = polys, *b2 = b + polysN ; b < b2 ; ++b) {
      free (b->planes);
   }
   free (polys);

   polysN = 0;
   polys  = NULL;
}
