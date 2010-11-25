/** \file tag_plasma_box.c
  * Helpers to read and apply filters composed from many boxes. File is
  * included directly into tag_plasma.c and is separated for readability only.
  */

#include <assert.h>

/// Plane in space 'A*x + B*y + C*z + D = 0'.
typedef struct
{
  vec3d_t  center, 			// Position of the center.
           dir1, 			// Half-width along direction 1.
           dir2, 			// Half-width along direction 2.
           dir3; 			// Half-width along direction 3.
} box_t;

static int    boxesN = 0;		///< Number of polyhedrons.
static box_t *boxes = NULL;		///< Polyhedrons.

// ---------------------------------------------------------------------------
/// Reads parameters of planes in uniformly indented block of config file.
// ---------------------------------------------------------------------------
static void
box_read (FILE *fp)
{
  while (cfg_isOption (fp))							// Scans all '>' indented data blocks.
  {
    boxes = (box_t*) realloc (boxes, (++boxesN)*sizeof (box_t));			// Adds new block of planes.
    box_t *b = boxes + boxesN - 1;
    assert (boxes);

    ENSURE (3 == fscanf (fp, "> %le %le %le %*[^\n] ", &b->center.v.x, &b->center.v.y, &b->center.v.z),
               "Cannot read position of center.");
    ENSURE (3 == fscanf (fp, "> %le %le %le %*[^\n] ", &b->dir1.v.x, &b->dir1.v.y, &b->dir1.v.z), "Cannot read direction 1.");
    ENSURE (3 == fscanf (fp, "> %le %le %le %*[^\n] ", &b->dir2.v.x, &b->dir2.v.y, &b->dir2.v.z), "Cannot read direction 2.");
    ENSURE (3 == fscanf (fp, "> %le %le %le %*[^\n] ", &b->dir3.v.x, &b->dir3.v.y, &b->dir3.v.z), "Cannot read direction 3.");

    double cmp = 0.001*(mc_have_x*h1*h1 + mc_have_y*h2*h2 + mc_have_z*h3*h3);
    ENSURE (vd_dot (b->dir1, b->dir1) > cmp, "Too small direction-1 vector(%s).", vd_print (b->dir1));
    ENSURE (vd_dot (b->dir2, b->dir2) > cmp, "Too small direction-1 vector(%s).", vd_print (b->dir2));
    ENSURE (vd_dot (b->dir3, b->dir3) > cmp, "Too small direction-1 vector(%s).", vd_print (b->dir3));

    b->dir1 = vd_scale (b->dir1, 1.0/vd_dot (b->dir1, b->dir1));
    b->dir2 = vd_scale (b->dir2, 1.0/vd_dot (b->dir2, b->dir2));
    b->dir3 = vd_scale (b->dir3, 1.0/vd_dot (b->dir3, b->dir3));
  }
}

// ---------------------------------------------------------------------------
/// Takes an array with markers and keeps only markers inside of the bounded region.
// ---------------------------------------------------------------------------
static int
box_filter (marker_t *p, int N, box_t *b)
{
  for (marker_t *end = p + N ; p < end ; ++p)
  {
    vec3d_t diff = vd_sub (vd (p->x, p->y, p->z), b->center);			// Distance between particle and center.

    if (fabs (vd_dot (b->dir1, diff)) > 1 || 					// Removes particles outside of the box.
        fabs (vd_dot (b->dir2, diff)) > 1 ||
        fabs (vd_dot (b->dir3, diff)) > 1)
    {
      *p = *(--end);								// Replaces deleted particles by untested one.
      --p;									// Unrolls counter to test the new marker.
      --N;
    }
  }
  return N;
}

// ---------------------------------------------------------------------------
/// Removes loaded filters to prepare the tag for a new use.
// ---------------------------------------------------------------------------
static void
box_free (void)
{
  boxesN = 0;
  free (boxes);
  boxes = NULL;
}
