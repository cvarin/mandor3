/** \file tag_plasma_spheres.c
  * Helpers to read and apply filters composed from many spherical cavities or
  * balls. File is included directly into tag_plasma.c and is separated for
  * readability only.
  */

// ---------------------------------------------------------------------------
/// Sphere in space:
///   if 'r' > 0 than all particles inside  of sphere are removed
///   if 'r' < 0 than all particles outside of sphere are removed.
// ---------------------------------------------------------------------------
typedef struct {
   double  x, y, z, r;
} sphere_t;

static int       spheresN = 0;		///< Number of spheres.
static sphere_t *spheres  = NULL;	///< Spheres.

// ---------------------------------------------------------------------------
/// Reads parameters of spheres in uniformly indented block of config file.
// ---------------------------------------------------------------------------
static void
sphere_read (FILE *fp)
{
   // Scans all '>' indented data blocks.
   while (cfg_isOption (fp)) {
      spheres = (sphere_t*) realloc (spheres, (++spheresN)*sizeof (sphere_t));
      assert (spheres);
      sphere_t *s = spheres + spheresN - 1;
      int loaded_params = fscanf (fp, "> %le %le %le %le %*[^\n] ",
                                  &s->x, &s->y, &s->z, &s->r);
      ENSURE (4 == loaded_params, "cannot read all parameters of the sphere");
   }
}

// ---------------------------------------------------------------------------
/// Takes an array with markers and keeps only markers inside of the region.
// ---------------------------------------------------------------------------
static int
sphere_filter (marker_t *p, int N, sphere_t *s)
{
   const double r2 = s->r*s->r;
   for (marker_t *end = p + N ; p < end ; ++p) {
      double d2 = (p->x - s->x)*(p->x - s->x) +
                  (p->y - s->y)*(p->y - s->y) +
                  (p->z - s->z)*(p->z - s->z);

      // Sign of radius chooses which part to keep.
      if ((r2 - d2)*s->r < 0) {
         *p = *(--end);		// Replaces deleted particles by untested one.
         --p;			// Unrolls counter to test the new marker.
         --N;
      }
   }
   return N;
}

// ---------------------------------------------------------------------------
/// Removes loaded filters to reset the tag for next use.
// ---------------------------------------------------------------------------
static void
sphere_free (void)
{
   free (spheres);
   spheresN = 0;
   spheres  = NULL;
}
