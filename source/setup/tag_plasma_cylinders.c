/** \file tag_plasma_cylinders.c
  * Helpers to read and apply filters composed from many cylindrical channels
  * or wires. File is included directly into tag_plasma.c and is separated for
  * readability only.
  */

// ---------------------------------------------------------------------------
/// Cylinder in space:
///   if 'r' > 0 than all particles inside  of cylinder are removed;
///   if 'r' < 0 than all particles outside of cylinder are removed.
// ---------------------------------------------------------------------------
typedef struct {
   vec3d_t  pos, axis;
   double r;
} cylinder_t;

static int         cylindersN = 0;		///< Number of cylinders.
static cylinder_t *cylinders  = NULL;		///< Cylinders.

// ---------------------------------------------------------------------------
/// Reads parameters of cylinders in '>' indented block of config file.
// ---------------------------------------------------------------------------
static void
cylinder_read (FILE *fp)
{
   // Scans all '>' indented data blocks.
   while (cfg_isOption (fp)) {
      // Adds new cylinder.
      cylinders = (cylinder_t*) realloc (cylinders,
                                         (++cylindersN)*sizeof (cylinder_t));
      assert (cylinders);
      // Reads new cylinder.
      cylinder_t *c = cylinders + cylindersN - 1;
      ENSURE (3 == fscanf (fp, "> %le %le %le %*[^\n] ",
                           &c->pos.v.x, &c->pos.v.y, &c->pos.v.z),
              "cannot read axis position of cylinder");
      ENSURE (3 == fscanf (fp, "> %le %le %le %*[^\n] ",
                           &c->axis.v.x, &c->axis.v.y, &c->axis.v.z),
              "cannot read axis direction of cylinder");
      ENSURE (1 == fscanf (fp, "> %le %*[^\n] ", &c->r),
              "cannot read radius of cylinder");

      // Normalizes an amplitude of the axis vector.
      ENSURE (vd_dot (c->axis, c->axis) > 0.01,
              "too small axis vector(%s)", vd_print (c->axis));
      c->axis = vd_scale (c->axis, 1.0/sqrt (vd_dot (c->axis, c->axis)));
   }
}

// ---------------------------------------------------------------------------
/// Takes an array with markers and keeps only markers inside of the bounded region.
// ---------------------------------------------------------------------------
static int
cylinder_filter (marker_t *p, int N, cylinder_t *c)
{
   const double r2 = c->r*c->r;
   for (marker_t *end = p + N ; p < end ; ++p) {
      vec3d_t diff = vd_sub (vd (p->x, p->y, p->z), c->pos);
      diff = vd_cross (diff, c->axis);

      // Sign of radius chooses which part to keep.
      if ((r2 - vd_dot (diff, diff))*c->r < 0) {
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
cylinder_free (void)
{
   free (cylinders);
   cylindersN = 0;
   cylinders  = NULL;
}
