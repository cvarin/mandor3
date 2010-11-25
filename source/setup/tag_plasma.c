/** \file tag_plasma.c
  * Creation of arbitrary shaped plasma.
  *
  * Filtering problem is tightly related to the creation of the arbitrary
  * shaped figure composed from many primitives. The problem is very nicely
  * described in SIGGRAPH paper on drawing constructions and involves parsing
  * of big boolean expressions like "(fig1 + fig2)*fig3 - fig4", where '+'
  * stays for 'or', '*' means 'and' and so on. I do not implement this
  * completely because of engine already reproduces all old primitives and
  * adds many new.
  *
  * First of all I'll rewrite all old tags (basically, remove all shaping and
  * sampling options from there because of it is done here already).
  *
  * Than I may compile object file, produse .so module and load it dynamically
  * to produce patch simular to my C_pathc in Python. It will permit to pose
  * any function of coordinates and velocities to be used to set density
  * distribution.
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "type_vector.h"

#include "log.h"
#include "misc_cfgReader.h"

#include "main.h"
#include "setup/plasma.h"

/// Parameters of the last plasma object we have allocated using tag_plasma.
int plasma_nx, plasma_ny, plasma_nz, plasma_layers, plasma_PPC;

#include "tag_plasma_box.c"		// Polyhedrons filters.
#include "tag_plasma_convex.c"		// Polyhedrons filters.
#include "tag_plasma_spheres.c"		// Spherical regions filters (balls, clusters, cavities, bumps on the surface).
#include "tag_plasma_cylinders.c"	// Cylinders regions filters (channels, wires, scratches on the surface).

// ---------------------------------------------------------------------------
/// Creates piece of plasma (velocity distribution to be set later).
// ---------------------------------------------------------------------------
double
tag_plasma (FILE *fp)
{
   // Cell sub-sampling and velocity space resolution.
   plasma_nx     = cfg_readInt (fp);
   plasma_ny     = cfg_readInt (fp);
   plasma_nz     = cfg_readInt (fp);
   plasma_layers = cfg_readInt (fp);

   // Reads all exctrusion elements.
   while (cfg_isParameter (fp)) {
      const char *element = cfg_readWord (fp);
      if (!strcmp (element, "planes")) {
         convex_read (fp);			// Reads set of planes.
      } else if (!strcmp (element, "spheres")) {
         sphere_read (fp);			// Reads set of spheres.
      } else if (!strcmp (element, "cylinders")) {
         cylinder_read (fp);			// Reads set of cylinders.
      } else if (!strcmp (element, "boxes")) {
         box_read (fp);				// Reads set of boxes.
      } else {
         DIE ("bad shape '%s' ('planes', 'spheres', or 'cylinders' expected)",
              element);
      }
   }

   // Deactivates axises if necessary.
   plasma_nx  = mc_have_x ? plasma_nx : 1;
   plasma_ny  = mc_have_y ? plasma_ny : 1;
   plasma_nz  = mc_have_z ? plasma_nz : 1;
   plasma_PPC = plasma_nx*plasma_ny*plasma_nz*plasma_layers;

   // Checks parameters.
   ENSURE (plasma_nx > 0 && plasma_ny > 0 && plasma_nz > 0 &&
           plasma_layers > 0,
           "bad nx(%d), ny(%d), nz(%d) or number of velocity layers (%d)",
           plasma_nx, plasma_ny, plasma_nz, plasma_layers);

   say ("tag_plasma: ");
   say ("  - %d x %d x %d in-cell spacing", plasma_nx, plasma_ny, plasma_nz);
   say ("  - %d velocity layers", plasma_layers);
   say ("  - %d particles per cell", plasma_PPC);

   marker_t cell[plasma_PPC], filtered[plasma_PPC];
   double dx = h1*mc_have_x/(double) plasma_nx,
          dy = h2*mc_have_y/(double) plasma_ny,
          dz = h3*mc_have_z/(double) plasma_nz;

   // 'p' enumerates particles in cell.
   int p = 0;
   for (int i = 0 ; i < plasma_nx ; ++i)
   for (int j = 0 ; j < plasma_ny ; ++j)
   for (int k = 0 ; k < plasma_nz ; ++k) {
      // Fills sample (NAN) to catch unitialized variables.
      marker_t m = {.x = (i + 0.5)*dx,
                    .y = (j + 0.5)*dy,
                    .z = (k + 0.5)*dz,
                    .vx = NAN,
                    .vy = NAN,
                    .vz = NAN,
                    .rho = NAN,
                    .qDivM = NAN};
      for (int l = 0 ; l < plasma_layers ; ++l) {
         cell[p]    = m;
         cell[p].vx = p++;	// Stores enumeration position.
      }
   }

   double memUsage = 0;
   plasma_newObject ();
   for (int I = cpu_min[0] ; I < cpu_max[0] + 1 - mc_have_x ; ++I)
   for (int J = cpu_min[1] ; J < cpu_max[1] + 1 - mc_have_y ; ++J)
   for (int K = cpu_min[2] ; K < cpu_max[2] + 1 - mc_have_z ; ++K) {
      // Adds cell data into current cell.
      for (int p = 0 ; p < plasma_PPC ; ++p) {
         filtered[p] = cell[p];
         filtered[p].x += I*h1*mc_have_x;
         filtered[p].y += J*h2*mc_have_y;
         filtered[p].z += K*h3*mc_have_z;
      }

      // Set of filter extrusions.
      int PPC = plasma_PPC;
      for (convex_t *b = polys, *b2 = b + polysN ; b < b2 ; ++b) {
         PPC = convex_filter (filtered, PPC, b);
      }

      for (sphere_t *s = spheres, *s2 = s + spheresN ; s < s2 ; ++s) {
         PPC = sphere_filter (filtered, PPC, s);
      }

      for (cylinder_t *c = cylinders, *c2 = c + cylindersN ; c < c2 ; ++c) {
         PPC = cylinder_filter (filtered, PPC, c);
      }

      for (box_t *b = boxes, *b2 = b + boxesN ; b < b2 ; ++b) {
         PPC = box_filter (filtered, PPC, b);
      }

      // Adds particles to the storage.
      for (int p = 0 ; p < PPC && (!memEstimateOnly) ; ++p) {
         *plasma_marker () = filtered[p];
      }

      memUsage += PPC;
   }

   say ("  - %.3e particles on cpu %d", memUsage, cpu_here);

   // Turns off all filters.
   box_free      ();
   convex_free   ();
   cylinder_free ();
   sphere_free   ();

   return memUsage*sizeof (marker_t);						// Returns memory usage.
}
