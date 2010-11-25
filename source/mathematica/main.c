/** \file main.c
  * Converts binary file to the txt-form for Mathematica.
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <assert.h>

#include <unistd.h>
#include <dirent.h>

#include "type_mesh.h"

#include "log.h"
#include "misc_units.h"
#include "misc_partition.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "IO_tecplot.h"

static meshDouble_t rho;
static meshVec_t    E, H, J;

static int  rec_start, rec_end, rec_step;

// Choice of output units: 0 <=> r0, t0;
//                         1 <=> microns, femtoseconds.
static int    use_um_and_fs = -1;

static reg_t  subdomain = {
   .min = {INT_MIN, INT_MIN, INT_MIN},
   .max = {INT_MAX, INT_MAX, INT_MAX},
   .cpu = 0,
};

// ---------------------------------------------------------------------------
/// Checks if domain parameters were updated and meshes should be resized.
// ---------------------------------------------------------------------------
static void
main_resize_meshes (const reg_t *reg)
{
   static reg_t old = mc_reg_initBad;

   // Returns if regions are the same.
   if (reg_isInside(reg, &old)
   &&  reg_isInside(&old, reg)) {
      return;
   }

   old = *reg;
   mf_reg_collapse (&old);

   mesh_resize (mcast_mesh (&E),   old.min[0], old.min[1], old.min[2],
                                   old.max[0], old.max[1], old.max[2]);
   mesh_resize (mcast_mesh (&H),   old.min[0], old.min[1], old.min[2],
                                   old.max[0], old.max[1], old.max[2]);
   mesh_resize (mcast_mesh (&J),   old.min[0], old.min[1], old.min[2],
                                   old.max[0], old.max[1], old.max[2]);
   mesh_resize (mcast_mesh (&rho), old.min[0], old.min[1], old.min[2],
                                   old.max[0], old.max[1], old.max[2]);

   // Fills all nodes with QNaNs to catch uninitialized boundary nodes.
   memset (E.storage,   -1, E.size);
   memset (H.storage,   -1, H.size);
   memset (J.storage,   -1, J.size);
   memset (rho.storage, -1, rho.size);

   partition_init ();
}

// ---------------------------------------------------------------------------
/// Reads configuration file and returns title for the zone(s).
// ---------------------------------------------------------------------------
static void
main_init_all (int argc, char *argv[])
{
   parameter_enterMPI (argc, argv, 0);
   parameter_load ();

   mesh_allocate (mcast_mesh (&E),   0, 0, 0, 0, 0, 0, "E", mc_vec3D_t);
   mesh_allocate (mcast_mesh (&H),   0, 0, 0, 0, 0, 0, "H", mc_vec3D_t);
   mesh_allocate (mcast_mesh (&J),   0, 0, 0, 0, 0, 0, "J", mc_vec3D_t);
   mesh_allocate (mcast_mesh (&rho), 0, 0, 0, 0, 0, 0, "ρ", mc_double);

   const int lastMeshNum = tecIO_recordsTotal () - 1;
   units_load ();

   FILE *fp  = cfg_open ("diag_math.cfg", "rt", __func__);
   rec_start = cfg_readInt (fp);
   rec_end   = cfg_readInt (fp);
   rec_step  = cfg_readInt (fp);

   // Reads and clamps to use as array index.
   use_um_and_fs = cfg_readInt (fp) == 1;

   // Gets optional subdomain.
   if (cfg_isOption (fp)) {
      subdomain.min[0] = cfg_readOptInt (fp);
      subdomain.max[0] = cfg_readOptInt (fp);
      subdomain.min[1] = cfg_readOptInt (fp);
      subdomain.max[1] = cfg_readOptInt (fp);
      subdomain.min[2] = cfg_readOptInt (fp);
      subdomain.max[2] = cfg_readOptInt (fp);
   }
   fclose (fp);

   mf_reg_collapse (&subdomain);

   if (rec_end < 0 || rec_end > lastMeshNum) {
      rec_end = lastMeshNum;
   }

   ENSURE (rec_start >= 0 && rec_end >= rec_start && rec_step > 0,
           "bad set of record's (from %d to %d with step %d)",
           rec_start, rec_end, rec_step);
}

// Protection against predefined IBM stuff under AIX.
#undef hz

// ---------------------------------------------------------------------------
/// Removes old files.
// ---------------------------------------------------------------------------
static void
main_clear_folder (void)
{
   DIR           *folder;
   struct dirent *item;

   int error = chdir ("./output/fields");
   assert (!error);

   int  n = 0;
   folder = opendir (".");
   while ((item = readdir (folder)) != NULL) {
      remove (item->d_name);
      ++n;
   }
   closedir (folder);
   chdir ("../..");

   say ("\nOld files are erased.\n");
}


// ---------------------------------------------------------------------------
/// Diagnostic entry point.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
   main_init_all (argc, argv);
   main_clear_folder ();

   // Chooses units for output.
   double A0     = units (mc_A0)/units (mc_E0),
          rho_cr = units (mc_ne_critical)*mc_CGS_e/units (mc_rho0),
          r0     = 1.0,
          t0     = 1.0;
   if (use_um_and_fs) {
      t0 = units (mc_femtosecond);
      r0 = units (mc_micron);
      say ("[t] = 'fs', [r] = 'μm'");
   } else {
      say ("[t] = 't0' = λ/c, [r] = 'r0' = λ");
   }
   say ("[E] = [H] = 'A0', [ρ] = '|e⋅ne_crit|', [J] = [a.u]\n");

   ENSURE (cpu_total == 1, "diagnostic isn't parallel");

   for (int record = rec_start ; record <= rec_end ; record += rec_step) {
      reg_t reg;
      reg.min[0] = fmax (dmn_min[0], subdomain.min[0]) + 0.1;
      reg.min[1] = fmax (dmn_min[1], subdomain.min[1]) + 0.1;
      reg.min[2] = fmax (dmn_min[2], subdomain.min[2]) + 0.1;
      reg.max[0] = fmin (dmn_max[0], subdomain.max[0]) + 0.1;
      reg.max[1] = fmin (dmn_max[1], subdomain.max[1]) + 0.1;
      reg.max[2] = fmin (dmn_max[2], subdomain.max[2]) + 0.1;

      main_resize_meshes (&reg);

      double time;
      tecIO_load (record, &time, &reg, &E, &H, &J, &rho);
      say_doing ("rec: %d, t: %.3f, mesh: %d x %d x %d",
                 record,
                 time,
                 reg.max[0] - reg.min[0] + 1,
                 reg.max[1] - reg.min[1] + 1,
                 reg.max[2] - reg.min[2] + 1);

      FILE *fp = cfg_open (_("output/fields/rec_%08d.dat", record),
                           "wt",
                           "math.out");
      fprintf (fp, "x y z Ex Ey Ez E² Hx Hy Hz H² jx jy jz ρ\n");
      for (int i = reg.min[0] ; i <= reg.max[0] ; ++i)
      for (int j = reg.min[1] ; j <= reg.max[1] ; ++j)
      for (int k = reg.min[2] ; k <= reg.max[2] ; ++k) {
         vec3D_t e = mv_v(&E,   i, j, k),
                 h = mv_v(&H,   i, j, k),
                 c = mv_v(&J,   i, j, k);
         double  r = mv_f(&rho, i, j, k);
         fprintf (fp, "%le %le %le ", i*h1/r0, j*h2/r0, k*h3/r0);
         fprintf (fp, "% le % le % le %le ", e.x/A0, e.y/A0, e.z/A0,
                  (e.x*e.x + e.y*e.y + e.z*e.z)/(A0*A0));
         fprintf (fp, "% le % le % le %le ", h.x/A0, h.y/A0, h.z/A0,
                  (h.x*h.x + h.y*h.y + h.z*h.z)/(A0*A0));
         fprintf (fp, "% le % le % le ", c.x, c.y, c.z);
         fprintf (fp, "% le\n", r/rho_cr);
      }
      fclose (fp);
   }

   say ("\nDone.");

   return EXIT_SUCCESS;
}
