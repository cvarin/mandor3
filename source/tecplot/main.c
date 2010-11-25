/** \file tecplot.c
  * Exporter of the binary files to the tecplot format. Supports slices and subdomains to reduce data size.
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <float.h>

#include <unistd.h>			// To have 'chdir'.
#include <dirent.h>			// To do folder scanning.

#include "type_mesh.h"
#include "type_CFile.h"

#include "log.h"
#include "misc_units.h"
#include "misc_partition.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "IO_tecplot.h"

meshDouble_t rho;			///< Charge density.
meshVec_t    E;				///< Vector fields.
int 	     recStart, recEnd, recStep;	///< Set of records to process.

int    use_microns_and_fs = -1;		///< Output units: r₀/t₀ or μm/fs.
double uFemtosecond,			///< Dimensionless units: fs/t₀,
       uRhoCrit,			///< 			  ρ_crit/ρ₀,
       uMicron,				///< 			  μm/r₀,
       uA0;				///< 			  A₀/E₀,

/// Clipping subdomain.
enum  { eUnknown = -1, eDomain = 0, eSlice, eSubdomain } tecplotRegType;
static reg_t  subdomain = {
   .min = {INT_MIN, INT_MIN, INT_MIN},
   .max = {INT_MAX, INT_MAX, INT_MAX},
   .cpu = 0,
};

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// ---------------------------------------------------------------------------
/// Checks if domain parameters were updated and meshes should be resized.
// ---------------------------------------------------------------------------
static void
tecplot_reconfigureMeshes (const reg_t *reg)
{
   // Flag to check if parameters were updated.
   static reg_t old = mc_reg_initBad;

   // Returns if regions are the same.
   if (reg_isInside(reg, &old)
   &&  reg_isInside(&old, reg)) {
      return;
   }

   old = *reg;			// Makes copy.
   mf_reg_collapse (&old);	// Deactivates axises.

   mesh_resize (mcast_mesh (&E),   old.min[0], old.min[1], old.min[2],
                                   old.max[0], old.max[1], old.max[2]);
   mesh_resize (mcast_mesh (&rho), old.min[0], old.min[1], old.min[2],
                                   old.max[0], old.max[1], old.max[2]);

   // Fills all nodes with QNaNs to catch uninitialized boundary nodes.
   memset (E.storage,   -1, E.size);
   memset (rho.storage, -1, rho.size);

   partition_init ();
}

// ---------------------------------------------------------------------------
/// Reads configuration file and returns title for the zone(s).
// ---------------------------------------------------------------------------
static void
tecplot_prepareEnvironment (int argc, char *argv[])
{
   // Initializes log-file and MPI subsystem.
   parameter_enterMPI (argc, argv, 0);

   // Allocates empty meshes.
   mesh_allocate (mcast_mesh (&E),   0, 0, 0, 0, 0, 0, "tecplot_vector", mc_vec3D_t);
   mesh_allocate (mcast_mesh (&rho), 0, 0, 0, 0, 0, 0, "tecplot_scalar", mc_double);

   const int lastMeshNum = tecIO_recordsTotal () - 1;
   units_load ();

   // Reads which records to export.
   FILE *fp = cfg_open ("diag_tecplot.cfg", "rt", __func__);
   recStart = cfg_readInt (fp);
   recEnd   = cfg_readInt (fp);
   recStep  = cfg_readInt (fp);

   // Reads flag and normalizes it to use as array index.
   use_microns_and_fs = (cfg_readInt (fp) == 1);

   const char *regTypeName = cfg_readWord (fp);
   tecplotRegType = cfg_identifyWord (regTypeName, "domain",    eDomain,
                                                   "slice",     eSlice,
                                                   "subdomain", eSubdomain,
                                                   mc_cfgTermGuesses);
   switch (tecplotRegType) {
   case -1:
   default:
      DIE ("bag diag_tecplot.cfg: unknown parameter '%s' (expect 'domain', "
           "'slice' or 'subdomain')", regTypeName);
      break;

   case eDomain:													// Domain is domain, no additional info needed.
      break;

   case eSlice:
      {
         int axis = cfg_readOptInt (fp);
         ENSURE(axis >= 0 && axis < 3,
                "bag axis: use 0(x), 1(y) or 2(z), not '%d'", axis);
         subdomain.min[axis] = subdomain.max[axis] = cfg_readOptInt (fp);
      }
      break;

   case eSubdomain:
      subdomain.min[0] = cfg_readOptInt (fp);
      subdomain.min[1] = cfg_readOptInt (fp);
      subdomain.min[2] = cfg_readOptInt (fp);
      subdomain.max[0] = cfg_readOptInt (fp);
      subdomain.max[1] = cfg_readOptInt (fp);
      subdomain.max[2] = cfg_readOptInt (fp);
      break;
   }
   fclose (fp);

   mf_reg_collapse (&subdomain);

   // Sets range of records (max or given).
   if (recEnd < 0 || recEnd > lastMeshNum) {
      recEnd = lastMeshNum;
   }

   ENSURE (recStart >= 0 && recEnd >= recStart && recStep > 0,
           "bad record set: start = %d, end = %d, step = %d",
           recStart, recStep, recEnd);

   // Sets output units.
   uA0 = units (mc_A0)/units (mc_E0);
   if (use_microns_and_fs) {
      uFemtosecond = units (mc_femtosecond);
      uMicron = units (mc_micron);
   } else {
      uFemtosecond = uMicron = 1;
   }
   uRhoCrit = units (mc_ne_critical)*mc_CGS_e/units (mc_rho0);
}

// Protection from predefined IBM's stuff on AIX.
#undef hz

// ---------------------------------------------------------------------------
/// Outputs coordinate field in block format.
// ---------------------------------------------------------------------------
static void
export_coordinate (FILE *fp, const reg_t *reg,
                   const double hx, const double hy, const double hz)
{
   for (int k = reg->min[2] ; k <= reg->max[2] ; ++k)
   for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
   for (int i = reg->min[0] ; i <= reg->max[0] ; ++i) {
      fprintf (fp, "%.5e\n", i*hx + j*hy + k*hz);
   }
}

// ---------------------------------------------------------------------------
/// Outputs scalar field in block format.
// ---------------------------------------------------------------------------
static void
export_scalar (FILE *fp, const reg_t *reg, meshDouble_p rho,
               const double scale, double *min, double *max)
{
   double Min = DBL_MAX,
          Max = DBL_MIN;
   for (int k = reg->min[2] ; k <= reg->max[2] ; ++k)
   for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
   for (int i = reg->min[0] ; i <= reg->max[0] ; ++i) {
      double val = mv_f(rho, i, j, k)*scale;
      fprintf (fp, "%.3e\n", val);
      Min = MIN(val, Min);
      Max = MAX(val, Max);
   }
   *min = Min;
   *max = Max;
}

// ---------------------------------------------------------------------------
/// Outputs component of vector field in block format.
// ---------------------------------------------------------------------------
static void
export_vecComponent (FILE *fp, const reg_t *reg, meshVecI_p E,
                     const int component, const double scale,
                     double *min, double *max)
{
   double Min = DBL_MAX,
          Max = DBL_MIN;
   for (int k = reg->min[2] ; k <= reg->max[2] ; ++k)
   for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
   for (int i = reg->min[0] ; i <= reg->max[0] ; ++i) {
      double val = mv_fi(E, i, j, k, component)*scale;
      fprintf (fp, "%.3e\n", val);
      Min = MIN(val, Min);
      Max = MAX(val, Max);
   }
   *min = Min;
   *max = Max;
}

// ---------------------------------------------------------------------------
/// Outputs squared absolute value of vector field in block format.
// ---------------------------------------------------------------------------
static void
export_vecSquare (FILE *fp, const reg_t *reg, meshVec_p E, const double scale, double *min, double *max)
{
   double Min = DBL_MAX,
          Max = DBL_MIN;
   for (int k = reg->min[2] ; k <= reg->max[2] ; ++k)
   for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
   for (int i = reg->min[0] ; i <= reg->max[0] ; ++i) {
      double val = scale*(mv_fx(E, i, j, k)*mv_fx(E, i, j, k)
                        + mv_fy(E, i, j, k)*mv_fy(E, i, j, k)
                        + mv_fz(E, i, j, k)*mv_fz(E, i, j, k));
      fprintf (fp, "%.3e\n", val);
      Min = MIN(val, Min);
      Max = MAX(val, Max);
   }
   *min = Min;
   *max = Max;
}

// ---------------------------------------------------------------------------
/// Removes old files.
// ---------------------------------------------------------------------------
static void
main_clearFolder (void)
{
   DIR           *folder;
   struct dirent *item;

   const int found = chdir ("./output/tecplot");	// >>>
   assert (found == 0 && "cannot chdir to './output/tecplot'");

   // Opens filelist and removes all files.
   int n = 0;
   folder = opendir (".");
   while ((item = readdir (folder)) != NULL) {
      remove (item->d_name);
      ++n;
   }
   closedir (folder);
   chdir ("../..");					// <<<

   if (n) say ("tecplot.out: all exported earlier files are erased.\n");
}

// ---------------------------------------------------------------------------
/// Checks that older files are saved with proper units.
// ---------------------------------------------------------------------------
static void
main_checksOldUnits (void)
{
   assert (use_microns_and_fs >= 0 && "input parameters are uninitialized");

   CFile_t *registry = CF_openUpdate ("output/tecplot/registry.txt");
   int unitsAreDifferent = 1;

   // Checks if units are changed.
   if (!CF_findChunk (registry, ">>> units frame")) {
      int oldUnits;
      CF_scan (registry, "%d", &oldUnits);
      if (oldUnits == use_microns_and_fs) {
         unitsAreDifferent = 0;
      }
   }

   if (unitsAreDifferent) {
      // Switches to new units.
      CF_close (registry);

      main_clearFolder ();

      registry = CF_openWrite ("output/tecplot/registry.txt");
      CF_openChunk (registry, ">>> units frame");
      CF_print (registry, "%d\tUnits system:\n"
                          "\to 0 => numerical, r_0/t_0\n"
                          "\to 1 => laser, microns/fs\n",
                          (use_microns_and_fs > 0));
      CF_closeChunk (registry);
      CF_close (registry);
   }
}

// ---------------------------------------------------------------------------
/// Diagnostic entry point.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
   // Prepares all environment frames.
   tecplot_prepareEnvironment (argc, argv);

   ENSURE (cpu_total == 1, "tecplot out is not yet parallel");

   // Removes old registry if asked.
   for (int i = 1 ; i < argc ; ++i) {
      if (!strcmp (argv[i], "--new-export") || !strcmp (argv[i], "-n")) {
         main_clearFolder ();
      }
      if (!strcmp (argv[i], "--help") || !strcmp (argv[i], "-h")) {
         printf (
"Usage: tecplot.out -n[--new-export] -h[--help]\n\n"
"Notes:\n"
"    * input parameters are taken from 'diag_tecplot.cfg'. \n"
"    * output files are in './output/tecplot' folder.\n\n"
"WARNING: output files are big, so if a record was once exported,\n"
"         old file is simply kept untouched. Use -n (--new_export) \n"
"         option to overwrite the output file.\n\n"
"Normally scripts vTec1D.sh, vTec2D.sh and vTec3D.sh take care \n"
"about opening proper layout and generating menu for data browser.\n");
         exit (0);
      }
   }

   say ("Start to generate tecplot data file...");

   main_checksOldUnits ();
   CFile_t *registry = CF_openUpdate ("output/tecplot/registry.txt");

   static const int output1D = (mc_have_x + mc_have_y + mc_have_z == 1);
   for (int record = recStart ; record <= recEnd ; record += recStep) {
      parameter_load ();

      // Packs loaded parameters of the domain.
      reg_t reg;
      reg.min[0] = MAX(dmn_min[0], subdomain.min[0]);
      reg.min[1] = MAX(dmn_min[1], subdomain.min[1]);
      reg.min[2] = MAX(dmn_min[2], subdomain.min[2]);
      reg.max[0] = MIN(dmn_max[0], subdomain.max[0]);
      reg.max[1] = MIN(dmn_max[1], subdomain.max[1]);
      reg.max[2] = MIN(dmn_max[2], subdomain.max[2]);

      // Fits meshes to the domain and updates everything.
      tecplot_reconfigureMeshes (&reg);

      // Checks if record was created already.
      if (!CF_findChunk (registry, ">>> record %08d", record)) {
         reg_t tmp;

         // Simple test if region is the same.
         CF_scan (registry, "%d%*[^\n]", &tmp.min[0]);
         CF_scan (registry, "%d%*[^\n]", &tmp.min[1]);
         CF_scan (registry, "%d%*[^\n]", &tmp.min[2]);
         CF_scan (registry, "%d%*[^\n]", &tmp.max[0]);
         CF_scan (registry, "%d%*[^\n]", &tmp.max[1]);
         CF_scan (registry, "%d%*[^\n]", &tmp.max[2]);
         if (reg_isInside (&reg, &tmp)
         &&  reg_isInside (&tmp, &reg)) {
            continue;
         }

         CF_deleteChunk (registry, ">>> record %08d", record);
      }

      char _unitSpace[2][32] = {"[<greek>l</greek><sub>0</sub>]",
                                "[<greek>m</greek>m]"},
           _unitTime [2][20] = {"[t<sub>0</sub>]",
                                "[fs]"};
      const char * const cSpace  = _unitSpace[(use_microns_and_fs > 0)];
      const char * const cTime   = _unitTime[(use_microns_and_fs > 0)];
      const char * const cField  = "[a<sub>0</sub>]";
      const char * const cField2 = "[a<sub>0</sub><sup>2</sup>]";
      const char * const cRho    = "[|e|n<sub>e crit</sub>]";

      // Loads time and E.
      double time;
      tecIO_load (record, &time, &reg, &E, NULL, NULL, NULL);

      // Opens output file (for 1D case we append to the same file).
      FILE *fp = cfg_open (_("output/tecplot/tecplot_%08d.dat", record*(!output1D)),
                           (!output1D || record == recStart) ? "wt"
                                                             : "at",
                           "tecplot.out");
      fprintf (fp, "variables = \"x %s\", \"y %s\", \"z %s\", "
                   "\"E<sub>x</sub> %s\", "
                   "\"E<sub>y</sub> %s\", "
                   "\"E<sub>z</sub> %s\", "
                   "\"E<sup>2</sup> %s\", "
                   "\"H<sub>x</sub> %s\", "
                   "\"H<sub>y</sub> %s\", "
                   "\"H<sub>z</sub> %s\", "
                   "\"H<sup>2</sup> %s\", "
                   "\"j<sub>x</sub> [a.u.]\", "
                   "\"j<sub>y</sub> [a.u.]\", "
                   "\"j<sub>z</sub> [a.u.]\", "
                   "\"<greek>r</greek> %s\"\n",
                   cSpace, cSpace, cSpace,
                   cField, cField, cField, cField2,
                   cField, cField, cField, cField2,
                   cRho);

      fprintf (fp, "zone t=\"Rec %d (%.3f%s)\", "
                   "i = %d, j = %d, k = %d, f = block\n",
                   record, time/uFemtosecond, cTime,
                   reg.max[0] - reg.min[0] + 1,
                   reg.max[1] - reg.min[1] + 1,
                   reg.max[2] - reg.min[2] + 1);

      export_coordinate (fp, &reg, h1/uMicron, 0, 0);
      export_coordinate (fp, &reg, 0, h2/uMicron, 0);
      export_coordinate (fp, &reg, 0, 0, h3/uMicron);
      say_doing ("record #%d, time = %f, mesh exported ...", record, time);

      double minEx, minEy, minEz, minE2,
             maxEx, maxEy, maxEz, maxE2;
      export_vecComponent (fp, &reg, mcast_meshVecI (&E), mc_x, 1.0/uA0, &minEx, &maxEx);
      export_vecComponent (fp, &reg, mcast_meshVecI (&E), mc_y, 1.0/uA0, &minEy, &maxEy);
      export_vecComponent (fp, &reg, mcast_meshVecI (&E), mc_z, 1.0/uA0, &minEz, &maxEz);
      export_vecSquare (fp, &reg, &E, 1.0/(uA0*uA0), &minE2, &maxE2);
      say_doing ("record #%d, time = %f, E exported ...", record, time);

      // Loads H.
      double minHx, minHy, minHz, minH2, maxHx, maxHy, maxHz, maxH2;
      tecIO_load (record, &time, &reg, NULL, &E, NULL, NULL);
      export_vecComponent (fp, &reg, mcast_meshVecI (&E), mc_x, 1.0/uA0, &minHx, &maxHx);
      export_vecComponent (fp, &reg, mcast_meshVecI (&E), mc_y, 1.0/uA0, &minHy, &maxHy);
      export_vecComponent (fp, &reg, mcast_meshVecI (&E), mc_z, 1.0/uA0, &minHz, &maxHz);
      export_vecSquare (fp, &reg, &E, 1.0/(uA0*uA0), &minH2, &maxH2);
      say_doing ("record #%d, time = %f, H exported ...", record, time);

      // Loads J.
      double minJx, minJy, minJz, minRho, maxJx, maxJy, maxJz, maxRho;
      tecIO_load (record, &time, &reg, NULL, NULL, &E, NULL);
      export_vecComponent (fp, &reg, mcast_meshVecI (&E), mc_x, 1.0/uA0, &minJx, &maxJx);
      export_vecComponent (fp, &reg, mcast_meshVecI (&E), mc_y, 1.0/uA0, &minJy, &maxJy);
      export_vecComponent (fp, &reg, mcast_meshVecI (&E), mc_z, 1.0/uA0, &minJz, &maxJz);

      // Loads rho.
      tecIO_load (record, &time, &reg, NULL, NULL, NULL, &rho);
      export_scalar (fp, &reg, &rho, 1.0/uRhoCrit, &minRho, &maxRho);
      fclose (fp);

      // Updates registry.
      CF_openChunk (registry, ">>> record %08d", record);
      CF_print (registry, "%d   \ti-min\n", reg.min[0]);
      CF_print (registry, "%d   \tj-min\n", reg.min[1]);
      CF_print (registry, "%d   \tk-min\n", reg.min[2]);
      CF_print (registry, "%d   \ti-max\n", reg.max[0]);
      CF_print (registry, "%d   \tj-max\n", reg.max[1]);
      CF_print (registry, "%d   \tk-max\n", reg.max[2]);
      CF_print (registry, "[%+.4e, %+.4e]   \tEx range\n", minEx, maxEx);
      CF_print (registry, "[%+.4e, %+.4e]   \tEy range\n", minEy, maxEy);
      CF_print (registry, "[%+.4e, %+.4e]   \tEz range\n", minEz, maxEz);
      CF_print (registry, "[%+.4e, %+.4e]   \tE2 range\n", minE2, maxE2);
      CF_print (registry, "[%+.4e, %+.4e]   \tHx range\n", minHx, maxHx);
      CF_print (registry, "[%+.4e, %+.4e]   \tHy range\n", minHy, maxHy);
      CF_print (registry, "[%+.4e, %+.4e]   \tHz range\n", minHz, maxHz);
      CF_print (registry, "[%+.4e, %+.4e]   \tH2 range\n", minH2, maxH2);
      CF_print (registry, "[%+.4e, %+.4e]   \tjx range\n", minJx, maxJx);
      CF_print (registry, "[%+.4e, %+.4e]   \tjy range\n", minJy, maxJy);
      CF_print (registry, "[%+.4e, %+.4e]   \tjz range\n", minJz, maxJz);
      CF_print (registry, "[%+.4e, %+.4e]   \trho range\n", minRho, maxRho);
      CF_closeChunk (registry);

      say_doing ("record #%d, time = %f is ready", record, time);
   }

   registry = CF_defrag (registry);
   CF_close (registry);

   say ("Done.");

   return EXIT_SUCCESS;
}
