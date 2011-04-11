/** \file main.c
  * Setup - creation of initial state using multi-pass approach and external
  * config file to set the sequence of passes.
  */

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>		// For folder scanning.

#define mc_activateProfiler 1	///< \internal Activates profiler
#define setup_main_header	///< \internal Prevents setup/main.h from being included.

#include <mpi.h>
#include <math.h>

#include "type_marker.h"

#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_partition.h"
#include "misc_parameters.h"
#include "misc_markerPlacer.h"

#include "log.h"
#include "commandLine.h"

#include "IO_sys.h"				// Disk IO.
#include "IO_tecplot.h"
#include "IO_names.h"
#include "IO_fileMap.h"

#include "setup/tag.h"
#include "setup/plasma.h"

/// Enum for memory consumption statistic.
enum {
   mMesh,
   mGhost,
   mParticles
};

int           memEstimateOnly;			///< Flag singnals if we are in \b estimate mode (no allocations, only RAM usage counts and tag parameters checks).
static double memLayout[3] = {0, 0, 0};		///< Estimated amount of memory consumed (OS overhead might be bigger by 1 - 10%).

static meshDouble_t rho;			///< Global meshes.
static meshVec_t    E, H, J;

// ---------------------------------------------------------------------------
/// Sets fundamental parameters of the domain: sizes, BC and partitioning.
/// XXX Remove 2 nodes wide ghost cells (only one node is necessary).
// ---------------------------------------------------------------------------
static void
main_setupDomain (FILE *fp)
{
   int    tags[3]          = {TAG_UNITS, TAG_MESH, TAG_BOUNDARY};
   void (*func[3]) (FILE*) = {tag_units, tag_mesh, tag_boundary};

   // Forces first three tags to be 'units', 'mesh', 'boundary', and processes
   // them in this order.
   for (int i = 0 ; i < 3 ; ++i) {
      ENSURE (!feof (fp) && getTag (fp) == tags[i],
              "first tags must be [units]/[mesh]/[boundary].");
      MPI_Barrier (MPI_COMM_WORLD);
      func[i] (fp);
   }

   parameter_dump ();
   partition_init ();
   partition_show ();

   // Allocates meshes.
   if (!memEstimateOnly) {
      mesh_allocate (mcast_mesh (&E),   cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2, cpu_max[0] + 2, cpu_max[1] + 2, cpu_max[2] + 2, "setup_E",   mc_vec3D_t);
      mesh_allocate (mcast_mesh (&H),   cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2, cpu_max[0] + 2, cpu_max[1] + 2, cpu_max[2] + 2, "setup_H",   mc_vec3D_t);
      mesh_allocate (mcast_mesh (&J),   cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2, cpu_max[0] + 2, cpu_max[1] + 2, cpu_max[2] + 2, "setup_J",   mc_vec3D_t);
      mesh_allocate (mcast_mesh (&rho), cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2, cpu_max[0] + 2, cpu_max[1] + 2, cpu_max[2] + 2, "setup_rho", mc_double);

      mf_mesh_clean (&E);
      mf_mesh_clean (&H);
      mf_mesh_clean (&J);
      mf_mesh_clean (&rho);
   }

   // Bodies of the meshes (rho accounted).
   memLayout[mMesh] = (sizeof (double)*2.0 + sizeof (vec3D_t)*3)
                      *(mc_have_x*(cpu_max[0] - cpu_min[0] + 4) + 1)
                      *(mc_have_y*(cpu_max[1] - cpu_min[1] + 4) + 1)
                      *(mc_have_z*(cpu_max[2] - cpu_min[2] + 4) + 1);
}

// ---------------------------------------------------------------------------
/// Converts string into '2.34e5 (1G 10M 123K)' form (one can use up to 8 last strings).
// ---------------------------------------------------------------------------
static char*
main_memString (double size, char *buff)
{
   *buff = 0;
   double radix = (1 << 10),
	  factor = radix*radix*radix*radix;			// Starts from terabytes.
   char symbol[5] = {'b', 'K', 'M', 'G', 'T'}, *slot = buff;			// Size prefixes.
   for (int digit = 4 ; digit >= 0 ; --digit, factor /= radix) {		// Computes all digits.
      if ((int) fmod (size/factor, radix) > 0)
         sprintf (slot, "%d%c ", (int) fmod (size/factor, radix), symbol[digit]);
      slot += strlen (slot);
   }
   sprintf (slot, "(%.3e)", size);						// Adds size in exponential form.
   return buff;
}

// ---------------------------------------------------------------------------
/// Reports expected memory usage and exits if we are in \c "estimate" mode.
// ---------------------------------------------------------------------------
static void
main_memReport (void)
{
   const int sx = mc_have_x*(cpu_max[0] - cpu_min[0] + 4) + 1,	// Calculates sizes of domain (in nodes).
             sy = mc_have_y*(cpu_max[1] - cpu_min[1] + 4) + 1,
             sz = mc_have_z*(cpu_max[2] - cpu_min[2] + 4) + 1;

   // Ghost cells for rho, for E/H (2 layers), for J (4 layers), 2 sides for each direction.
   double totalMem[3] = {0, 0, 0};
   memLayout[mGhost] = (2.0*sizeof (double) + 2.0*2*sizeof (vec3D_t) +  4*sizeof (vec3D_t))*
                        2*(mc_have_x*sy*sz + mc_have_y*sx*sz + mc_have_z*sx*sy);
   MPI_Allreduce (memLayout, totalMem, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	// Calculates total memory usage.

   char b1[200], b2[200];
   say ("\nMemory usage (node/total):");
   say ("  o mesh: %s \t/ %s", main_memString (memLayout[mMesh], b1), main_memString (totalMem[mMesh], b2));
   say ("  o particles: %s \t/ %s", main_memString (memLayout[mParticles], b1), main_memString (totalMem[mParticles], b2));
   say ("  o ghost: %s \t/ %s", main_memString (memLayout[mGhost], b1), main_memString (totalMem[mGhost], b2));
   say ("  o all: %s \t/ %s", main_memString (memLayout[0] + memLayout[1] + memLayout[2], b1),
                              main_memString (totalMem[0] + totalMem[1] + totalMem[2], b2));
}

// ---------------------------------------------------------------------------
/// Clears file left from old runs (done by root cpu only, uses bin/sh - Windows needs reimplementation).
// ---------------------------------------------------------------------------
static void
main_clearFolders (void)
{
   if (cpu_here || memEstimateOnly)
      return;

   // Removes all intermediate config files (setup should generate new ones).
   ENSURE (system (NULL), "shell is not avaliable.");
   system ("rm -f .gaussSpot.cfg .EM_sources/TFSF/TFSF.cfg "
                 ".EM_sources/TFSF/openFaces.cfg ERROR.*");
   system ("for f in ./binData/*.* ; do rm -f $f ; done");
}

// ---------------------------------------------------------------------------
/// Entry point for \b setup.out module.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
   ENSURE (argc >= 3,
           "Usage: %s <configuration file name> {estimate|create}", argv[0]);

   cl_import ("./source/import.pl", argc, argv, 2);
   parameter_enterMPI (argc, argv, 0);

   memEstimateOnly = cfg_identifyWord (argv[2], "estimate", 1,
                                                "create",   0,
                                                mc_cfgTermGuesses);
   ENSURE (memEstimateOnly >= 0,
           "usage: %s <configuration file name> <estimate|create>",
           argv[0]);

   main_clearFolders ();

   say ("Config file: %s", argv[1]);

   FILE *fp = cfg_open (argv[1], "rt", "setup.out");

   // Reads first three key tags (domain, boundary conditions, units).
   main_setupDomain (fp);

   while (!feof (fp)) {
      int tag = getTag (fp);
      MPI_Barrier (MPI_COMM_WORLD);
      switch (tag) {
         case TAG_PLASMA:		memLayout[mParticles] += tag_plasma (fp);		break;
         case TAG_DF_UNIFORM:		tag_DF_uniform (fp);					break;
         case TAG_TWO_STREAM:		memLayout[mParticles] += tag_twoStream (fp);		break;
         case TAG_FOIL:			memLayout[mParticles] += tag_foil (fp);			break;
         case TAG_PHOTOELECTRONS:	memLayout[mParticles] += tag_photoelectrons (fp);	break;
         case TAG_MAXWELL:		memLayout[mParticles] += tag_maxwell (fp);		break;
         case TAG_CLUSTER:		memLayout[mParticles] += tag_cluster (fp);		break;
         case TAG_TRIANGLE_PRIZM:	memLayout[mParticles] += tag_trianglePrizm (fp);	break;

         case TAG_TFSF:			tag_TFSF (fp);						break;
         case TAG_TFSF_FACES:		tag_TFSFOpenFaces (fp);					break;
         case TAG_SRC_MIRROR:		memLayout[mMesh] += tag_laser (fp);			break;

         case TAG_EM_POINT:		tag_point (fp, &E, &H);					break;
         case TAG_EM_WAVE:		tag_EMWave (fp, &E, &H);				break;
         case TAG_EM_GAUSS_SPOT:	tag_gaussSpot (fp, argv[1]);				break;
         case TAG_EM_PLANE_PULSE:	tag_planePulse (fp, &E, &H);				break;
         case TAG_EM_RESONATOR:		tag_EMResonator (fp, &E);				break;

         case TAG_SEED_PITS:		tag_seed_PITS (fp, &E);					break;
         case TAG_SEED_WEIBEL:		tag_seedWeibel (fp, &H);				break;
         case TAG_RING_DF:		tag_ringDF (fp);					break;
         case TAG_PLASMA_WAVE:		tag_plasmaWave (fp);					break;
         case TAG_VELOCITY_SHIFT:	tag_velocityShift (fp);					break;
         case TAG_CURRENT_BALANCER:	tag_currentBalancer (fp);				break;
         case TAG_MEAN_V_NOISE:		tag_meanVNoise (fp);					break;
         case TAG_SCALES:		tag_scales (fp);					break;
         case TAG_SCISSORS:		tag_scissors (fp);					break;
         case TAG_GRADIENT:		tag_gradient (fp);					break;

      case TAG_EOF:		// End of file is reached, no tags found.
	 break;

      case TAG_UNKNOWN_TAG:	// Unknown tag - terminating, setup is useless anyway now.
	 DIE ("unknown tag '%s' (bad syntax?)", getLastTagName ());

      default:
	 DIE ("unregistered switch value '%d' (%s), update sources.",
	       tag, getLastTagName ());
      }
   }
   fclose (fp);

   main_memReport ();								// Reports the memory consuming.

   if (!memEstimateOnly) {							// Reports memory usage and quits.
      say ("Evaluating charge density...");					// Charge density (for visualization).
      plasma_rho (&rho);

      say ("Writing tecplot data-point...");
      if (!cpu_here)								// Master saves new units.
         units_save ();

      fileMap_init (mc_fileMap_setup);
      tecIO_setRecordNum (0);							// Creates start tecplot visualization record.
      tecIO_saveFields (Time, mcast_meshVec_RO (&E), mcast_meshVec_RO (&H));
      tecIO_saveCurrents (mcast_meshVec_RO (&J), mcast_meshDouble_RO (&rho));

      say ("Writing start check-point...");
      sysIO_setRecordNum (0);							// Creates start checkpoint.
      sysIO_save (Time, mcast_meshVec_RO (&E), mcast_meshVec_RO (&H));
      plasma_save (IO_plasmaName (0, cpu_here));				// Saves plasma into checkpoint #0.
#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
      FILE *fp = cfg_open ("binData/tracer.inf", "wt", __func__);
      fprintf (fp, "%le \tτ (time step)\n",        tau);
      fprintf (fp, "%d  \tnumber of trajectories", 0);
      fclose (fp);
#endif
   }

   say ("All done.\nBye :)");

   return EXIT_SUCCESS;
}
