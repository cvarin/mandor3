/** \file tag_laser.c
  * Interface to initialize soft sources of specific form.
  *
  * \sa tag_laser.h, em_mirror.h, em_mirror.c, em_sources.c.
  */

#include <math.h>
#include <dirent.h>		// For folder scanning.
#include <unistd.h>		// For chdir.
#include <limits.h>		// For chdir.

#include "type_mesh.h"

#include "main.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"

#include "scpic/em_sources.h"

static int first = 1;

// ---------------------------------------------------------------------------
/// Removes all files in folder '.EM_sources/laser'.
// ---------------------------------------------------------------------------
static void
remove_old_data (void)
{
   DIR    *folder = opendir (".EM_sources/laser");
   ENSURE (folder, "cannot access focused laser data folder");

   struct dirent *file;
   while ((file = readdir (folder))) {
      remove (_(".EM_sources/laser/%s", file->d_name));
   }
   closedir (folder);

   SAY_DEBUG ("Removed all files in folder '.EM_sources/laser'.");
}

// ---------------------------------------------------------------------------
/// Extracts parameters and caches precomputed soft source of EM waves.
// ---------------------------------------------------------------------------
double
tag_laser (FILE *fp)
{
   if (first && !memEstimateOnly)
      remove_old_data ();
   first = 0;

   ENSURE (mc_have_x && mc_have_y && mc_have_z,
           "Stratton-Chu integrals are written for 3D case only");

   mirror_t src = read_mirror_parameters (fp);
   print_mirror_parameters (src);

   if (!memEstimateOnly) {
      add_new_mirror (src);
   }

   return (dmn_max[1] + 1)*(dmn_max[2] + 1)*2.0*sizeof (wave_x_t);
}
