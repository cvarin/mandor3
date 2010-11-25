/** \file tag_TFSF.c
  * Function to read \b [TFSF] tag parameters from config file and prepare
  * hidden config file for \b core.out (main calculation module).
  *
  * \sa   tag_TFSF.h, tag_laser.h.
  * \todo Link on main documentation page for \b em_TFSF.c.
  */

#include <math.h>
#include <dirent.h>		// For folder scanning.
#include <unistd.h>		// For chdir.
#include <limits.h>		// For chdir.

#include "type_mesh.h"

#include "core/em_TFSF.h"	///< To get all IDs of interface modes.

#include "main.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"

enum {BLEND_NEW, BLEND_ADD};

static int first = 1;

// ---------------------------------------------------------------------------
/// Removes all files in folder '.EM_sources/TFSF'.
// ---------------------------------------------------------------------------
static void
tag_removeAllTFSFFiles (void)
{
   ENSURE (!chdir (".EM_sources/TFSF"), "cannot chdir into the TFSF folder");

   DIR           *folder = opendir (".");
   struct dirent *item;
   while ((item = readdir (folder)) != NULL) {
      remove (item->d_name);
   }
   closedir (folder);
   chdir ("../..");

   SAY_DEBUG ("All files in folder '.EM_sources/TFSF' are removed.");
}

// ---------------------------------------------------------------------------
/// Checks if previously written data are consistent with current parameters.
// ---------------------------------------------------------------------------
static int
reg_is_the_same (const reg_t *reg)
{
    FILE *fp = fopen (".EM_sources/TFSF/consistencyCheck.dat", "rt");
    if (!fp)
        return 1;

    int    di = cfg_readInt (fp),
           dj = cfg_readInt (fp),
           dk = cfg_readInt (fp);
    double dt = cfg_readDouble (fp),
           dx = cfg_readDouble (fp),
           dy = cfg_readDouble (fp),
           dz = cfg_readDouble (fp);

    return    reg->max[0] - reg->min[0] == di
           && reg->max[1] - reg->min[1] == dj
           && reg->max[2] - reg->min[2] == dk
           &&           fabs (dt - tau) < 1e-13*tau
           && mc_have_x*fabs (dx - h1)  < 1e-13*h1
           && mc_have_y*fabs (dy - h2)  < 1e-13*h2
           && mc_have_z*fabs (dz - h3)  < 1e-13*h3;
}

// ---------------------------------------------------------------------------
/// TFSF setup entry point.
// ---------------------------------------------------------------------------
void
tag_TFSF (FILE *fp)
{
    int   mode = -1;
    reg_t reg  = mc_reg_init;

    ENSURE (first, "second [TFSF] tag  is not allowed");
    first = 0;

    reg.min[0] = cfg_readInt (fp); // Reads interfaced region bounding box.
    reg.max[0] = cfg_readInt (fp);
    reg.min[1] = cfg_readInt (fp);
    reg.max[1] = cfg_readInt (fp);
    reg.min[2] = cfg_readInt (fp);
    reg.max[2] = cfg_readInt (fp);
    mf_reg_collapse (&reg);

    const char *mode_str = cfg_readWord (fp);

    // Gets mode ID to use in the code.
    mode = cfg_identifyWord (mode_str,
                             "record",       mc_TFSF_rec,
                             "playForward",  mc_TFSF_playForward,
                             "playBackward", mc_TFSF_playBackward,
                             mc_cfgTermGuesses);
    ENSURE (mode != -1, "unknown mode '%s'", mode_str);

    if ( (mc_have_x && (reg.min[0] <= dmn_min[0] + 5 || reg.max[0] >= dmn_max[0] - 5)) ||
         (mc_have_y && (reg.min[1] <= dmn_min[1] + 5 || reg.max[1] >= dmn_max[1] - 5)) ||
         (mc_have_z && (reg.min[2] <= dmn_min[2] + 5 || reg.max[2] >= dmn_max[2] - 5)) ) {
        reg_t domain = {
            {dmn_min[0], dmn_min[1], dmn_min[2]},
            {dmn_max[0], dmn_max[1], dmn_max[2]}
        };
        DIE ("interfaced region '%s' is too big for domain '%s'",
             reg_printRanges (&reg), reg_printRanges (&domain));
    }

    int   blend_mode   = BLEND_NEW,
          filter_min   = INT_MIN,
          filter_max   = INT_MAX,
          shift_frames = 0;;
    if (cfg_isOption (fp)) {
        const char *mode = cfg_readOptWord (fp);
        if (!strcmp (mode, "add")) {
            blend_mode = BLEND_ADD;
            ENSURE (reg_is_the_same (&reg),
                       "new TFSF region differs from already written one");
        }
        filter_min = cfg_readOptInt (fp);
        filter_max = cfg_readOptInt (fp);
        if (cfg_isOption (fp))
            shift_frames = cfg_readOptInt (fp);
    }

    say ("tag_TFSF:\n  - config file 'TFSF.cfg' is created.");
    say ("  - approximate crossing time of bbox in X direction is %d time "
              "steps,", (int) ((reg.max[0] - reg.min[0] + 1)*h1*units (mc_r0)/
                               (mc_CGS_c*tau*units (mc_t0))));
    say ("  - interfaced region in mesh steps/dimensionless units/[μm]:");
    say ("      %s [nodes],", reg_printRanges (&reg));
    say ("      [%.3f, %.3f]x[%.3f, %.3f]x[%.3f, %.3f] [r0],",
                reg.min[0]*h1, reg.max[0]*h1,
                reg.min[1]*h2, reg.max[1]*h2,
                reg.min[2]*h3, reg.max[2]*h3);
    say ("      [%.3f, %.3f]x[%.3f, %.3f]x[%.3f, %.3f] [μm].",
              reg.min[0]*h1/units (mc_micron), reg.max[0]*h1/units (mc_micron),
              reg.min[1]*h2/units (mc_micron), reg.max[1]*h2/units (mc_micron),
              reg.min[2]*h3/units (mc_micron), reg.max[2]*h3/units
(mc_micron));
    say ("  - interface mode is '%s'.", mode_str);
    if (blend_mode == BLEND_ADD) {
        say ("  - signal will be added the already recorded.");
        say ("  - sensitive part of the interface: i∊[%d %d].",
                  filter_min, filter_max);
    }
    say ("  Please double-check output of the 'core.out' with REAL "
              "parameters of the TFSF interface in use.");

    if (memEstimateOnly)
        return;

    // Master node updates all files.
    if (!cpu_here) {
        if (mode == mc_TFSF_rec && blend_mode == BLEND_NEW)
            tag_removeAllTFSFFiles ();

        FILE *fp = cfg_open (".EM_sources/TFSF/TFSF.cfg", "wt", __func__);
        fprintf (fp, "@ %5d       TFSF interface bounding box: reg.min[0].\n", reg.min[0]*mc_have_x);
        fprintf (fp, "@ %5d       TFSF interface bounding box: reg.max[0].\n", reg.max[0]*mc_have_x);
        fprintf (fp, "@ %5d       TFSF interface bounding box: reg.min[1].\n", reg.min[1]*mc_have_y);
        fprintf (fp, "@ %5d       TFSF interface bounding box: reg.max[1].\n", reg.max[1]*mc_have_y);
        fprintf (fp, "@ %5d       TFSF interface bounding box: reg.min[2].\n", reg.min[2]*mc_have_z);
        fprintf (fp, "@ %5d       TFSF interface bounding box: reg.max[2].\n", reg.max[2]*mc_have_z);
        fprintf (fp, "@ %5d       Interface regime (%s).\n", mode, mode_str);
        fprintf (fp, "\n\n   This file is automatically created by setup.out "
                     "(see output/setup_*.log files).\n");
        fclose (fp);

        fp = cfg_open (".EM_sources/TFSF/filter.cfg", "wt", __func__);
        fprintf (fp, "@ % 5d       I0.\n", filter_min);
        fprintf (fp, "@ % 5d       I1.\n", filter_max);
        fprintf (fp, "@ % 5d       Shift frames.\n", shift_frames);
        fclose (fp);
    }
}
