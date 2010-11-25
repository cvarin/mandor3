/*
 * Phasedot. Simple converter of the binary data file with DF into plt file with \vec r, \gamma\vec v.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "type_marker.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "IO_sys.h"

#define XAxis 0				///< Predefined const for X axis.
#define YAxis 1				///< Predefined const for Y axis.
#define ZAxis 2				///< Predefined const for Z axis.

static int  markerN = 0;		///< Number of markers passed to the output.

static int     filter_types = 0;
static double *filter_qDivM = NULL;

// ---------------------------------------------------------------------------
/// Prints all possible charge to mass ratios.
// ---------------------------------------------------------------------------
static void
filter_prepare (void)
{
  filter_types = markers_dumpQDivM (&filter_qDivM);

  for (int p = 0 ; p < filter_types ; ++p)
    say ("Component %d has q/M = %.3e", p, filter_qDivM[p]);
}

// ---------------------------------------------------------------------------
/// Filters all particles loaded and keeps only one component we want.
// ---------------------------------------------------------------------------
static void
filter_keep (double qDivM)
{
  filter_types = markers_dumpQDivM (&filter_qDivM);

  if (fabs (qDivM) > 1e8)												// Switch to keep all particles.
    return;

  int p = 0;
  for ( ; p < filter_types ; ++p)											// Gets type-ID of the particles.
    if (fabs (qDivM - filter_qDivM[p]) < 1e-3*fabs (filter_qDivM[p]))
      break;

  if (p >= filter_types)
  {
    filter_prepare ();
    error ("filter_keep: cannot define type-ID of the component with q/M = %e (all components are listed above).", qDivM);
  }

  markerN = 0;														// Removes unnecessary pages.
  for (int ID = markerChapter_first () ; ID >= 0 ; markerChapter_next (&ID))
  {
    #pragma set woff 1343												// To avoid warning on SGI.
    markerIterator_t page;
    #pragma reset woff 1343
    for (markerPage_first (ID, &page) ; page.df ; markerPage_next (&page))
      if (p != page.type)
        markerPage_setN (&page, 0);
  }
}

// ---------------------------------------------------------------------------
/// Entry point for \b phasedot diagnostic.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
  parameter_enterMPI (argc, argv, 0, NULL);										// Initializes log-file and MPI frame.

  const int lastCheckPoint = sysIO_recordsTotal () - 1;									// Loads last record number.
  units_load ();													// Loads all units.

  FILE *fp = cfg_open ("diag_phasedot.cfg", "rt", "tecplot.out");							// Reads parameters of the data processing.
  int start = cfg_readInt (fp);
  int end = cfg_readInt (fp);
  int step = cfg_readInt (fp);
  int writePLTFile = cfg_readInt (fp);											///< \todo Banish this dummy to hell.
  int unitsMicronFemtosec = cfg_readInt (fp) == 1;									// Clamps to use as array index.
  double qDivM = cfg_readDouble (fp);											// Gets desired q/M ratio.
  fclose (fp);

  if (end < 0 || end > lastCheckPoint)											// Sets range of records (max or given).
    end = lastCheckPoint;

  if (start < 0 || start > end)
    error ("phasedot.cfg: bad start record number.");

  if (step <= 0)
    error ("phasedot.cfg: bad step.");

  double femtosecond = 1, micron = 1;											// Chooses units.
  if (unitsMicronFemtosec)
  {
    femtosecond = units (mc_femtosecond);
    micron = units (mc_micron);
  }

  msg_setRefreshTime (0.2);
  say ("Start to generate phase-plane(dot) data file...\n");

  fp = cfg_open ("output/phasedot.dat", "wt", "phasedot.out");
  fprintf (fp, "variables = x, y, z, `gv_x, `gv_y, `gv_z, q\n");
  for (int record = start ; record <= end ; record += step)
  {
    double time;
    int ID, revision;
    markerIterator_t page;

    sysIO_parameters (record, &revision, &time, NULL, NULL, NULL);							// Loads DF.
    parameter_load (revision);
    partition_init ();

    reg_t reg;
    reg.min[0] = dmn.imin - 10;
    reg.min[1] = dmn.jmin - 10;
    reg.min[2] = dmn.kmin - 10;
    reg.max[0] = dmn.imax + 10;
    reg.max[1] = dmn.jmax + 10;
    reg.max[2] = dmn.kmax + 10;
    sysIO_loadPlasma (record, &time, &reg);

//     filter_prepare ();
    filter_keep (qDivM);

    markerN = 0;													// Finds the total number of particles.
    for (ID = markerChapter_first () ; ID >= 0 ; markerChapter_next (&ID))
      for (markerPage_first (ID, &page) ; page.df ; markerPage_next (&page))
        markerN += page.N;

    static const char titles[2][17] = {"t[t<sub>0</sub>]\0", "t[fs]\0"};
    fprintf (fp, "zone t=\"%s=%.3f\", i = %d, f = point\n", titles[unitsMicronFemtosec], time/femtosecond, markerN);

    for (ID = markerChapter_first () ; ID >= 0 ; markerChapter_next (&ID))
      for (markerPage_first (ID, &page) ; page.df ; markerPage_next (&page))
      {
        marker_t *f = page.df;
        for (int p = 0 ; p < page.N ; p++)
          fprintf (fp, "%.4e %.4e %.4e %.4e %.4e %.4e %.4e\n", f[p].x/micron, f[p].y/micron, f[p].z/micron, f[p].vx, f[p].vy, f[p].vz, f[p].rho);
      }

    say_doing ("time = %.4f (%d particles)", time, markerN);
  }

  fclose (fp);

  return EXIT_SUCCESS;
}
