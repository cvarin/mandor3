/** \file main.c
  * Converter of the dumped EM field into spectral energy distribution.
  */

#define mc_convertTimeLimit		(55*60)		///< Time limit on convertion (55 minutes).

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "frame.h"
#include "type_reg.h"
#include "type_CFile.h"

#include "log.h"
#include "timer.h"
#include "misc_units.h"
#include "misc_cfgReader.h"

#include "spectr_process.h"

static int recStart, recEnd, recStep;			///< First, last and interval between records to process.
static int unitsMicronFemtosec;				///< Type of the units of the output data.
static int NKx, NKy, NKz;				///< Spectrum domain limits: [-NKx, NKx] x .. x [-NKz, NKz].

static reg_t subdomain = {				///< Subdomain where to take data for fft (covers everything).
  {INT_MIN, INT_MIN, INT_MIN},
  {INT_MAX, INT_MAX, INT_MAX}
};

static double A0;					///< Unit of E.
static double femtosecond;				///< Unit of time.

static timeTick_t programStart;				///< Time to count execution time from.
static int removeDumps = 0;				///< Flag to clear dumped files.

// ---------------------------------------------------------------------------
/// Reads parameters of the config file.
// ---------------------------------------------------------------------------
static void
main_readConfig (void)
{
  FILE *fp = cfg_open ("diag_spectr.cfg", "rt", __func__);								// Reads parameters of the record set to draw.
  recStart = cfg_readInt (fp);
  recEnd = cfg_readInt (fp);
  recStep = cfg_readInt (fp);
  unitsMicronFemtosec = cfg_readInt (fp) == 1;										// Clamps to use as array index.
  NKx = cfg_readInt (fp);
  NKy = cfg_readInt (fp);
  NKz = cfg_readInt (fp);												// Reads region of the spectral domain to fill.
  if (cfg_isOption (fp))
  {
    const char *optionName = cfg_readOptWord (fp);
    int optionType = cfg_identifyWord (optionName, "subdomain", 0, mc_cfgTermGuesses);

    if (optionType == -1)												// Checks if option is unknown.
      error ("main_readConfig: unknown option '%s'", optionName);

    subdomain.min[0] = cfg_readOptInt (fp);
    subdomain.min[1] = cfg_readOptInt (fp);
    subdomain.min[2] = cfg_readOptInt (fp);
    subdomain.max[0] = cfg_readOptInt (fp);
    subdomain.max[1] = cfg_readOptInt (fp);
    subdomain.max[2] = cfg_readOptInt (fp);
  }
  mf_reg_collapse (&subdomain);												// Accounts for 1D/2D cases.
  fclose (fp);

  units_load ();													// Loads all units.
  A0 = units (mc_A0)/units (mc_E0);
  femtosecond = (unitsMicronFemtosec) ? units (mc_femtosecond) : 1;

  fp = cfg_open ("binData/spectr.N", "rt", __func__);									// Updates total number of records.
  int lastDumpNum = cfg_readInt (fp) - 1;										// Updates recEnd.
  fclose (fp);

  recEnd = (recEnd < 0 || lastDumpNum < recEnd) ? lastDumpNum : recEnd;

  if (recStart < 0 || recStep <= 0 || recStart > recEnd)								// Checks set of the records.
    error ("diag_tecplot.cfg: bad start (%d)/step(%d)/end(%d) combination.", recStart, recStep, recEnd);
}

// ---------------------------------------------------------------------------
/// Converts all dumps into spectral energy density distributions.
// ---------------------------------------------------------------------------
static void
main_convertDumps (double timeLimit)
{
  say_doing ("Start to convert dump files ...");

  double estimate = 0;
  for (int rec = recStart ; rec <= recEnd ; rec += recStep)								// Converts dumps to spectral density data.
  {
    timeTick_t start, end;
    time_get (&start);
    spectr_convert (rec, NKx, NKy, NKz, &subdomain, removeDumps);
    time_get (&end);

    double convertTime = time_elapsed (&start, &end);
    estimate = (0.3*estimate + 0.7*convertTime);									// Weighted iteration time.
    if (time_elapsed (&programStart, &end) + 2*estimate > timeLimit)
      error ("main_convert: out of time for FFT convertion - restart diagnostic to continue.");
    say_doing ("Record %d is converted in %.3e sec ...", rec, convertTime);
  }

  say_doing ("Dump files are converted.");
}

// ---------------------------------------------------------------------------
/// Opens prepared by spectr_processDump() file and gets parameters (spectral domain size, etc).
// ---------------------------------------------------------------------------
static void
main_readGlobals (CFile_t *file, int *NKx, int *NKy, int *NKz, double *time, double *Lx, double *Ly, double *Lz)
{
  int found = (CF_findChunk (file, "spectral domain size") == 0);
  assert (found);

  CF_read (NKx, sizeof (int), 1, file);
  CF_read (NKy, sizeof (int), 1, file);
  CF_read (NKz, sizeof (int), 1, file);
  CF_read (time, sizeof (double), 1, file);
  CF_read (Lx, sizeof (double), 1, file);
  CF_read (Ly, sizeof (double), 1, file);
  CF_read (Lz, sizeof (double), 1, file);
}

// ---------------------------------------------------------------------------
/// Opens prepared by spectr_processDump() file and prints it in tecplot block format (\b warning: energy is rescaled to units [A0^2]).
// ---------------------------------------------------------------------------
static void
main_exportBlock (CFile_t *file, const char *chunk, FILE *fp)
{
  int found = (CF_findChunk (file, chunk) == 0);
  assert (found);

  const float scale = 1.0/(A0*A0);
  for (int i = - NKx ; i <= NKx ; i++)
    for (int j = - NKy ; j <= NKy ; j++)
      for (int k = - NKz ; k <= NKz ; k++)
      {
        float W;
        CF_read (&W, sizeof (float), 1, file);
        fprintf (fp, "%.4e\n", W*scale);
      }
}

// ---------------------------------------------------------------------------
/// Entry point for spectr.out diagnostic.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
  time_get (&programStart);												// Gets start time.

  msg_openLogFile ("output/spectr.log", 0);
  msg_setRefreshTime (0.2);												// Updates screen logger's frequency.

  for (int i = 1 ; i < argc ; ++i)											// Checks if disk space clearing is requested.
    if (!strcmp (argv[i], "--remove-dumps") || !strcmp (argv[i], "-r"))
    {
      removeDumps = 1;
      say ("%s: diagnostic will remove all dump files.", argv[0]);
    }

  main_readConfig ();													// Loads configuration.

  timeTick_t now;
  time_get (&now);
  main_convertDumps (mc_convertTimeLimit - time_elapsed (&programStart, &now));						// Converts dumps to spectral density data.

  FILE *fp = cfg_open ("output/spectr.dat", "wt", __func__);								// Opens output file.
  fprintf (fp, "variables = mz, my, mx, \"E<sub>x k</sub><sup>2</sup> [a<sub>0</sub><sup>2</sup>]\", "			// Order is reversed due to tecplot order.
                                       "\"E<sub>y k</sub><sup>2</sup> [a<sub>0</sub><sup>2</sup>]\", "
                                       "\"E<sub>z k</sub><sup>2</sup> [a<sub>0</sub><sup>2</sup>]\", "
                                       "\"H<sub>x k</sub><sup>2</sup> [a<sub>0</sub><sup>2</sup>]\", "
                                       "\"H<sub>y k</sub><sup>2</sup> [a<sub>0</sub><sup>2</sup>]\", "
                                       "\"H<sub>z k</sub><sup>2</sup> [a<sub>0</sub><sup>2</sup>]\"\n");

  double estimate = 0;
  regList_t zones = mc_regList_init;											// List of meshes (to share, cpu = zone).
  for (int rec = recStart, zone = 1 ; rec <= recEnd ; rec += recStep, ++zone)						// Exports spectral density data to tecplot.
  {
    double time, Lx, Ly, Lz;

    timeTick_t start;													// Remembers start.
    time_get (&start);

    CFile_t *file = CF_openRead ("binData/spectr_%d", rec);
    assert (file);
    main_readGlobals (file, &NKx, &NKy, &NKz, &time, &Lx, &Ly, &Lz);							// Gets parameters of the converted dump.
    fprintf (fp, "zone t=\"t[%s]=%.3f\", i = %d, j = %d, k = %d, f = block\n",  					// Saves title.
     	     (unitsMicronFemtosec) ? "fs" : "t<sub>0</sub>", time/femtosecond, 2*NKz + 1, 2*NKy + 1, 2*NKx + 1);

    int l = -1;														// Looks if mesh with this sizes is written.
    for (reg_t *mesh = zones.list, *end = mesh + zones.N ; mesh < end ; ++mesh)
      if ( ((mesh->min[0] == - NKx && mesh->max[0] == NKx) || (!mc_have_x)) &&
           ((mesh->min[1] == - NKy && mesh->max[1] == NKy) || (!mc_have_y)) &&
           ((mesh->min[2] == - NKz && mesh->max[2] == NKz) || (!mc_have_z)) )
      {
        l = mesh->barcode;
        break;
      }

    if (l == -1)												// Saves zone blocks.
    {
      for (int i = -NKx ; i <= NKx ; ++i)										// Saves mesh node k-coordinates.
        for (int j = -NKy ; j <= NKy ; ++j)
          for (int k = -NKz ; k <= NKz ; ++k)
            fprintf (fp, "%d\n", k);

      for (int i = -NKx ; i <= NKx ; ++i)										// Saves mesh node j-coordinates.
        for (int j = -NKy ; j <= NKy ; ++j)
          for (int k = -NKz ; k <= NKz ; ++k)
            fprintf (fp, "%d\n", j);

      for (int i = -NKx ; i <= NKx ; ++i)										// Saves mesh node i-coordinates.
        for (int j = -NKy ; j <= NKy ; ++j)
          for (int k = -NKz ; k <= NKz ; ++k)
            fprintf (fp, "%d\n", i);

      reg_t reg = {{-NKx, -NKy, -NKz}, {+NKx, +NKy, +NKz}};
      reg.barcode = zone;
      regList_add (&zones, &reg);
    }
    else
    {
      fprintf (fp, "varsharelist = ([1-3]=%d)\n", l);									// Reuses mesh node coordinates.
    }

    main_exportBlock (file, "E.x", fp);											// Exports all scaled energies.
    main_exportBlock (file, "E.y", fp);
    main_exportBlock (file, "E.z", fp);
    main_exportBlock (file, "H.x", fp);
    main_exportBlock (file, "H.y", fp);
    main_exportBlock (file, "H.z", fp);
    CF_close (file);

    say ("Record %d (t = %.3e [t0]): domain size [%e, %e, %e].", rec, time, Lx*mc_have_x, Ly*mc_have_y, Lz*mc_have_z);

    time_get (&now);
    time = time_elapsed (&start, &now);											// Gets iteration time.
    estimate = (0.3*estimate + 0.7*time);										// Averages iteration time.
    if (time_elapsed (&programStart, &now) + 1.2*estimate > mc_convertTimeLimit)					// Checks limit.
    {
      fclose (fp);
      error ("spectr.out: out of time for output (consider using PBS for long job or restart if some time was absorbed by fft).");
    }
  }
  fclose (fp);

  say ("Done.");

  return 0;
}
