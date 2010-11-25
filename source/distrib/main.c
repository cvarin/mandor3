/** \file main.c
  * \brief Diagnostic to build distribution function (histogram) of particles.
  *
  * \todo Make it 1D/2D/3D/.. using general way somehow.
  */

/** \mainpage Distribution function diagnostic module of the "Mandor".
  * This diagnostic is used to specify the set of particles to analyse and than to
  * build the distribution function of the particles under interest. Distribution
  * function is a 2D histogram with number of bins and region to draw passed through
  * config-file \b diag_distrib.cfg. See main.c for details.
  */

#include <math.h>
#include <stdlib.h>

#include <mpi.h>

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

#include "IO_sys.h"
#include "IO_names.h"

#define mc_nameLength		64						///< Variable/unit name length.
#define mc_pqrN			11						///< Number of possible coordinates.
#define mc_chunkSize 		50000						///< Number of markers to cache from disk for postprocessing.

#define mc_fastEnergyTableStepInversed	(1000.0)
#define mc_fastEnergyTableLimit	        (1.0)
static double *fastEnergyTable = NULL;						///< Look-up table for \f$ (\gamma - 1) \f$ calculation.

static marker_t chunk[mc_chunkSize];						///< Markers cache.

/// Function to extract value from the marker (coordinate, energy, velocity).
typedef double (*param_f) (const marker_t *p);

/// Structure to hold all parameters of the general coordinate (unit, id, name for Tecplot, function-generator, etc).
typedef struct
{
  const char name[mc_nameLength];						///< Name of the variable.
  const      param_f value;							///< Function to generate value for distribution.
  double     unitValue;								///< Unit for value.
  char       unitName[mc_nameLength];						///< Unit name for value.
} pqrVar_t;

static int    planeActivated = 0;						///< Flag to activate plane slicing.
static double planeA = 0, planeB = 0, planeC = 0, planeD = 1;			///< Default plane - all particles are over it.

/// Time (dimensionless) and unit of time.
static double time, unitTime;
static char unitTimeName[mc_nameLength];					///< Name of the time units ("fs"/"t_0"/<whatever>).

/// Parameters of the record set.
static int diagStart, diagEnd, diagStep;

/// DF agruments / number of bins.
static int pVar, qVar, rVar, pN, qN, rN;

/// DF region and histogram's body.
static double pMin, pMax, qMin, qMax, rMin, rMax;
static float  *diagHist = NULL;

static int    markerN = 100111412;						///< Total number of the markers sended to histogram accumulator.
static double filter_qDivM = 1e21;						///< Filter to grab desired \b q/M.

static double gamma2eV = 0;							///< Current factor of energy conversion \f$ (\gamma - 1) \to \varepsilon\, [eV]\f$.

static double param_x (const marker_t *p) { return p->x; }			///< Returns x.
static double param_y (const marker_t *p) { return p->y; }			///< Returns y.
static double param_z (const marker_t *p) { return p->z; }			///< Returns z.
static double param_vx (const marker_t *p) { return p->vx; }			///< Returns vx.
static double param_vy (const marker_t *p) { return p->vy; }			///< Returns vy.
static double param_vz (const marker_t *p) { return p->vz; }			///< Returns vz.
static double param_e  (const marker_t *p) { return gamma2eV*(sqrt (1.0 + p->vx*p->vx + p->vy*p->vy + p->vz*p->vz) - 1.0); }	///< Returns energy in [eV].
static double param_thetaZ  (const marker_t *p) { return acos (p->vz/sqrt (p->vx*p->vx + p->vy*p->vy + p->vz*p->vz + 1e-99)); }	///< Angle between speed and OZ axis.
static double param_phiZ  (const marker_t *p) { return atan2 (p->vx, p->vy); }	///< Spherical coordinate angle around OZ axis.
static double param_zero  (const marker_t *p) { return 0.0; }			///< Returns \b zero (to use on degenerated axises).

/// Speed-optimized energy caclulation (in [eV]).
static double
param_eFast (const marker_t *p)
{
  double p2 = p->vx*p->vx + p->vy*p->vy + p->vz*p->vz;

  if (p2 < mc_fastEnergyTableLimit)						// (gamma*v)^2 < 5 => look-up table.
  {
    double s = p2*mc_fastEnergyTableStepInversed;
    int i = (int) (s + 1) - 1;
    s -= i;

    return gamma2eV*(fastEnergyTable[i]*(1 - s) + fastEnergyTable[i+1]*s);
  }

  return gamma2eV*(sqrt (1.0 + p2) - 1.0); 					// Fallback to exact expression.
}

/// All defined mappers (grouped and partially initialized).
pqrVar_t pqrFrame[mc_pqrN] =
{
  {"x", param_x, 1.0, "?"},
  {"y", param_y, 1.0, "?"},
  {"z", param_z, 1.0, "?"},
  {"v<sub>x</sub>", param_vx, 1.0, "?"},
  {"v<sub>y</sub>", param_vy, 1.0, "?"},
  {"v<sub>z</sub>", param_vz, 1.0, "?"},
  {"<greek>e</greek>", param_e, 1.0, "?"},
  {"<greek>e</greek><sub>*</sub>", param_eFast, 1.0, "?"},
  {"<greek>q</greek><sub>Z</sub>", param_thetaZ, 1.0, "?"},
  {"<greek>j</greek><sub>Z</sub>", param_phiZ, 1.0, "?"},
  {"-", param_zero, 1.0, "?"}
};

// ---------------------------------------------------------------------------
/// Gets \b q/M and prints it if new (helps to select proper value of parameter in diag_distrib.cfg).
// ---------------------------------------------------------------------------
static void
register_qDivM (double qDivM)
{
  static int     N = 0;
  static double *qDivMs = NULL;

  gamma2eV = fabs (mc_CGS_m*mc_CGS_c*mc_CGS_c/(qDivM*mc_CGS_eV));		// Saves value to use later.

  for (int t = 0 ; t < N ; ++t)
    if (fabs (qDivM - qDivMs[t]) < 1e-5)
      return;

  say ("New q/M (%e) registered (mc^2 = %.4e eV)...", qDivM, gamma2eV);

  qDivMs = (double *) realloc (qDivMs, (++N)*sizeof (double));
  qDivMs[N-1] = qDivM;
}

// ---------------------------------------------------------------------------
/// Simple wrapper for the unit initalization.
// ---------------------------------------------------------------------------
static void
pqr_setUnits (int p, double unit, const char *name)
{
  if (p < 0 || p >= mc_pqrN || fabs (unit) < 1e-10 || strlen (name) >= mc_nameLength)
    error ("pqr_setUnits: bad general coordinate id (%d) or unit (%e, '%s') or too long name (> %d characters).", p, unit, name, mc_nameLength);

  pqrFrame[p].unitValue = unit;							// Sets units to use.
  strncpy (pqrFrame[p].unitName, name, mc_nameLength);
  pqrFrame[p].unitName[mc_nameLength-1] = 0;

#if mc_nameLength < 10
  Compile time check of the length of the unit name (um, fs, n_cr, etc) is failed.
#endif
}

// ---------------------------------------------------------------------------
/// Gets parameters line in the form "{@, >} '<VarName>' in [<min>, <max>], <number of bins> bins" and returns values
/// of the parameter id, range and number of bins.
// ---------------------------------------------------------------------------
static void
main_parseAxis (char *line, double *min, double *max, int *parameterID, int *numberOfBins)
{
  char *arg = line + 1, *term;

  if (! (arg = strstr (arg, "'")))						// Looks for the first '.
    error ("main_parseAxis: bad format (%s).", line);
  ++arg;

  if (! (term = strstr (arg, "'")))						// Looks for the second '.
    error ("main_parseAxis: bad format.");
  *term = 0;

  if ((*parameterID = cfg_identifyWord (arg, "x", 0, "y", 1, "z", 2, "vx", 3, "vy", 4, "vz", 5, "e", 6, "e_fast", 7, "thetaZ", 8, "phiZ", 9, mc_cfgTermGuesses)) == -1)
    error ("main_parseAxis: bad parameter (%s).", arg);
  arg = ++term;

  if (!(arg = strstr (arg, "[")))
    error ("main_parseAxis: bad format.");
  if (!(term = strstr (arg, ",")))
    error ("main_parseAxis: bad format.");
  *term = 0;
  *min = atof (++arg);
  arg = ++term;

  if (!(term = strstr (arg, "]")))
    error ("main_parseAxis: bad format.");
  *term = 0;
  *max = atof (arg);
  arg = term + 2;

  *numberOfBins = abs (atoi (arg));
}

// ---------------------------------------------------------------------------
/// Opens config-file and checks parameters after loading. \todo Make parallel version (root reads and sends).
// ---------------------------------------------------------------------------
static void
distrib_setup (void)
{
  FILE *fp = cfg_open ("diag_distrib.cfg", "rt", "tecplot.out");		// Reads parameters of the request.
  diagStart = cfg_readInt (fp);
  diagEnd = cfg_readInt (fp);
  diagStep = cfg_readInt (fp);
  int diagMicronFemtosec = cfg_readInt (fp) == 1;				// Clamps to use as array index.

  filter_qDivM = cfg_readDouble (fp);						// Filter q/M setup.

  qN = rN = 0;									// 1D run by default.
  qVar = rVar = mc_pqrN - 1;
  qMin = qMax = rMin = rMax = 0;

  main_parseAxis (cfg_readLine (fp), &pMin, &pMax, &pVar, &pN);			// Gets P-axis parameters.
  if (cfg_isOption (fp))
    main_parseAxis (cfg_readLine (fp), &qMin, &qMax, &qVar, &qN);		// Gets Q-axis parameters.
  if (cfg_isOption (fp))
    main_parseAxis (cfg_readLine (fp), &rMin, &rMax, &rVar, &rN);		// Gets R-axis parameters.

  if (cfg_isParameter (fp))
  {
    const char *word = cfg_readWord (fp);
    if (!strcmp (word, "filter:plane"))						// Reads parameters of the slicer plane.
    {
      double x = cfg_readOptDouble (fp);
      double y = cfg_readOptDouble (fp);
      double z = cfg_readOptDouble (fp);
      planeA = cfg_readOptDouble (fp);
      planeB = cfg_readOptDouble (fp);
      planeC = cfg_readOptDouble (fp);
      planeD = - x*planeA - y*planeB - z*planeC;
      planeActivated = 1;
    }
  }
  fclose (fp);

  int lastCheckPoint = sysIO_recordsTotal () - 1;
  if (diagEnd < 0 || diagEnd > lastCheckPoint)
    diagEnd = lastCheckPoint;

  if (diagStart < 0 || diagStart > diagEnd || diagStep <= 0)
    error ("diag_distrib.cfg: bad set of records.");

  if (pN <= 0)
    error ("diag_distrib.cfg: bad histogram sizes.");

  say ("%dD data set requested.", (int)((pN != 0) + (qN != 0) + (rN != 0)));

  fp = cfg_open ("tmp/distrib_to_script", "wt", __func__);
  fprintf (fp, "%d\n  That is dimensionality of the generated dataset for vDistrib.sh to call proper layout.\n", (int)((pN != 0) + (qN != 0) + (rN != 0)));
  fclose (fp);

  units_load ();								// Sets units.

  if (diagMicronFemtosec)							// Chooses units to use.
  {
    unitTime = units (mc_femtosecond);
    strncpy (unitTimeName, "fs", mc_nameLength);
    unitTimeName[mc_nameLength-1] = 0;

    pqr_setUnits (0, units (mc_micron), "<greek>m</greek>m");
    pqr_setUnits (1, units (mc_micron), "<greek>m</greek>m");
    pqr_setUnits (2, units (mc_micron), "<greek>m</greek>m");
  }
  else
  {
    unitTime = 1.0;
    strncpy (unitTimeName, "t<sub>0</sub>", mc_nameLength);
    unitTimeName[mc_nameLength-1] = 0;

    pqr_setUnits (0, 1.0, "<greek>l</greek>");
    pqr_setUnits (1, 1.0, "<greek>l</greek>");
    pqr_setUnits (2, 1.0, "<greek>l</greek>");
  }
  pqr_setUnits (3, 1.0, "c");
  pqr_setUnits (4, 1.0, "c");
  pqr_setUnits (5, 1.0, "c");
  pqr_setUnits (6, 1.0, "eV");
  pqr_setUnits (7, 1.0, "eV");
  pqr_setUnits (8, mc_pi/180.0, "<sup>o</sup>");
  pqr_setUnits (9, mc_pi/180.0, "<sup>o</sup>");

  pMin *= pqrFrame[pVar].unitValue;							// Gets limits in dimensionless units.
  pMax *= pqrFrame[pVar].unitValue;
  qMin *= pqrFrame[qVar].unitValue;
  qMax *= pqrFrame[qVar].unitValue;
  rMin *= pqrFrame[rVar].unitValue;
  rMax *= pqrFrame[rVar].unitValue;
  planeD *= units (mc_micron);								// [micron] to [r0] transformation.

  say ("Domain is:\n  - %s in [%e, %e], %d bins", pqrFrame[pVar].name, pMin, pMax, pN);
  if (qN)
    say ("  - %s in [%e, %e], %d bins", pqrFrame[qVar].name, qMin, qMax, qN);
  if (rN)
    say ("  - %s in [%e, %e], %d bins", pqrFrame[rVar].name, rMin, rMax, rN);

  int tableN = mc_fastEnergyTableLimit*1.05*mc_fastEnergyTableStepInversed + 5;	// Size of the look-up table with small extra.
  if (!(fastEnergyTable = (double *) malloc (tableN*sizeof (double))))
    error ("distrib_setup: cannot allocate look-up table (%d elements, %.3f Kb) for fast gamma.", tableN, 1.0*tableN*sizeof (double)/1024.0);
  for (int i = 0 ; i < tableN ; ++i)
  {
    fastEnergyTable[i] = sqrt (1.0 + (1.0*i)/mc_fastEnergyTableStepInversed) - 1.0;
    SAY_DEBUG ("%d => %e", i, (1.0*i)/mc_fastEnergyTableStepInversed);
  }

  say ("Look-up table for fast gamma computation is allocated, %d elements, %.3f Kb", tableN, 1.0*tableN*sizeof (double)/1024.0);

  if (planeActivated)
  {
    double alpha = planeD/(planeA*planeA + planeB*planeB + planeC*planeC);

    say ("Filtering plane is activated:\n  - n = (%.3e, %.3e, %.3e)\n  - r = (%.3e, %.3e, %.3e) [r0] or (%.3e, %.3e, %.3e) [micron].",
              planeA, planeB, planeC, alpha*planeA, alpha*planeB, alpha*planeC,
              alpha*planeA/units (mc_micron), alpha*planeB/units (mc_micron), alpha*planeC/units (mc_micron));
  }
  say ("Start to generate distribution function data file...\n");
  msg_setRefreshTime (0.2);
}

// ---------------------------------------------------------------------------
/// Removes particles which are under the slicing plane (returns number of particles left).
// ---------------------------------------------------------------------------
static int
filter_sliceOff (int N)
{
  if (!planeActivated)								// Checks if plane exists.
    return N;

  int pos = 0;
  for (int l = 0, last = N ; l < N ; ++l)
  {
    if (planeA*chunk[pos].x + planeB*chunk[pos].y + planeC*chunk[pos].z + planeD >= 0)
      ++pos;
    else
      chunk[pos] = chunk[--last];
  }

  return pos;
}

// ---------------------------------------------------------------------------
/// Processor of the chunk.
// ---------------------------------------------------------------------------
static int
distrib_loader (int N, double qDivM)
{
  register_qDivM (qDivM);

  if (fabs (filter_qDivM) < 1e10 && fabs (qDivM - filter_qDivM) > 1e-5)		// Checks if we want to skip this q/M.
    return 0;

  N = filter_sliceOff (N);							// Applies filter plane.

  const param_f pFunc = pqrFrame[pVar].value;					// Gets axis mappers.
  const param_f qFunc = pqrFrame[qVar].value;
  const param_f rFunc = pqrFrame[rVar].value;
  const double dp = (pN) ? (pMax - pMin)/pN : 1e10;				// Calculates steps of the histogram.
  const double dq = (qN) ? (qMax - qMin)/qN : 1e10;
  const double dr = (rN) ? (rMax - rMin)/rN : 1e10;

  for (int l = 0 ; l < N ; ++l, ++markerN)
  {
    double sp = (pFunc (chunk + l) - pMin)/dp + 1,
           sq = (qFunc (chunk + l) - qMin)/dq + 1,
           sr = (rFunc (chunk + l) - rMin)/dr + 1;
    int    p = (int) (sp + 1000) - 1000, q = (int)(sq + 1000) - 1000, r = (int)(sr + 1000) - 1000;
    if (p >= 0 && p <= pN + 1 && q >= 0 && q <= qN + 1 && r >= 0 && r <= rN + 1)
    {
      sp -= p;
      sq -= q;
      sr -= r;
      const double weight = fabs (chunk[l].rho);
#define mv_gist(p, q, r) diagHist[(p) + (pN + 3)*(q) + (pN + 3)*(qN + 3)*(r)]
      mv_gist(p,   q,   r) += (1 - sp)*(1 - sq)*(1 - sr)*weight;
      mv_gist(p+1, q,   r) += sp*(1 - sq)*(1 - sr)*weight;
      mv_gist(p,   q+1, r) += (1 - sp)*sq*(1 - sr)*weight;
      mv_gist(p+1, q+1, r) += sp*sq*(1 - sr)*weight;
      mv_gist(p,   q,   r+1) += (1 - sp)*(1 - sq)*sr*weight;
      mv_gist(p+1, q,   r+1) += sp*(1 - sq)*sr*weight;
      mv_gist(p,   q+1, r+1) += (1 - sp)*sq*sr*weight;
      mv_gist(p+1, q+1, r+1) += sp*sq*sr*weight;
#undef mv_gist
    }
  }

  return 0;
}

// ---------------------------------------------------------------------------
/// Opens file and processes all markers.
// ---------------------------------------------------------------------------
static int
distrib_processRecord (int record, int cpu)
{
  SAY_DEBUG ("Processing record %d/cpu %d ...", record, cpu);
  FILE *fp = cfg_open (IO_nameCpuRec(sysName_DF_full, cpu, record), "rb", __func__);// Opens file.
  while (1)
  {
    int N;
    double qDivM;

    fread (&N, sizeof (int), 1, fp);						// Reads header.
    if (!N)									// Checks terminator.
      break;
    fread (&qDivM, sizeof (double), 1, fp);

    for (int pos = 0 ; pos < N ; pos += mc_chunkSize)
    {
      int size = (pos + mc_chunkSize > N) ? N - pos : mc_chunkSize;
      fread (chunk, size, sizeof (marker_t), fp);				// Reads chunk of markers.
      if (distrib_loader (size, qDivM))						// Sends chunk to the pipeline.
        return 1;

      say_doing ("Processing data (cpu %d): %.3f%%...", cpu, 100.0*pos/N);
    }
  }
  fclose (fp);
  return 0;
}

// ---------------------------------------------------------------------------
/// Function used to terminate parallel program.
// ---------------------------------------------------------------------------
static void
main_terminator (int code)
{
  MPI_Abort (MPI_COMM_WORLD, -1);
}

// ---------------------------------------------------------------------------
/// Entry point for diag_distrib diagnostic.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
  parameter_enterMPI (argc, argv, 0, NULL);					// Initializes log-file and MPI.
  error_hookTerminator (main_terminator, cpu_here);				// Links parallel terminator.

  distrib_setup ();								// Sets parameters using config file.

  say ("Generating distribution function data file...\n");
  say ("Warning: charges are supposed to be '+e' or '-e' - no multiionized ions in energy calculations!");

  if (!(diagHist = (float*) malloc (sizeof (float)*(pN + 3)*(qN + 3)*(rN + 3))))// Allocates histogram array.
    error ("%s: cannot allocate memory.", __FILE__);

  FILE *fp = NULL;								// Output file.
  float *globalHist = NULL;
  if (!cpu_here)
  {
    // Allocates histogram array for MPI_Reduce.
    if (!(globalHist = (float*) malloc (sizeof (float)*(pN + 3)*(qN + 3)*(rN + 3))))
      error ("%s: cannot allocate memory.", __FILE__);

    fp = cfg_open ("output/distrib.dat", "wt", __FILE__);			// Root opens output file.
    fprintf (fp, "variables = p, ");
    if (qN)
      fprintf (fp, "q, ");
    if (rN)
      fprintf (fp, "r, ");
    fprintf (fp, "\"f(p,q)\"\n");

    fprintf (fp, "DATASETAUXDATA planeType=\"p=%s[%s]", pqrFrame[pVar].name, pqrFrame[pVar].unitName);
    if (qN)
      fprintf (fp, ", q=%s[%s]", pqrFrame[qVar].name, pqrFrame[qVar].unitName);
    if (rN)
      fprintf (fp, ", r=%s[%s]", pqrFrame[rVar].name, pqrFrame[rVar].unitName);
    fprintf (fp, "\"\n");
  }

  const double dp = (pN) ? (pMax - pMin)/pN : 1e10;				// Calculates steps of the histogram.
  const double dq = (qN) ? (qMax - qMin)/qN : 1e10;
  const double dr = (rN) ? (rMax - rMin)/rN : 1e10;
  for (int record = diagStart ; record <= diagEnd ; record += diagStep)
  {
    int revision, cpuN, globMarkerN;

    sysIO_parameters (record, &revision, &time, &cpuN, NULL, NULL);		// Loads DF.
    parameter_load (revision);

    memset (diagHist, 0, sizeof (float)*(pN + 3)*(qN + 3)*(rN + 3));		// Cleans histogram bins.

    markerN = 0;
    for (int cpu = cpu_here ; cpu < cpuN ; cpu += cpu_total)			// Collect input from all cpus of subset.
      if (distrib_processRecord (record, cpu))
        error ("%s: cannot process file for record %d and cpu %d on cpu %d", __FILE__, record, cpu, cpu_here);

    MPI_Reduce (diagHist, globalHist, (pN + 3)*(qN + 3)*(rN + 3), MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce (&markerN, &globMarkerN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (cpu_here)
      continue;

    // Writes histogram to the file.
    fprintf (fp, "zone t=\"t=%.3f[%s]\", i = %d", time/unitTime, unitTimeName, pN + 1);	// Prints zone header.
    if (qN)
      fprintf (fp, ", j = %d", qN + 1);
    if (rN)
      fprintf (fp, ", k = %d", rN + 1);
    fprintf (fp, ", f = point\n");

    for (int p = 0 ; p <= pN + 2 ; ++p)						// Summation pass for reduced dimensionality.
    {
      if (!rN)									// Accumulates output in r = 1 layer.
      {
        for (int q = 1 ; q <= qN + 1 ; ++q)
        {
          double sum = 0;
          for (int r = 0 ; r <= rN + 2 ; ++r)
          {
            sum += globalHist[p + q*(pN + 3) + r*(pN + 3)*(qN + 3)];
            globalHist[p + q*(pN + 3) + r*(pN + 3)*(qN + 3)] = 0;
          }
          globalHist[p + q*(pN + 3) + 1*(pN + 3)*(qN + 3)] = sum;
        }
      }

      if (!qN)									// Accumulates output in q = 1 layer.
      {
        double sum = 0;
        for (int q = 1 ; q <= qN + 1 ; ++q)
        {
          sum += globalHist[p + q*(pN + 3) + 1*(pN + 3)*(qN + 3)];
          globalHist[p + q*(pN + 3) + 1*(pN + 3)*(qN + 3)] = 0;
        }
        globalHist[p + 1*(pN + 3) + 1*(pN + 3)*(qN + 3)] = sum;
      }
    }

    const double pUnit = pqrFrame[pVar].unitValue, qUnit = pqrFrame[qVar].unitValue, rUnit = pqrFrame[rVar].unitValue;
    const double capacity = 1.0/(((pN) ? dp*pUnit : 1)*((qN) ? dq*qUnit : 1)*((rN) ? dr*rUnit : 1));
    for (int r = 1 ; r <= rN + 1 ; ++r)
      for (int q = 1 ; q <= qN + 1 ; ++q)
        for (int p = 1 ; p <= pN + 1 ; ++p)
        {
          fprintf (fp, "%.4e ", (pMin + dp*(p - 1))/pUnit);
          if (qN)
            fprintf (fp, "%.4e ", (qMin + (q - 1)*dq)/qUnit);
          if (rN)
            fprintf (fp, "%.4e ", (rMin + (r - 1)*dr)/rUnit);
          fprintf (fp, "%.4e\n", globalHist[p + q*(pN + 3) + r*(pN + 3)*(qN + 3)]*capacity);
        }

    say ("time = %.4f (%d/%d particles)", time, markerN, globMarkerN);
  }

  if (!cpu_here)
  {
    fclose (fp);								// Exit clean-up.
    free (globalHist);
  }
  free (diagHist);
  return 0;
}
