/** \file main.c
  * Simple filter - loads all markers, pipes them through filter, and prints
  * the result as text.
  *
  * Filter is defined in 'distiller.cfg'; Python frontend converts it into
  * 'dll.c', recompiles this executable, and executes the one.
  */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>

#include <math.h>

#include "misc_units.h"

#include "log.h"
#include "dll.h"
#include "IO_names.h"

/// Output file buffer; big one may make flushes rare and uncoherent and help
/// HDD to do the job better.
/// XXX --- test it on 3Gb dataset and remove if difference is < 5%.
#define MC_BUFFER_SIZE (1<<20)

/// Some small offset to ensure quality of floating number comparisons.
#define EPS     (1e-5)

// Output gate for particles of given type.
typedef struct {
   double  qDivM;
   FILE   *fp;
} gate_t;

// Incoming plasma storage.
static int       record  = 0;
static marker_t *plasma  = NULL;
static long int  plasmaN = 0;

// Outgoing routes (each component to each file).
static gate_t *gate  = NULL;
static int     gateN = 0;

// ---------------------------------------------------------------------------
/// Adds/reopens output stream to redirect outgoing particles by kind.
// ---------------------------------------------------------------------------
static FILE*
open_gate (double qDivM, int record)
{
   int i = 0;
   while (i < gateN && fabs (qDivM - gate[i].qDivM) >=
                                 EPS*(fabs (qDivM) + fabs (gate[i].qDivM)))
      ++i;

   if (i == gateN) {
      gate = (gate_t*) realloc (gate, (++gateN)*sizeof (gate_t));
      ENSURE (gate, "cannot allocate gate %d", i);
      gate[i].qDivM = qDivM;
      // XXX uninitialized gate[i].fp?
   }

   // XXX why do I test this?
   if (record >= 0) {
      const char *name = _("output/markers/r%03ds%02d.txt", record, i);
      // XXX replace with CFG_OPEN.
      gate[i].fp = fopen (name, "wt");
      ENSURE (gate[i].fp, "cannot open file '%s'", name);
      ENSURE (0 == setvbuf (gate[i].fp, NULL, _IOFBF, MC_BUFFER_SIZE),
              "cannot setup file buffering");
   }

   return gate[i].fp;
}

// ---------------------------------------------------------------------------
/// Finds output stream or creates one.
// ---------------------------------------------------------------------------
static inline void
save_marker (marker_t *p)
{
   FILE *fp = NULL;
   for (int i = 0 ; i < gateN && !fp ; ++i)
      if (fabs (p->qDivM - gate[i].qDivM) < EPS*(fabs (p->qDivM) +
                                                 fabs (gate[i].qDivM)))
         fp = gate[i].fp;

   if (!fp)
      fp = open_gate (p->qDivM, record);

   fprintf (fp, "%le %le %le %le %le %le %e\n",
            p->x, p->y, p->z, p->vx, p->vy, p->vz, p->rho);
}


// ---------------------------------------------------------------------------
/// Loads plasma, filters it, and sends to the output system.
// ---------------------------------------------------------------------------
static void
process_record (int apply_filter)
{
   int cpu = 0;
   while (1) {
      say_doing ("rec: %2d, cpu: %2d", record, cpu);

      // Probes file for reading.
      FILE *fp = fopen (IO_plasmaName (record, cpu++), "rb");
      if (!fp)
         break;

      // Loads data into temprorary allocated buffer.
      fread (&plasmaN, sizeof (long int), 1, fp);
      if (plasmaN) {
         plasma = (marker_t *) realloc (plasma, plasmaN*sizeof (marker_t));
         ENSURE (plasma, "cannot allocate %d markers.", plasmaN);
         fread (plasma, sizeof (marker_t), plasmaN, fp);
      }
      fclose (fp);

      // Converts units to microns.
      double factor = 1.0/units (mc_micron);
      for (int i = 0 ; i < plasmaN ; ++i) {
         plasma[i].x *= factor;
         plasma[i].y *= factor;
         plasma[i].z *= factor;
      }

      if (apply_filter)
         plasmaN = distiller_filter (plasma, plasmaN);

      for (int i = 0 ; i < plasmaN ; ++i)
         save_marker (plasma + i);
   }
}

// ---------------------------------------------------------------------------
/// Entry point for C core of \b distiller diagnostic.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
   ENSURE(argc == 5,
          "Usage: %s rec_first rec_last rec_step use_filter", argv[0]);

   int rec0 = fmax (0,           atoi (argv[1])) + 0.1,
       rec1 = fmin (INT_MAX - 3, atoi (argv[2])) + 0.1,
       step = fmax (1,           atoi (argv[3])) + 0.1,
       use_filter = atoi (argv[4]);

   // Loads units to acquire [r0] -> [micron] factor.
   units_load ();

   // Loads component parameters if they were defined before.
   FILE *fp = fopen ("output/markers/map.txt", "rt");
   if (fp) {
      double qDivM;
      while (1 == fscanf (fp, "%le ", &qDivM))
         open_gate (qDivM, -1);
      fclose (fp);
   }

   for (record = rec0 ; record <= rec1 ; record += step) {
      // Checks if record exists.
      if (access (IO_plasmaName (record, 0), R_OK))
         break;

      // Reopens streams for the components (new record => new file name).
      for (int i = 0 ; i < gateN ; ++i)
         open_gate (gate[i].qDivM, record);

      // Filters markers.
      process_record (use_filter);

      // Closes all open streams.
      for (int i = 0 ; i < gateN ; ++i)
         fclose (gate[i].fp);
   }

   // Saves mapping between 'component number' and 'q/M'.
   fp = fopen ("output/markers/map.txt", "wt");
   for (int i = 0 ; i < gateN ; ++i)
      fprintf (fp, "%+.14le\n", gate[i].qDivM);
   fclose (fp);

   say ("Done");

   return EXIT_SUCCESS;
}
