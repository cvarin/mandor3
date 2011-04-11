/** \file main.c
  * Simple filter: loads all markers, pipes them through filter, and saves
  * the result in ascii files.
  *
  * Filter is defined in 'distiller.cfg'; Python frontend converts it into
  * 'dll.c', recompiles this executable, and executes.
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>

#include <assert.h>

#include <math.h>

#include "misc_units.h"
#include "misc_cfgReader.h"

#include "log.h"
#include "dll.h"
#include "IO_names.h"

/// Robust comparison of floating point numbers.
#define EQUAL(a, b) (fabs((a) - (b)) < 1e-5*(fabs(a) + fabs(b)))

// Output gate for particles of given type.
typedef struct {
   double   qDivM;
   FILE    *fp;
   long int N;
} gate_t;

// Incoming plasma storage.
static marker_t *plasma  = NULL;
static long int  plasmaN = 0;

// Outgoing routes (each component to each file).
static gate_t   *gate    = NULL;
static int       gateN   = 0;

// ---------------------------------------------------------------------------
/// Defines the kind of particles (grows the gate storage if necessary).
// ---------------------------------------------------------------------------
static inline gate_t*
mapping_find (double qDivM)
{
   int i = 0;
   while (i < gateN) {
      if (EQUAL(qDivM, gate[i].qDivM))
	 return gate + i;
      ++i;
   }
   gate = (gate_t*) realloc (gate, sizeof(gate_t)*(++gateN));
   assert (gate && "cannot allocated filtering arrays");

   memset (gate + i, 0, sizeof(gate_t));
   gate[i].qDivM = qDivM;

   return gate + i;
}

// ---------------------------------------------------------------------------
/// Saves mapping between marker type and integer id.
// ---------------------------------------------------------------------------
static void
mapping_save_and_close (void)
{
   FILE *fp = cfg_open("output/markers/map.txt", "wt", __func__);
   for (int i = 0 ; i < gateN ; ++i) {
      if (gate[i].fp)
	 fclose (gate[i].fp);
      fprintf(fp, "%+.14le\n", gate[i].qDivM);
   }
   fclose(fp);

   if (gate)
      free(gate);
   gate  = NULL;
   gateN = 0;
}

// ---------------------------------------------------------------------------
/// Restores previous mapping between marker type and integer id.
// ---------------------------------------------------------------------------
static void
mapping_load (void)
{
   assert (!gate && "'gate' array is dirty (must be empty)");

   FILE *fp = fopen ("output/markers/map.txt", "rt");
   if (fp) {
      double qDivM;
      while (1 == fscanf (fp, "%le ", &qDivM)) {
         mapping_find (qDivM);
      }
      fclose (fp);
   }
}

// ---------------------------------------------------------------------------
/// Converts all raw files into a legacy vtk files.
// ---------------------------------------------------------------------------
static void
_save_into_raw_files (void)
{
   while (plasmaN) {
      const double qDivM = plasma->qDivM;
      gate_t      *g     = mapping_find (qDivM);
      if (!g->fp) {
	 const char *name = _("output/markers/.defrag_%03d", g - gate);
	 g->fp = cfg_open (name, "wb", __func__);
      }

      // Sorts all markers into the end of plasma array.
      marker_t *p = plasma,
	       *s = plasma + plasmaN;
      while (p < s) {
	 if (EQUAL (p->qDivM, qDivM)) {
	    marker_t tmp = *p;
	    *p = *(--s);
	    *s = tmp;
	    p--;
	 }
	 p++;
      }

      // Writes collected markers.
      long int N = plasmaN - (s - plasma);
      ENSURE (N, "there are must be at least the first marker to write");
      ENSURE (N == fwrite (s, sizeof(marker_t), N, g->fp),
	       "cannot write %ld markes of type %d(q/m = %le)",
	       N, g - gate, qDivM);
      plasmaN -= N;
   }
}

// ---------------------------------------------------------------------------
/// Extracts filtered and separated components from Mandor checkpoint.
// ---------------------------------------------------------------------------
static void
defrag_into_raw_files (int record, int apply_filter)
{
   mapping_load ();

   int cpu = 0;
   while (1) {
      say_doing ("rec: %2d, cpu: %2d", record, cpu);

      // Probes if checkpoint has no files left.
      FILE *fp = fopen (IO_plasmaName(record, cpu++), "rb");
      if (!fp)
	 break;

      // Loads data fragment.
      fread(&plasmaN, sizeof(long int), 1, fp);
      if (plasmaN) {
         plasma = (marker_t *) realloc (plasma, plasmaN*sizeof(marker_t));
         ENSURE (plasma, "cannot allocate %d markers.", plasmaN);
         fread  (plasma, sizeof(marker_t), plasmaN, fp);
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

      _save_into_raw_files ();
   }
   mapping_save_and_close ();
}

// ---------------------------------------------------------------------------
/// Converts raw files (one file ─ one component) into legacy vtk format.
// ---------------------------------------------------------------------------
static void
convert_raw_to_vtk (int record)
{
   mapping_load ();

   // Computes total number of markers (legacy vtk, required in advance).
   long int N = 0;
   for (int i = 0 ; i < gateN ; ++i) {
      const char *name = _("output/markers/.defrag_%03d", i);
      gate[i].fp = cfg_open (name, "rb", __func__);
      fseek (gate[i].fp, 0L, SEEK_END);
      gate[i].N = ftell (gate[i].fp)/sizeof(marker_t);
      N += gate[i].N;
   }

   FILE *vtk = cfg_open (_("output/markers/r%03d.vtk", record),
			 "wt", __func__);
   fprintf (vtk, "# vtk DataFile Version 3.0\n"
		 "Particles distribution.\n"
		 "ASCII\n"
		 "\n"
		 "DATASET POLYDATA\n"
		 "POINTS %ld float\n", N);
   marker_t p;
   for (int i = 0 ; i < gateN ; ++i) {
      rewind (gate[i].fp);
      for (long int j = 0 ; j < gate[i].N ; ++j) {
	 fread (&p, sizeof(p), 1, gate[i].fp);
	 fprintf (vtk, "%le %le %le\n", p.x, p.y, p.z);
      }
   }

   fprintf(vtk, "\n"
                "VERTICES 1 %ld %ld\n", N + 1, N);
   for (long int i = 0 ; i < N ; ++i) {
      fprintf (vtk, "%ld\n", i);
   }

   fprintf(vtk, "\n"
	        "POINT_DATA %ld\n"
	        "VECTORS speed float\n", N);
   for (int i = 0 ; i < gateN ; ++i) {
      rewind (gate[i].fp);
      for (long int j = 0 ; j < gate[i].N ; ++j) {
	 fread (&p, sizeof(p), 1, gate[i].fp);
	 fprintf (vtk, "%le %le %le\n", p.vx, p.vy, p.vz);
      }
   }

   fprintf(vtk, "\n"
	        "SCALARS charge float 1\n"
	        "LOOKUP_TABLE default\n");
   for (int i = 0 ; i < gateN ; ++i) {
      rewind (gate[i].fp);
      for (long int j = 0 ; j < gate[i].N ; ++j) {
	 fread (&p, sizeof(p), 1, gate[i].fp);
	 fprintf (vtk, "%le\n", p.rho);
      }
   }

   fprintf(vtk, "\n"
	        "SCALARS qDivM float 1\n"
	        "LOOKUP_TABLE default\n");
   for (int i = 0 ; i < gateN ; ++i) {
      for (long int j = 0 ; j < gate[i].N ; ++j) {
	 fprintf (vtk, "%le\n", gate[i].qDivM);
      }
   }

#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
   fprintf(vtk, "\n"
	        "SCALARS marker_id int 1\n"
	        "LOOKUP_TABLE default\n");
   for (int i = 0 ; i < gateN ; ++i) {
      rewind (gate[i].fp);
      for (long int j = 0 ; j < gate[i].N ; ++j) {
	 fread (&p, sizeof(p), 1, gate[i].fp);
	 fprintf (vtk, "%d\n", p.id);
      }
   }
#endif

   fclose (vtk);

   mapping_save_and_close ();
}

// ---------------------------------------------------------------------------
/// Converts raw files (one file ─ one component) into legacy vtk format.
// ---------------------------------------------------------------------------
static void
convert_raw_to_csv (int record)
{
   mapping_load ();

   for (int i = 0 ; i < gateN ; ++i) {
      const char *in   = _("output/markers/.defrag_%03d", i),
                 *out  = _("output/markers/r%03ds%02d.csv", record, i);
      FILE       *fin  = cfg_open (in,  "rb", __func__),
                 *fout = cfg_open (out, "wt", __func__);
      marker_t p;
      while (1 == fread (&p, sizeof(marker_t), 1, fin)) {
	 fprintf (fout, "% le % le % le % le % le % le %e\n",
			p.x, p.y, p.z, p.vx, p.vy, p.vz, p.rho);
      }
      fclose (fin);
      fclose (fout);
   }

   mapping_save_and_close ();
}

// ---------------------------------------------------------------------------
/// Entry point for C backend of the 'distiller' diagnostic.
// ---------------------------------------------------------------------------
int
main (int argc, char *argv[])
{
   ENSURE (argc == 5,
           "Usage: %s rec_first rec_last rec_step use_filter", argv[0]);

   int rec0 = fmax(0,           atoi(argv[1])) + 0.1,
       rec1 = fmin(INT_MAX - 3, atoi(argv[2])) + 0.1,
       step = fmax(1,           atoi(argv[3])) + 0.1,
       use_filter = atoi(argv[4]);

   // Loads [r0] to [micron] conversion factor.
   units_load ();

   for (int record = rec0 ; record <= rec1 ; record += step) {
      // Checks if record exists.
      if (access (IO_plasmaName(record, 0), R_OK))
         break;

      defrag_into_raw_files (record, use_filter);
      convert_raw_to_vtk    (record);
      convert_raw_to_csv    (record);
   }

   say("Done");

   return EXIT_SUCCESS;
}
