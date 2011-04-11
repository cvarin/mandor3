/** \file mark_to_trace.c
  * Takes given record and marks requested particles (making their id negative)
  * for particle tracer.
  *
  * Currently it is a hack to evaluate different variants of tracers.
  */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include <math.h>

#include <assert.h>

#include "type_marker.h"

#include "log.h"
#include "misc_cfgReader.h"

static double tau       = NAN;
static int    traces_N  = 0,
              to_mark_N = 0,
             *to_mark   = NULL;


// ---------------------------------------------------------------------------
/// Comparison function to sort the 'to_mark' array.
// ---------------------------------------------------------------------------
static int
cmp_int (const void *a, const void *b)
{
   int x = *(int*) a,
       y = *(int*) b;
   return (x > y) - (x < y);
}

// ---------------------------------------------------------------------------
/// Takes given checkpoint file and marks the id of selected markers.
// ---------------------------------------------------------------------------
static void
mark_file (FILE *fp)
{
   rewind (fp);

   long int N;
   fread (&N, sizeof(N), 1, fp);

   qsort (to_mark, to_mark_N, sizeof(to_mark[0]), cmp_int);
   while (N--) {
      marker_t p;
      fread (&p, sizeof(p), 1, fp);
      if (bsearch (&p.id, to_mark, to_mark_N, sizeof(to_mark[0]), cmp_int)) {
	 p.id = -(++traces_N);
	 fseek  (fp, -sizeof(p), SEEK_CUR);
	 fwrite (&p,  sizeof(p), 1, fp);
      }
   }
}

// ---------------------------------------------------------------------------
/// Parses paraview generated csv file and loads the ids of selected markers.
// ---------------------------------------------------------------------------
static void
load_id_from_csv_file (const char *filename)
{
   FILE *fp = cfg_open (filename, "rt", __func__);

   // Parses the header.
   char header[1000+1] = {[1000] = 0};
   ENSURE (fscanf (fp, "%1000[^\n] ", header), "cannot read header");
   char *id_field = strstr(header, "\"marker_id\"");
   ENSURE(id_field, "cannot find 'marker_id' field in file '%s'", filename);

   // Generates format string which ignores all field but the 'marker_id'.
   const char *format = "";
   for (char *p = header ; p < id_field ; ++p)
      if (*p == ',')
	 format = _("%%*le,%s", format);
   format = _("%s%%d%%*[^\n] ", format);

   // Loads ids.
   int id;
   while (1 == fscanf(fp, format, &id)) {
      if (id > 0) {
	 to_mark = (int*) realloc(to_mark, (to_mark_N + 1)*sizeof(int));
	 to_mark[to_mark_N++] = id;
      }
   }

   fclose (fp);
}

// ████████████████████████████████████████████████████████████████████████████
// ████████████████████████████████████████████████████████████████████████████
int
main (int argc, char **argv)
{
   if (argc <= 1) {
      say ("Usage: %s file1.csv [file2.csv] ...", argv[0]);
      return EXIT_FAILURE;
   }

#if !(defined(ACTIVATE_TRACER) && ACTIVATE_TRACER)
   SAY_WARNING ("Tracer is disabled at compile time.\n"
                "What can you do:\n"
                "  + update the Makefile (uncomment/add -DACTIVATE_TRACER=1)\n"
                "  + recompile the code\n"
                "  + repeat your run\n"
                "  + use this diagnostic again\n"
                "\n"
                "Good luck :-)");
   return EXIT_FAILURE;
#endif

   // Loads global parameters.
   FILE *fp = cfg_open ("binData/tracer.inf", "rt", "tracer2vtk");
   fscanf (fp, "%le %*[^\n]", &tau);
   fscanf (fp, "%d  %*[^\n]", &traces_N);
   fclose (fp);

   // Loads all IDs we want to mark.
   for (int f = 1 ; f < argc ; ++f) {
      say                   ("File: %s", argv[f]);
      load_id_from_csv_file (argv[f]);
   }

   // Opens and counts all input streams.
   int rec         = 0,
       old_tracesN = traces_N;
   for (int cpu = 0 ; ; ++cpu) {
      const char *name = _("binData/plasma_%06d_%03d.bin", rec, cpu);
      FILE       *fp   = fopen (name, "r+");
      if (!fp) break;

      say_doing ("File: %s", name);
      mark_file (fp);
      fclose    (fp);
   }
   say ("new trajectories: %d", traces_N - old_tracesN);

   // Saves updated number of selected markers.
   fp = cfg_open ("binData/tracer.inf", "wt", "mark_trajectory");
   fprintf (fp, "%le \tτ (time step)\n",        tau);
   fprintf (fp, "%d  \tnumber of trajectories", traces_N);
   fclose  (fp);

   say ("Done");
   return EXIT_SUCCESS;
}
