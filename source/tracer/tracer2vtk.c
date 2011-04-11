/** \file tracer2vtk.c
  * Converts binary files 'binData/tracer_<cpu_num>.bin' produced by core.out
  * into legacy vtk format.
  *
  * That is hack to benchmark different variants of tracers.
  */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include <unistd.h>			// For 'chdir'.
#include <dirent.h>			// For folder scanning.

#include <math.h>

#include <assert.h>

#include "type_marker.h"

#include "log.h"
#include "misc_cfgReader.h"

static double tau          = NAN;

static int    input_file_N = 0,
              traces_N     = 0;

static FILE **input_file   = NULL,
            **traces_file  = NULL;

// ---------------------------------------------------------------------------
/// Converts given file to legacy vtk format: extracts separate data chunks
/// into tmp files and than merges all of them using 'append_to_file'.
// ---------------------------------------------------------------------------
static void
export_to_vtk_line (FILE *raw, FILE *vtk)
{
   // Computes number of particles in the 'raw' bag.
   fseek (raw, 0L, SEEK_END);
   long int N = ftell (raw)/sizeof(marker_t);

   fprintf (vtk, "# vtk DataFile Version 3.0\n"
		 "Particles distribution.\n"
		 "ASCII\n"
		 "\n"
		 "DATASET POLYDATA\n"
		 "POINTS %ld float\n", N);
   marker_t p;
   rewind (raw);
   for (long int i = 0 ; i < N ; ++i) {
      fread (&p, sizeof(p), 1, raw);
      fprintf (vtk, "%le %le %le\n", p.x, p.y, p.z);
   }

   fprintf(vtk, "\n"
                "LINES 1 %ld %ld\n", N + 1, N);
   for (long int i = 0 ; i < N ; ++i) {
      fprintf (vtk, "%ld\n", i);
   }

   fprintf(vtk, "\n"
	        "POINT_DATA %ld\n"
	        "VECTORS speed float\n", N);
   rewind (raw);
   for (long int i = 0 ; i < N ; ++i) {
      fread (&p, sizeof(p), 1, raw);
      fprintf (vtk, "%le %le %le\n", p.vx, p.vy, p.vz);
   }
}


// ████████████████████████████████████████████████████████████████████████ //
// ████████████████████████████████████████████████████████████████████████ //
int
main (int argc, char **argv)
{
   // Loads global parameters.
   FILE *fp = cfg_open ("binData/tracer.inf", "rt", "tracer2vtk");
   fscanf (fp, "%le %*[^\n]", &tau);
   fscanf (fp, "%d %*[^\n]",  &traces_N);
   fclose (fp);

   // Opens and counts all input streams.
   while (1) {
      FILE *fp = fopen (_("binData/tracer_%03d.bin", input_file_N), "rb");
      if (!fp)
	 break;

      input_file_N += 1;
      input_file = (FILE**) realloc (input_file, sizeof(FILE*)*input_file_N);
      assert (input_file && "cannot allocate memory for input file handlers");
      input_file[input_file_N-1] = fp;
   }
   say ("Defragmenting %d binary files into %d trajectories...", input_file_N,
								 traces_N);

   // Loads time associated with all incoming frames.
   double tmp_time[input_file_N];
   for (int cpu = 0 ; cpu < input_file_N ; ++cpu) {
      fread (tmp_time + cpu, sizeof(tmp_time[0]), 1, input_file[cpu]);
   }

   // Opens output file streams.
   traces_file = (FILE **) malloc (sizeof(FILE*)*traces_N);
   assert (traces_file && "cannot allocate array of output file pointers");
   for (int trace = 0 ; trace < traces_N ; ++trace) {
      traces_file[trace] = cfg_open (_("output/tracer/%03d.bin", trace),
				     "w+", __func__);
   }

   double cur_time = 0,
	  max_time;
   do {
      max_time = -1;
      for (int cpu = 0 ; cpu < input_file_N ; ++cpu) {
	 if (fabs(tmp_time[cpu] - cur_time) < 0.1*tau) {
	    while (1) {
	       marker_t p;
	       fread (&p, sizeof(marker_t), 1, input_file[cpu]);
	       if (!p.id)
		  break;
// 	       fprintf (traces_file[-p.id-1], "% e % e % e % e % e % e % e\n",
// 			cur_time, p.x, p.y, p.z, p.vx, p.vy, p.vz);
	       fwrite (&p, sizeof(marker_t), 1, traces_file[-p.id-1]);
	    }

	    // Update the time of the next frame.
	    int pending = fread (tmp_time + cpu, sizeof(double), 1,
				 input_file[cpu]);
	    if (!pending)
	       tmp_time[cpu] = -3*tau;
	 }
	 max_time = fmax(max_time, tmp_time[cpu]);
      }
      cur_time += tau;
   } while (max_time >= -0.1*tau);

   // Converts binary tmp files to ascii vtk format and builds the mapping.
   FILE *info = cfg_open ("output/tracer/info", "wt", __func__);
   fprintf (info, "trajectory q/M charge\n");
   for (int trace = 0 ; trace < traces_N ; ++trace) {
      marker_t p;
      rewind  (traces_file[trace]);
      fread   (&p, sizeof(marker_t), 1, traces_file[trace]);
      fprintf (info, "% 3d % le % le\n", trace, p.qDivM, p.rho);

      FILE *vtk = cfg_open (_("output/tracer/%03d.vtk", trace),
			    "wt", __func__);
      export_to_vtk_line (traces_file[trace], vtk);
      fclose (traces_file[trace]);
      fclose (vtk);
   }

   say ("Done");
   return EXIT_SUCCESS;
}
