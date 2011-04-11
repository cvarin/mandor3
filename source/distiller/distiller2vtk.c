/** \file distiller2vtk.c
  * Converts txt files produced by distiller to legacy vtk format.
  *
  * That is hack to ease a transition to paraview/mayavi2/visit tools.
  */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include <unistd.h>			// For 'chdir'.
#include <dirent.h>			// For folder scanning.

#include <math.h>

#include <assert.h>

#include "log.h"
#include "misc_cfgReader.h"

/// See 'append_to_file' function.
#define MERGE_BUFFER_SIZE	(1<<10)

//    FILE *fsrc = cfg_open(src, "rt", __func__);
//    char  buffer[MERGE_BUFFER_SIZE+1];
//    while (!feof(fsrc)) {
//       buffer[0] = buffer[MERGE_BUFFER_SIZE] = 0;
//       fgets(buffer, MERGE_BUFFER_SIZE, fsrc);
//       fputs(buffer, fdest);
//    }

//    int   fetched;
//    while (!feof(fsrc)) {
//       fetched = fread(buffer, 1, MERGE_BUFFER_SIZE, fsrc);
//       fwrite(buffer, 1, fetched, fdest);
//    }


// ---------------------------------------------------------------------------
/// Appends data from the 'src' file to the 'fdest' and deletes the source.
// ---------------------------------------------------------------------------
static void
append_to_file (FILE *fdest, const char *src)
{
   FILE *fsrc = cfg_open(src, "rt", __func__);
   char  buffer[MERGE_BUFFER_SIZE+1];
   int   fetched;
   while (!feof(fsrc)) {
      fetched = fread(buffer, 1, MERGE_BUFFER_SIZE, fsrc);
      fwrite(buffer, 1, fetched, fdest);
   }
   fclose(fsrc);
   unlink(src);
}

// ---------------------------------------------------------------------------
/// Converts given file to legacy vtk format: extracts separate data chunks
/// into tmp files and than merges all of them using 'append_to_file'.
// ---------------------------------------------------------------------------
static void
convert_to_vtk (const char *name_in, const char *name_out)
{
   FILE *fin     = cfg_open(name_in,   "rt", __func__),
        *fpoints = cfg_open(".points", "wt", __func__),
        *fspeeds = cfg_open(".speeds", "wt", __func__),
        *fcharge = cfg_open(".charge", "wt", __func__);

   int  N = 0;
   char x[31], y[31], z[31], u[31], v[31], w[31], rho[31];
   while (7 == fscanf(fin, "%30s %30s %30s %30s %30s %30s %30s ",
			   x, y, z, u, v, w, rho)) {
      fprintf(fpoints, "%s %s %s\n", x, y, z);
      fprintf(fspeeds, "%s %s %s\n", u, v, w);
      fprintf(fcharge, "%s\n", rho);
      N++;
   }
   fclose(fin);
   fclose(fpoints);
   fclose(fspeeds);
   fclose(fcharge);
   fflush(NULL);

   FILE *res = cfg_open(name_out, "wt", __func__);
   fprintf(res, "# vtk DataFile Version 2.0\n"
		"Particles distribution.\n"
		"ASCII\n"
		"\n"
		"DATASET UNSTRUCTURED_GRID\n"
		"POINTS %d float\n", N);
   append_to_file(res, ".points");

   fprintf(res, "\n\n"
                "CELL_TYPES 2\n"
	        "1\n"
	        "1\n"
	        "\n"
	        "POINT_DATA %d\n"
	        "VECTORS speed float\n", N);
   append_to_file(res, ".speeds");

   fprintf(res, "\n\n"
                "SCALARS charge float 1\n"
	        "LOOKUP_TABLE default\n");
   append_to_file(res, ".charge");

   fclose(res);
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
int
main (void)
{
   const int found = chdir("./output/markers");
   assert(found == 0 && "cannot chdir to './output/markers'");

   say("Converting data to vtk");

   struct dirent **items;
   int             itemsN = scandir(".", &items, NULL, alphasort);
   for (int i = 0 ; i < itemsN ; ++i) {
      struct dirent *item = items[i];
      // Cheap test that it is a marker file created at 'distiller/main.c:65'.
      char ext[4];
      int  record, specie;
      if (3 == sscanf(item->d_name, "r%ds%d.%3s", &record, &specie, ext)
      &&  !strcmp(ext, "txt")) {
	 say_doing("%s", item->d_name);
	 convert_to_vtk(item->d_name, _("s%02dr%03d.vtk", specie, record));
      }
      free(item);
   }
   free(items);
   chdir("../..");

   say("Done");

   return EXIT_SUCCESS;
}
