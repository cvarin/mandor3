/** \file main_test.c
  * \brief Entry point for small test projects (I keep name \b "test.out" only to keep \b Makefile constant).
  */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "misc_units.h"

#include "misc_parameters.h"

// #include "socket_avalanch.c"			// Avalanch test.
// #include "placer_performance.c"		// Placer test.
// #include "realToInt.c"				// Real to integer convertion test.
// #include "distribution_mapper.c"		// Distribution mapper test.

// #include "test_FFT.c"				// Distribution mapper test.

#include "log.h"
#include "type_CFile.h"

static void
createStartChunks (void)
{
  CFile_t *fp;
  fp = CF_openWrite ("output/test_chunk_%06d.dat", 10);
  CF_close (fp);

  fp = CF_openUpdate ("output/test_chunk_%06d.dat", 10);
  CF_close (fp);

  for (int chunk = 0 ; chunk < 10 ; ++chunk)
  {
    int data[10];
    for (int x = 0 ; x < 10 ; ++x)
      data[x] = x;

    say ("Filling %d", chunk);
    CFile_t *fp = CF_openUpdate ("output/test_chunk_%06d.dat", 10);
    CF_openChunk (fp, "Chunk #%02d", chunk);
    CF_write (data, sizeof (int), 10, fp);
    CF_closeChunk (fp);
    CF_close (fp);
  }
}

static void
deleteFewChunks (void)
{
  CFile_t *fp = CF_openUpdate ("output/test_chunk_%06d.dat", 10);
  for (int chunk = 1 ; chunk < 10 ; chunk += 2)
  {
    say ("Deleting with chunk %d", chunk);
    CF_deleteChunk (fp, "Chunk #%02d", chunk);
  }
  fp = CF_defrag (fp);
  CF_close (fp);
}

static void
reportChunks (void)
{
  CFile_t *fp = CF_openRead ("output/test_chunk_%06d.dat", 10);
  assert (fp);
  for (int chunk = 9 ; chunk >= 0 ; --chunk)
  {
    int data[10];
    say ("  --------  %d  -----------   ", chunk);

    if (CF_findChunk (fp, "Chunk #%02d", chunk))
    {
      say ("cannot find chunk 'Chunk #%02d'.", chunk);
      continue;
    }
    CF_read (data, sizeof (int), 10, fp);
    for (int x = 0 ; x < 10 ; ++x)
      say ("%d/%d", data[x], data[x] - x);
  }
  CF_close (fp);
}

int
main (int argc, char *argv[])
{
  parameter_enterMPI (argc, argv, 0, NULL);

  createStartChunks ();
  say ("Done making chunk :-)))");
  getchar ();

  deleteFewChunks ();
  say ("Done deleting chunks :-)))");

  reportChunks ();

  CFile_t *fp = CF_openUpdate ("output/test_chunk_%06d.dat", 10);
  for (int chunk = 1 ; chunk < 10 ; chunk += 2)
  {
    int data[10];
    for (int x = 0 ; x < 10 ; ++x)
      data[x] = 2*x;

    say ("Adding chunk %d", chunk);
    CF_openChunk (fp, "Chunk #%02d", chunk);
    CF_write (data, sizeof (int), 10, fp);
    CF_closeChunk (fp);
  }
  CF_close (fp);

  reportChunks ();

  return 0;
}
