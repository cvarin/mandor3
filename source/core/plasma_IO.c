/** \file plasma_IO.c
  * Loads particles from hard-drive and sends them to the host cpus. Key
  * assumption is that HDD is order(s) of magnitude slower than parallel
  * exchange speed, so exchange traffic doesn't matter even if a cpu will load
  * wrong data to send it out.
  */

#include <mpi.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <assert.h>

#include "log.h"
#include "plasma.h"
#include "IO_names.h"
#include "misc_MPItags.h"
#include "misc_parameters.h"

/// Buffer size is 2^18 which is about 250Kb particles. Heap is increased by
/// much larger steps (5Mb of particles) to avoid many expensive reallocs.
#define BUFFER_N	(1<<18)
#define HEAP_INCREMENT	(5<<20)

// "Global" data.
static long int    *count_send = NULL,
                   *count_recv = NULL;

/// Information about particles we are loading now.
static marker_t    *buffer           = NULL;
static FILE 	   *file             = NULL;	///< File we are reading from.
static long int     particles_left   = 0,
                    particles_loaded = 0;

// ---------------------------------------------------------------------------
/// Estimates size of the plasma storage as
///   'max file size + correction for another number of the nodes',
/// and preallocates the plasma storage.
// ---------------------------------------------------------------------------
static void
_prealloc_plasma (int record)
{
   long int size = 0;

   int cpu;
   for (cpu = 0 ; cpu < INT_MAX - 10 ; ++cpu) {
      FILE *fp = fopen (IO_plasmaName (record, cpu), "rb");
      if (!fp)
         break;

      fseek (fp, 0L, SEEK_END);
      size = fmax (size, ftell (fp)/sizeof (marker_t) + 1);
      fclose (fp);
   }

   // Makes sure we didn't miss any file(s), i.e. quitted on the abscent file.
   assert (cpu < INT_MAX - 10);

   // Accounts for different number of nodes and adds small extra.
   size = 1.05*size*cpu/(double) cpu_total;

   // Allocates room for incoming data.
   if (size > countAll) {
      free (plasma);
      plasma = (marker_t*) malloc (size*sizeof (marker_t));
      ENSURE (plasma, "cannot allocate storage for %d particles", size);

      countAll = size;
   }

   // Empties 'core' and 'shell'.
   countCore  = countAll;
   countShell = 0;
}

// ---------------------------------------------------------------------------
/// Loads plasma from file (returns updated number of currently opened file).
// ---------------------------------------------------------------------------
static int
_1_open_file (int record, int file_number)
{
   if (!particles_left) {
      if (file) {
         fclose (file);
         file_number += cpu_total;
      }

      file = fopen (IO_plasmaName (record, file_number), "rb");
      if (!file) {
         particles_left = 0;
      } else {
         int loaded = fread (&particles_left, sizeof (long int), 1, file);
         ENSURE (1 == loaded, "unexpected end of file");
      }
   }

   return file_number;
}

// ---------------------------------------------------------------------------
/// Loads plasma from file (returns '1' if file exists).
// ---------------------------------------------------------------------------
static void
_2_load (int record)
{
   assert (buffer);

   particles_loaded = (particles_left <= BUFFER_N) ? particles_left
                                                   : BUFFER_N;
   if (particles_loaded) {
      int loaded = fread (buffer, sizeof (marker_t), particles_loaded, file);
      ENSURE (particles_loaded == loaded, "unexpected end of file");
   }

   particles_left -= particles_loaded;
}

// ---------------------------------------------------------------------------
/// Sorts loaded particles by host cpu: particles from [12412342341] becomes
/// [111222334444], where number is the number of host node.
// ---------------------------------------------------------------------------
static void
_3_sort (regList_t map)
{
   marker_t *start = buffer,
            *end   = buffer + particles_loaded;

   for (int host = map.N - 1 ; host >= 0 ; --host) {
      // Small offset ensures that particle on the boundary will go on one of
      // the sides, rather than stay unassigned.
      reg_t  r  = map.list[host];
      double x1 = (mc_have_x) ? (r.min[0] - 1.1*MC_EPS)*h1 : -DBL_MAX,
             y1 = (mc_have_y) ? (r.min[1] - 1.1*MC_EPS)*h2 : -DBL_MAX,
             z1 = (mc_have_z) ? (r.min[2] - 1.1*MC_EPS)*h3 : -DBL_MAX,
             x2 = (mc_have_x) ? (r.max[0] + 1.1*MC_EPS)*h1 : +DBL_MAX,
             y2 = (mc_have_y) ? (r.max[1] + 1.1*MC_EPS)*h2 : +DBL_MAX,
             z2 = (mc_have_z) ? (r.max[2] + 1.1*MC_EPS)*h3 : +DBL_MAX;

      // Groups all particles which will stay on the 'host' node.
      marker_t *tail = end,
                tmp;
      for (marker_t *p = start ; p < tail ; ++p) {
         if (p->x >= x1 && p->x <= x2
         &&  p->y >= y1 && p->y <= y2
         &&  p->z >= z1 && p->z <= z2) {
            // Swaps marker with a new untested markers taken from the tail.
            tmp    = *(--tail);
            *tail  = *p;
            *(p--) = tmp;	// Decreases 'p' to test incoming 'tmp' marker.
         }
      }

      // Counts how many particles will go to the cpu 'host'.
      count_send[host] = end - tail;
      end              = tail;
   }

   assert (start == end);
}

// ---------------------------------------------------------------------------
/// Communicates to define number of particles to exchange and prepares room
/// for incoming particles.
// ---------------------------------------------------------------------------
static void
_4_handshake (void)
{
   int res = MPI_Alltoall (count_send, 1, MPI_LONG,
                           count_recv, 1, MPI_LONG, MPI_COMM_WORLD);
   assert (res == MPI_SUCCESS);

   // Allocates room for incoming data.
   long int size = countShell;
   for (int i = 0 ; i < cpu_total ; ++i)
      size += count_recv[i];

   if (size > countAll) {
      long int new_heap = fmax (size, countAll + HEAP_INCREMENT) + 0.1;
      plasma = (marker_t*) realloc (plasma, new_heap*sizeof (marker_t));
      ENSURE (plasma, "cannot grow storage up to %d particles", new_heap);
      countAll = new_heap;
   }
}

// ---------------------------------------------------------------------------
/// Final exchange - data received directly into plasma storage.
// ---------------------------------------------------------------------------
static void
_5_exchange (void)
{
   MPI_Request wait_recv[cpu_total],
               wait_send[cpu_total];

   // Sends data using non-blocking MPI_Isend.
   marker_t *p = buffer;
   for (int i = 0 ; i < cpu_total ; ++i) {
      int res = MPI_Isend (p, count_send[i]*sizeof (marker_t), MPI_BYTE, i,
                           TAG_IO_BUFFER_DATA, MPI_COMM_WORLD, wait_recv + i);
      ENSURE (res == MPI_SUCCESS, "cannot isend %d markers to node %d",
                                  count_send[i], i);
      p += count_send[i];
   }
   assert (p == buffer + particles_loaded);

   // Starts receiving data directly into the storage.
   p = plasma + countShell;
   for (int i = 0 ; i < cpu_total ; ++i) {
      assert (p + count_recv[i] <= plasma + countAll);
      int res = MPI_Irecv (p, count_recv[i]*sizeof (marker_t), MPI_BYTE, i,
                           TAG_IO_BUFFER_DATA, MPI_COMM_WORLD, wait_send + i);
      assert (res == MPI_SUCCESS);
      p += count_recv[i];
   }
   countShell = p - plasma;

   // Completes data exchange.
   ENSURE (MPI_Waitall (cpu_total, wait_recv, MPI_STATUSES_IGNORE) == MPI_SUCCESS &&
           MPI_Waitall (cpu_total, wait_send, MPI_STATUSES_IGNORE) == MPI_SUCCESS,
           "cannot complete parallel exchange");
}

// ---------------------------------------------------------------------------
/// Loads plasma arrays and tries to have at least reserv bytes free.
// ---------------------------------------------------------------------------
void
plasma_load (int record)
{
   ENSURE (sizeof (marker_t)*BUFFER_N < INT_MAX - 10,
           "MPI_Isend@lammpi may die (too big buffer for single call)");

   // Resets plasma storage.
   _prealloc_plasma (record);

   buffer     = (marker_t *) malloc (BUFFER_N *sizeof (marker_t));
   count_send = (long int *) malloc (cpu_total*sizeof (long int));
   count_recv = (long int *) malloc (cpu_total*sizeof (long int));
   assert (buffer && count_send && count_recv);

   regList_t domain_map = partition_getMapping ();

   // Each cpu loads and distributes content of a file.
   int cpu = cpu_here;
   while (1) {
      cpu = _1_open_file (record, cpu);
      _2_load (record);
      _3_sort (domain_map);
      _4_handshake ();
      _5_exchange ();

      // Checks if all cpus have no job left.
      int job     = file || particles_left,
          the_job;
      MPI_Allreduce (&job, &the_job, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (!the_job)
         break;
   }

   assert (!file);

   free (buffer);
   free (count_send);
   free (count_recv);

   // Core is empty at this point.
   countCore = countAll;
}

// ---------------------------------------------------------------------------
/// Saves plasma arrays.
// ---------------------------------------------------------------------------
void
plasma_save (int record, int cpu)
{
    FILE *fp = fopen (IO_plasmaName (record, cpu), "wb");
    ENSURE (fp, "cannot open file '%s' for writing",
                IO_plasmaName (record, cpu));

    long int N = countAll - countCore + countShell;
    fwrite (&N, sizeof (long int), 1, fp);
    if (countShell)
        fwrite (plasma, sizeof (marker_t), countShell, fp);
    if (countAll - countCore)
        fwrite (plasma + countCore, sizeof (marker_t), countAll - countCore,
                fp);
    fclose (fp);
}
