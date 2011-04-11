/** \file plasma_parallel.c
  * Parallel exchange boundary conditions for the markers.
  *
  * Implementations of parallel particles exchange framework. Markers are sent
  * in two stages - fixed size message (number of markers) + data itself. To
  * overlap the exchange with computing, a non-blocking MPI_Isend & MPI_Irecv
  * are used.
  *
  * Interfaces:
  *   'comm_plasma_configure' analyzes domain decomposition and prepares
  *      list of gates used to send all particles to a new hosts. Should be
  *      called when domain partitioning is ready.
  *
  *   'comm_plasma_start' initiates exchange turn by posting MPI_Irecv for
  *      all incoming header messages. Should be called as soon as possible to
  *      maximize waiting period.
  *
  *   'comm_plasma_send' scans shell markers (see 'plasma.h') and isends
  *      headers and data blocks to the neighbour nodes. Outgoing data is
  *      stored in the same 'plasma' array, no additional buffer overheads.
  *
  *   'comm_plasma_test_inbox' checks if any header arrived and starts
  *      MPI_Irecv loading all incoming particles directly into the shell.
  *      When data fully arrives the shell will be up to date.
  *
  *   'comm_plasma_complete' waits for all messages to arrive, and finishes
  *      the data exchange (clears all requests, etc).
  *
  * Please note that with periodic boundary conditions things are complex - due
  * to a wrapping around a boundary, a number of neighbours may be large
  * (for 8-cpus, a single node can act as 8 neighbours at once, for example).
  */

#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>

#include "log.h"
#include "misc_MPItags.h"
#include "misc_partition.h"
#include "core/plasma.h"

/// State of the exchange.
enum { DONE = 0, WAIT_HEADER, WAIT_DATA };

// ---------------------------------------------------------------------------
/// Gate to collect particles goint to another node (see comm_plasma_send).
// ---------------------------------------------------------------------------
typedef struct
{
    double x1, y1, z1, x2, y2, z2;  ///< Region to grab the particles.
    double dx, dy, dz;		    ///< Periodic wrap correction, to add.
} gate_t;

// ---------------------------------------------------------------------------
/// Node is neighbour which has few gates through which particles may come into
/// the neighbours' domain.
// ---------------------------------------------------------------------------
typedef struct
{
    int state;		///< 0 - waiting for a header, 1 - waiting for a data.
    int cpu;		///< Cpu number of the neigbour.
    int sendN;		///< Number of markers to send to the neighbour.
    int recvN;		///< Number of markers to recv from the neighbour.

    ///< MPI requests for 2 stage (header + data) non-blocking exchange.
    MPI_Request rq_recvN, rq_recvData,
                rq_sendN, rq_sendData;

    int         gatesN;
    gate_t     *gates;
} neib_t;

// Parallel exchange framework - list of neighbours.
static int     neibsN = 0;
static neib_t *neibs  = NULL;

// ---------------------------------------------------------------------------
/// Creates parallel exchange infrastructure (note that each neighbour may have
/// few gates; e.g. 2 cpus and periodic boundary means 1 neighbour and two
/// gates).
/// XXX refactor after implementing region algebra and constructors.
// ---------------------------------------------------------------------------
void
comm_plasma_configure (void)
{
   // Composes list of neighbour nodes.
   int wrap[3] = { cpu_bc_min[0] == BC_SPLITTER && dmn_bc_min[0] == BC_PERIODIC,
                   cpu_bc_min[1] == BC_SPLITTER && dmn_bc_min[1] == BC_PERIODIC,
                   cpu_bc_min[2] == BC_SPLITTER && dmn_bc_min[2] == BC_PERIODIC };

   regList_t extMap  = mc_regList_init,
             baseMap = partition_getMapping ();

   reg_t domain = reg_vv (dmn_min, dmn_max);
   reg_buildMap  (&domain, &baseMap, &extMap, NULL, wrap);
   regList_clean (&baseMap);

   // Bbox to search for neighbours.
   reg_t bbox = { .min = { cpu_min[0] - 1,
                           cpu_min[1] - 1,
                           cpu_min[2] - 1 },
                  .max = { cpu_max[0] + 1,
                           cpu_max[1] + 1,
                           cpu_max[2] + 1 } };
   regList_t in_touch = mc_regList_init;
   reg_distributeOnMap (&bbox, &extMap, &in_touch);
   regList_clean       (&extMap);

   // Prepares exchange framework for all neigbours.
   ENSURE (neibsN == 0, "double initialization");
   neibs = (neib_t*) calloc (cpu_total, sizeof (neib_t));
   assert (neibs);

   // Sorts gates by cpu number of a receiver.
   for (int cpu = 0 ; cpu < cpu_total ; ++cpu) {
      // Searches for a gateways to a node 'cpu'.
      for (int i = 0 ; i < in_touch.N ; ++i) {
         reg_t *r = in_touch.list + i;

         // No exchange with myself.
         if (r->cpu == cpu_here || r->cpu != cpu)
               continue;

         gate_t g = { .x1 =  (r->min[0]*mc_have_x - MC_EPS)*h1,
                      .y1 =  (r->min[1]*mc_have_y - MC_EPS)*h2,
                      .z1 =  (r->min[2]*mc_have_z - MC_EPS)*h3,
                      .x2 =  (r->max[0]*mc_have_x + MC_EPS)*h1,
                      .y2 =  (r->max[1]*mc_have_y + MC_EPS)*h2,
                      .z2 =  (r->max[2]*mc_have_z + MC_EPS)*h3,
                      .dx = - r->wrap[0]*mc_have_x*h1,
                      .dy = - r->wrap[1]*mc_have_y*h2,
                      .dz = - r->wrap[2]*mc_have_z*h3 };

         // Adds structure to the list.
         neib_t *n = neibs + cpu;
         n->gates = (gate_t *) realloc (n->gates, (++n->gatesN)*
                                                   sizeof (gate_t));
         n->gates[n->gatesN-1] = g;
         n->cpu = cpu;
      }
   }

   regList_clean (&in_touch);

   // Filters only really existing neigbours.
   neibsN = 0;
   for (int i = 0 ; i < cpu_total ; ++i)
      if (neibs[i].gatesN)
         neibs[neibsN++] = neibs[i];
}

// ---------------------------------------------------------------------------
/// Starts non-blocking waiting for incoming headers.
// ---------------------------------------------------------------------------
void
comm_plasma_start (void)
{
   for (int i = 0 ; i < neibsN ; ++i) {
      neibs[i].state = WAIT_HEADER;
      MPI_Irecv (&neibs[i].recvN, 1, MPI_INT, neibs[i].cpu, TAG_MARKERS_N,
                 MPI_COMM_WORLD, &neibs[i].rq_recvN);
   }
}

// ---------------------------------------------------------------------------
/// Scans shell and schedules all outgoing particles.
// ---------------------------------------------------------------------------
void
comm_plasma_send (void)
{
   // Clears outgoing buffer.
   countSend = countCore - 1;

   // Each node grabs his particles.
   for (int i = 0 ; i < neibsN ; ++i) {
      // Remembers an end of outgoing block.
      long int origin = countSend;
      for (int j = 0 ; j < neibs[i].gatesN ; ++j) {
         gate_t *g = neibs[i].gates + j;
         for (marker_t *p = plasma ; p < plasma + countShell ; ++p) {
	    if (mc_have_x*(p->x - g->x1)*(g->x2 - p->x) >= 0
	     && mc_have_y*(p->y - g->y1)*(g->y2 - p->y) >= 0
             && mc_have_z*(p->z - g->z1)*(g->z2 - p->z) >= 0) {

	       // Optional periodic wrapping.
	       p->x += g->dx;
	       p->y += g->dy;
	       p->z += g->dz;

	       plasma[--countSend] = *p;
	       *p = plasma[--countShell];
	       --p;
	    }
         }
      }

      neibs[i].sendN = origin - countSend;
      MPI_Isend (&neibs[i].sendN, 1, MPI_INT, neibs[i].cpu,
                 TAG_MARKERS_N, MPI_COMM_WORLD, &neibs[i].rq_sendN);
      MPI_Isend (plasma + countSend, neibs[i].sendN*sizeof (marker_t),
                 MPI_BYTE, neibs[i].cpu, TAG_MARKERS_DATA,
                 MPI_COMM_WORLD, &neibs[i].rq_sendData);
   }
}

// ---------------------------------------------------------------------------
/// Checks if there are header(s) arrived and arranges data receiving process.
// ---------------------------------------------------------------------------
void
comm_plasma_test_inbox (void)
{
   int done;
   for (int i = 0 ; i < neibsN ; ++i) {
      if (neibs[i].state == WAIT_HEADER) {
         MPI_Test (&neibs[i].rq_recvN, &done, MPI_STATUS_IGNORE);
         if (done) {
            neibs[i].state = WAIT_DATA;
            MPI_Irecv (plasma + countShell,
                       neibs[i].recvN*sizeof (marker_t), MPI_BYTE,
                       neibs[i].cpu, TAG_MARKERS_DATA,
                       MPI_COMM_WORLD, &neibs[i].rq_recvData);
            countShell += neibs[i].recvN;
            ENSURE (countShell + 3 < countSend && countSend <= countCore,
                    "no room to receive %d markers from cpu %d",
                    neibs[i].recvN, neibs[i].cpu);
         }
      }

      if (neibs[i].state == WAIT_DATA) {
         MPI_Test (&neibs[i].rq_recvData, &done, MPI_STATUS_IGNORE);
         if (done)
            neibs[i].state = DONE;
      }
   }
}

// ---------------------------------------------------------------------------
/// Completes parallel data exchange turn.
// ---------------------------------------------------------------------------
void
comm_plasma_complete (int nanosecs)
{
   int waiting_for_isend = 1;
   while (waiting_for_isend) {
      waiting_for_isend = 0;
      for (int i = 0 ; i < neibsN ; ++i) {
         if (neibs[i].state != DONE)
            comm_plasma_test_inbox ();

         int sendN, sendData;
         MPI_Test (&neibs[i].rq_sendN,    &sendN,    MPI_STATUS_IGNORE);
         MPI_Test (&neibs[i].rq_sendData, &sendData, MPI_STATUS_IGNORE);

         // Sleeps to reduce cpu pressure and help MPI.
         if (!sendN || !sendData || neibs[i].state != DONE) {
            struct timespec delay = { .tv_sec  = 0,
                                      .tv_nsec = nanosecs };
            nanosleep (&delay, NULL);
            waiting_for_isend = 1;
         }
      }
   }
}

// ---------------------------------------------------------------------------
/// Dumps parallel exchange frame.
// ---------------------------------------------------------------------------
void
comm_plasma_dump (void)
{
   SAY_DEBUG ("Plasma exchange gates:");
   for (int i = 0 ; i < neibsN ; ++i) {
      SAY_DEBUG ("  o node %d:", neibs[i].cpu);
      for (int j = 0 ; j < neibs[i].gatesN ; ++j) {
         gate_t g = neibs[i].gates[j];
         SAY_DEBUG ("      [%.3e, %.3e, %.3e] - [%.3e, %.3e, %.3e] "
                                             "// [%.3e, %.3e, %.3e] ",
                    g.x1, g.y1, g.z1,
                    g.x2, g.y2, g.z2,
                    g.dx, g.dy, g.dz);
      }
   }
}
