/** \file parr_regLists.c
  * \brief Parallel syncronisation of the regList_t arrays over cluster.
  *
  * This library is used to collect requests (in form of regList_t elements) on the regions and than to syncronize all calls - each owner of the
  * claimed region will be given an information about who did request it and how to send it. This task is very common for parallel decompositions
  * of the ghost cells domains, parallel input/output and so on.
  */

#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "misc_MPItags.h"

#include "log.h"
#include "type_reg.h"

#include "misc_parameters.h"						///< \todo Replace by \b dimensions.h.

// ---------------------------------------------------------------------------
/// Creates empty array of the regList_t (\b cpu_total elements total).
// ---------------------------------------------------------------------------
regList_t*
syncReg_createArray (void)
{
  regList_t *array = (regList_t *) calloc (cpu_total, sizeof (regList_t));
  ENSURE (array, "out of memory");
  return array;
}

// ---------------------------------------------------------------------------
/// Adds claims to the syncronisation list.
// ---------------------------------------------------------------------------
void
syncReg_addClaims (regList_t* list, const regList_t *claims)
{
  reg_t *reg = claims->list;
  for (const reg_t * const end = claims->list + claims->N ; reg < end ; ++reg)
    regList_add (list + reg->cpu, reg);
}

// ---------------------------------------------------------------------------
/// \brief Takes requests generated on local side and distributes regions to the owners for futher processing.
/// Returns list of the requests received from the other members of cluster. Exchange is done in 2 rounds (syncronization of the number of regions and
/// exchange by the bodies themself). <b>Order of the regions in the each list is conserved</b>
// ---------------------------------------------------------------------------
regList_t*
syncReg_sync (regList_t *local)
{
  int *localN, *remoteN;												// Number of calls on local/remote side.
  MPI_Request *requests;												// Requests for sync routines.
  MPI_Status *statuses;													// Statuses for sync routines.
  regList_t *remote;

  localN   = (int*)          calloc (cpu_total, sizeof (int));									// Memory allocation.
  remoteN  = (int*)          calloc (cpu_total, sizeof (int));
  requests = (MPI_Request *) malloc (cpu_total*sizeof (MPI_Request));
  statuses = (MPI_Status *)  calloc (cpu_total, sizeof (MPI_Status));
  ENSURE (localN && remoteN && requests && statuses,
          "cannot allocate memory (cpu_total = %d)", cpu_total);
  remote = syncReg_createArray ();

  // Syncronzation: round one.
  for (int cpu = 0 ; cpu < cpu_total ; ++cpu)										// Prepares number of records ready for exchange.
    localN[cpu] = local[cpu].N;

  MPI_Alltoall (localN, 1, MPI_INT, remoteN, 1, MPI_INT, MPI_COMM_WORLD);						// First round: number of claims syncronisation.

  // Syncronzation: round two.
  for (int cpu = 0 ; cpu < cpu_total ; ++cpu)										// Starts non-blocking receiving.
  {
    requests[cpu] = MPI_REQUEST_NULL;
    if (remoteN[cpu])
    {
      remote[cpu].N = remoteN[cpu];
      remote[cpu].list = (reg_t*) malloc (remoteN[cpu]*sizeof (reg_t));
      MPI_Irecv (remote[cpu].list, remoteN[cpu]*sizeof (reg_t), MPI_BYTE, cpu, TAG_SYNC_REGLISTS, MPI_COMM_WORLD, requests + cpu);
    }
  }

  for (int cpu = 0 ; cpu < cpu_total ; ++cpu)										// Sends data (matching irecvs are posted).
    if (localN[cpu])
      MPI_Send (local[cpu].list, localN[cpu]*sizeof (reg_t), MPI_BYTE, cpu, TAG_SYNC_REGLISTS, MPI_COMM_WORLD);

  MPI_Waitall (cpu_total, requests, statuses);										// Waits for the end of second round.

  free (remoteN);													// Releases memory.
  free (localN);
  free (requests);
  free (statuses);

  return remote;
}
