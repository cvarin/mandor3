/** \file misc_socket.c
  * \brief Library to encapsulate low level MPI nonblocking exchange. Ensures that correct \b tags are used,
  * all sends/recvs are matched and all temporary information is never lost.
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "log.h"
#include "misc_socket.h"
#include "misc_partition.h"

#define MC_TAGS_STACK_SIZE (100)					///< Size of the stack (to make sure consistency check is robust).

static int channelTags[MC_TAGS_STACK_SIZE], channelTagsN = 0;		///< Stack of tags to ensure uniqueness.

// ---------------------------------------------------------------------------
/// Wrappers for MPI_Alloc_mem for computers which may not have them (MPI-1 computers).
// ---------------------------------------------------------------------------
static void
socket_allocateMem (socket_t *s, long long int size, const char *chName)
{
#ifndef MC_NO_MPI_ALLOC_MEM
  MPI_Alloc_mem (size, MPI_INFO_NULL, &(s->buffer));
  ENSURE (s->buffer,
          "[%s] cannot 'MPI_Alloc_mem' %.2f Kb", chName, size/1024.0);
#else
  s->buffer = malloc (size);
  ENSURE (s->buffer,
          "[%s] cannot 'malloc' %.2f Kb", chName, size/1024.0);
  SAY_DEBUG ("%s(%s): fallback (malloc) is used to allocate %.2f Kb.", __func__, chName, size/1024.0);
#endif
}

// ---------------------------------------------------------------------------
/// Wrappers for MPI_Free_mem for computers which may not have them (MPI-1 computers).
// ---------------------------------------------------------------------------
static void
socket_freeMem (socket_t *s)
{
  assert (s->buffer);
#ifndef MC_NO_MPI_ALLOC_MEM
  MPI_Free_mem (s->buffer);
#else
  free (s->buffer);
#endif
  s->buffer = 0;
}

// ---------------------------------------------------------------------------
/// Takes array of the buffer sizes and allocates list of sockets.
// ---------------------------------------------------------------------------
static void
channel_allocateSide (channel_t *ch, const int dir, int *bufferSizes)
{
  bufferSizes[cpu_here] = 0;												// Excludes me-to-myself exchange.

  for (int cpu = 0 ; cpu < cpu_total ; cpu++)										// Counts number of sockets required.
    ch->socketsN[dir] += (bufferSizes[cpu] != 0);

  if (ch->socketsN[dir])												// Allocates empty sockets.
  {
    ch->sockets[dir] = (socket_t *) calloc (ch->socketsN[dir], sizeof (socket_t));
    ch->requests[dir] = (MPI_Request *) calloc (ch->socketsN[dir], sizeof (MPI_Request));
  }

  socket_t *s = ch->sockets[dir];
  for (int cpu = 0 ; cpu < cpu_total ; ++cpu)										// Initializes sockets.
    if (bufferSizes[cpu])
    {
      s->cpu = cpu;
      s->tag = ch->tag;
      s->direction = dir;
      s->capacity = bufferSizes[cpu];
      s->request = ch->requests[dir] + ((int)(s - ch->sockets[dir]));
      socket_allocateMem (s, bufferSizes[cpu], ch->name);
      *(s->request) = MPI_REQUEST_NULL;
      ++s;
    }

  // Allocation of dummy arrays of the Requests (so Test/Wait search functions will have no need to check against NULL pointers).
  if (! ch->requests[dir])
  {
    ch->requests[dir] = (MPI_Request*) malloc (sizeof (MPI_Request));
    ch->requests[dir][0] = MPI_REQUEST_NULL;
  }
}

// ---------------------------------------------------------------------------
/// Checks that tag is not in use.
// ---------------------------------------------------------------------------
static int
channel_registerTag (int tag)
{
  for (int l = 0 ; l < channelTagsN ; ++l)
    if (tag == channelTags[l])
      return 1;

  channelTags[channelTagsN++] = tag;
  assert (channelTagsN < MC_TAGS_STACK_SIZE);
  return 0;
}

// ---------------------------------------------------------------------------
/// Checks that tag is not in use.
// ---------------------------------------------------------------------------
static void
channel_unregisterTag (int tag)
{
  for (int l = 0 ; l < channelTagsN ; ++l)
    if (tag == channelTags[l])
    {
      channelTags[l] = channelTags[--channelTagsN];
      return;
    };

  DIE ("tag '%d' was never registered", tag);
}

// ---------------------------------------------------------------------------
/// Opens channel and initializes all fields of the channel_t structure.
// ---------------------------------------------------------------------------
void
channel_open (channel_t *ch, int direction, int *swap_me)
{
  direction = (direction != 0);												// Clamps direction flag to {0, 1}.
  ENSURE (!channel_registerTag (ch->tag),
          "[%s] tag '%d' is used already", ch->name, ch->tag);

  int *swap_they = (int*) calloc (cpu_total, sizeof (int));								// Allocates second exchange list.
  MPI_Alltoall (swap_me, 1, MPI_INT, swap_they, 1, MPI_INT, MPI_COMM_WORLD);						// Exchanges by invitations.
  channel_allocateSide (ch, direction, swap_me);									// Allocates two array of sockets/requests.
  channel_allocateSide (ch, 1 - direction, swap_they);

  free (swap_they);													// Returns memory.
  ch->open = 1;														// Marks success.
}

// ---------------------------------------------------------------------------
/// Closes channel (uses \b MPI_Barrier and few tests to ensure that exchange is completed).
// ---------------------------------------------------------------------------
void
channel_close (channel_t *ch)
{
  if (! ch->open)
  {
//     msg_warning (mf_here, "channel '%s' is not open.", ch->name);
    return;
  }

  channel_unregisterTag (ch->tag);											// Unregisters tag.
  for (int side = 0 ; side < 2 ; ++side)
  {
    socket_t *s = ch->sockets[side];											// Returns buffers to the system.
    for (int i = 0 ; i < ch->socketsN[side] ; ++i, ++s)
    {
      ENSURE (*(s->request) == MPI_REQUEST_NULL,
              "[%s] %s socket is used by cpu %d.", ch->name, ((side) ? "send" : "recv"), s->cpu);
      socket_freeMem (s);
    }
  }

  #define MF_FREE(pntr) { if (!(pntr)) free (pntr); }									/// Safe deallocation.
  MF_FREE (ch->sockets[0]);												// Removes all sockets' arrays.
  MF_FREE (ch->sockets[1]);
  MF_FREE (ch->requests[0]);
  MF_FREE (ch->requests[1]);
  ch->socketsN[0] = ch->socketsN[1] = 0;
  #undef MF_FREE

  MPI_Barrier (MPI_COMM_WORLD);
  ch->open = 0;
}

// ---------------------------------------------------------------------------
/// Activates all sockets of given type (\b recv or \b send side of the channel).
// ---------------------------------------------------------------------------
void
channel_activateSide (const channel_t *ch, int direction)
{
  socket_t *s = ch->sockets[direction];
  if (direction)
  {
    for (int i = 0 ; i < ch->socketsN[direction] ; ++i, ++s)
    {
      /// \todo Remove after implementing channel-lock flag.
      ENSURE (!s->locked, "attempt to activate locked socket");
      MPI_Isend (s->buffer, s->capacity, MPI_BYTE, s->cpu, ch->tag, MPI_COMM_WORLD, s->request);
      s->locked = 1;
    }
  }
  else
  {
    for (int i = 0 ; i < ch->socketsN[direction] ; ++i, ++s)
    {
      /// \todo Remove after implementing channel-lock flag.
      ENSURE (!s->locked, "attempt to activate locked socket");
      MPI_Irecv (s->buffer, s->capacity, MPI_BYTE, s->cpu, ch->tag, MPI_COMM_WORLD, s->request);
      s->locked = 1;
    }
  }
}

// ---------------------------------------------------------------------------
/// Returns pointer to the socket responsible for exchange (sending or receiving is defined by \b direction) with particular \b cpu.
// ---------------------------------------------------------------------------
socket_t *
channel_findSocket (const channel_t *ch, int cpu, int direction)
{
  ENSURE (cpu != cpu_here, "[%s] do not connect to yourself", ch->name);

  int N = ch->socketsN[direction];
  socket_t *s = ch->sockets[direction];
  for (int i = 0 ; i < N ; ++i, ++s)
    if (s->cpu == cpu)
      return s;

  return NULL;
}

// ---------------------------------------------------------------------------
/// Checks that all channel operations (\b direction defines type to check) are completed (namely, that all sockets are unlocked).
// ---------------------------------------------------------------------------
int
channel_done (const channel_t *ch, int direction)
{
  int N = ch->socketsN[direction];											// Chooses set of sockets to check.
  socket_t *start = ch->sockets[direction], *end = start + N;

  for ( ; start < end ; ++start)
    if (start->locked)
      return 0;

  return 1;
}

// ---------------------------------------------------------------------------
/// \brief Seeks socket which has completed operation (\b irecv or \b isend). Used to process data in the arrival order.
/// <b>Unlocking of the socket must be done by client to free the socket</b>.
// ---------------------------------------------------------------------------
socket_t *
socket_seekTest (const channel_t *ch, int direction)
{
  int done, num;
  MPI_Status status;
  MPI_Testany (ch->socketsN[direction], ch->requests[direction], &num, &done, &status);
  if (done && num != MPI_UNDEFINED)
  {
    socket_t *s = ch->sockets[direction] + num;
    ENSURE (s->locked, "MPI_Testany pointed to the unlocked socket (cpu = %d, direction = %s).", s->cpu, (s->direction) ? "outcome" : "income");
    return s;
  }

  return NULL;
}

// ---------------------------------------------------------------------------
/// \brief Waits for any socket to complete operation (\b irecv or \b isend). Used to process data in the arrival order.
/// <b>Unlocking of the socket must be done by client to free the socket</b>. For performance reason tests of the ID may be omitted.
// ---------------------------------------------------------------------------
socket_t *
socket_seekWait (const channel_t *ch, int direction)
{
  int num;
  MPI_Status status;
  MPI_Waitany (ch->socketsN[direction], ch->requests[direction], &num, &status);
  if (num != MPI_UNDEFINED)
  {
    socket_t *s = ch->sockets[direction] + num;
    ENSURE (s->locked, "MPI_Waitany pointed to the unlocked socket (cpu = %d, direction = %s).", s->cpu, (s->direction) ? "outcome" : "income");
    return s;
  }
  return NULL;
}

// ---------------------------------------------------------------------------
/// Tries to unlock the socket to access/modify buffer (signal is given on fail).
// ---------------------------------------------------------------------------
int
socket_unlock (socket_t *s)
{
  if (!s->locked)
    return 0;

  int done;
  MPI_Status status;
  MPI_Test (s->request, &done, &status);
  if (!done)
    return 1;

  s->locked = 0;
  return 0;
}

// ---------------------------------------------------------------------------
/// Starts non-blocking operation.
// ---------------------------------------------------------------------------
void
socket_transfer (socket_t *s)
{
  if (s->locked || *(s->request) != MPI_REQUEST_NULL)
    DIE ("socket (cpu %d, dir %d) is trasferring data or locked", s->cpu, s->direction);

  if (s->direction)
    MPI_Isend (s->buffer, s->capacity, MPI_BYTE, s->cpu, s->tag, MPI_COMM_WORLD, s->request);
  else
    MPI_Irecv (s->buffer, s->capacity, MPI_BYTE, s->cpu, s->tag, MPI_COMM_WORLD, s->request);
  s->locked = 1;
}
