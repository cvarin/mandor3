/** \file socket(avalanch).c
  * \brief Simple test of the socket layer (\b MPI). Cpu ring is created and starts receiving, than
  * test message is launched and must pass entire ring closing connections behind.
  *
  * Response on the message consists of the message forwarding and terminator sending to the sender.
  * After receiving of the terminator all senders consider socket closed. Root node starts avalanch
  * and doesn't have to propagate the message. Bodies of the message are checked to briefly catch
  * the corruption of the message.
  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <mpi.h>

#include "log.h"
#include "misc_socket.h"
#include "misc_partition.h"
#include "misc_parameters.h"

#ifdef test_main
  ERROR: ONLY ONE TEST MAY BE PLUGGED AT THE TIME AND test_main is defined already!
#endif

#define test_main() test_avalanch()

// ---------------------------------------------------------------------------
/// Writes message to the end and start of the buffer.
// ---------------------------------------------------------------------------
static void
test_packMessage (int msg, char *buffer, int size)
{
  memset (buffer, 0, size);
  memcpy (buffer, &msg, sizeof (int));
  memcpy (buffer + size - sizeof (int), &msg, sizeof (int));
}

// ---------------------------------------------------------------------------
/// Forms ring, root starts sending message up, receiver confirms message by issuing a terminator as a repost.
// ---------------------------------------------------------------------------
static void
test_avalanch (void)
{
  const int msg = 23954, term = 32562, size = 128<<10;									// Message and terminator.
  const int cpuUp = (cpu_here + 1)%cpu_total, cpuDn = (cpu_here - 1 + cpu_total)%cpu_total;
  int *fingers = (int*) calloc (cpu_total, sizeof (int));
  fingers[cpuUp] = size;												// Sockets to the neighbours in ring.
  fingers[cpuDn] = size;

  channel_t channel = mf_channel_init(2351, "test:ring");								// Connection channel.
  channel_open (&channel, mc_channel_outcome, fingers);
  say ("Avalanch test is started with message size %d bytes.", size);

  double startTime = MPI_Wtime ();

  socket_transfer (channel_findSocket (&channel, cpuUp, mc_channel_income));						// Initiates receives.
  socket_transfer (channel_findSocket (&channel, cpuDn, mc_channel_income));

  if (!cpu_here)
  {
    socket_t *s = channel_findSocket (&channel, cpuUp, mc_channel_outcome);						// Root starts an avalanch.
    test_packMessage (msg, s->buffer, size);
    socket_transfer (s);
  }

  while (1)
  {
    socket_t *s = NULL;
    while ((s = socket_seekWait (&channel, mc_channel_income)) != NULL)
    {
      if (s->cpu == cpuDn)												// Checks message.
      {
        say ("Message is incoming from %d.", s->cpu);
        socket_unlock (s);												// Releases socket.
        if (*((int*)s->buffer) != msg || *((int*)(s->buffer + size - sizeof (int))) != msg)
          error ("test_socket: message is corrupted.");
        say ("Message is correct.");

        socket_t *repost = channel_findSocket (&channel, cpuDn, mc_channel_outcome);					// Sends terminator as a repost.
        test_packMessage (term, repost->buffer, size);
        socket_transfer (repost);

        if (cpu_here)
        {
          repost = channel_findSocket (&channel, cpuUp, mc_channel_outcome);						// Propagates an avalanch.
          test_packMessage (msg, repost->buffer, size);
          socket_transfer (repost);
        }
      }
      else
      {
        say ("Message is incoming from %d.", s->cpu);
        socket_unlock (s);												// Releases socket.
        if (*((int*)s->buffer) != term || *((int*)(s->buffer + size - sizeof (int))) != term)
          error ("test_socket: terminator is corrupted.");
        say ("Terminator received.");
      }
    }

    s = NULL;
    while ((s = socket_seekWait (&channel, mc_channel_outcome)) != NULL)						// Checks that send is completed.
    {
      say ("Send to cpu %d is completed.", s->cpu);
      socket_unlock (s);												// Releases socket.
    }

    if (channel_done (&channel, mc_channel_income) && channel_done (&channel, mc_channel_outcome))
      break;
  }
  channel_close (&channel);

  double finishTime = MPI_Wtime ();
  say ("Avalanch test is completed in %.3e sec.", finishTime - startTime);
}
