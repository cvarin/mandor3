/** \file misc_socket.h
  * \brief Headers and data structures for non-blocking exchange using IO sockets grouped into exchange channel.
  *
  * Main purpose is envelope MPI exchange by moving MPI-tags and exchange pattern into channel_open()
  * function. After creation of connection all sends and recvs matches so one can write to socket or
  * read from sockets using simple request/test/wait functionality provided. Few global functions
  * provides ability to launch all socket operations at once.
  *
  * That is low level collective non-blocking connection library. Goal is to take local requests to/from neighbours,
  * syncronize them across complex and group all two-side point-to-point connections (sockets) to channel. Big plus
  * is that all connections dependencies are satisfied and all blocking collective calls are used only at creation
  * time. Pointers to sockets are returned on demand and shoud be linked/wrapped to remove all seek-and-return accesses.
  * Each socket can be filled with data, activated to send/recv data and probed to test if data are avaliable. Few global
  * (channel level) calls can be used to tune exchange pattern on higher level.
  *
  * \b SGI notes:
  * - Memory for exchange is allocated using \b MPI_Alloc_mem (
  *   <a href="http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi?coll=linux&db=bks&srch=&fname=/SGI_Developer/MPT_MPI_PM/sgi_html/ch03.html">
  *     Basic MPI performance recomendations</a>).
  *  That means that <i> shell variable \b SMA_GLOBAL_ALLOC must be defined in PBS script</i> (otherwise task will be terminated with message "out of memory").
  * - Shell variable \b MPI_DSM_MUSTRUN is typically defined (<a href="http://techpubs.sgi.com/library/tpl/cgi-bin/getdoc.cgi/linux/bks/SGI_Developer/books/MPT_MPI_PM/sgi_html/ch06.html#Z1028132313jlb">Memory allocation suggestion for single copy exchange and hints on mem-allocations</a>).
  *
  * \note Exchange to myself is explicitly \b ignored during channel creation in the channel_open() function.
  *
  * \todo Check: implement lock-flag for entire side (like number of locked/unlocked sockets, for example) to remove cyclic checks for channel.
  * \todo Try to replace global start <b>for (i = 0 ; i < sockN ; ++i) MPI_IRecv ()</b> on the array operation: <b>for (..) MPI_Recv_init; ... MPI_Startall ()</b>.
  *       Peoples claim in for loop you have to wait for the first byte to arrive and latency of the usual loop is N*(latency of Irecv), in other works,
  *       not max{latency_i} but summ_i*latency_i.
  */

#ifndef MC_MISC_SOCKET_HEADER
#define MC_MISC_SOCKET_HEADER

#include <mpi.h>

#define mc_channel_income	0		/**< Recv type of the socket/operation; <b>must be 0 (used as array index).</b>		*/
#define mc_channel_outcome	1		/**< Send type of the socket/operation; <b>must be 1 (used as array index).</b>		*/
#define mc_channel_NULL		(-1)		/**< Channel \b ID which will never correspond to any real channel.			*/

/**
  * Socket (elementary point-to-pont connection structure).
  */
typedef struct
{
  int          cpu;				///< Number of the cpu on the other end.
  int          tag;				///< Tag of the message (the same like in ::channel_t structure).
  int          direction;			///< Direction of the data flow.
  int          locked;				///< Flag that data are used by client and buffer cannot be reused yet.
  int          capacity;			///< Capacity of the buffer.
  char        *buffer;				///< Buffer for data (used in non-blocking exchange).
  MPI_Request *request;				///< Non-blocking exchange request hosted by channel structure.
  void        *boss;				///< Socket is just a wrapper around buffer and this field is used to associate it with external control structure.
} socket_t;

/**
  * \brief Channel - all socket grouped together by the direction and type of the data flow.
  * All channels across the complex hold exact matches for input/output sockets for each other.
  */
typedef struct
{
  const int    tag;				///< Tag used to send/recv data.
  int          open;				///< Flag to say that channel is open (in use).
  int          socketsN[2];			///< Numbers of inbox/outbox sockets in the channel.
  MPI_Request *requests[2];			///< Arrays of income/outcome requests.
  socket_t    *sockets[2];			///< Inbox/outbox sockets' array'.
  const char  *name;				///< Name of the channel (good for debug/monitor purposes).
} channel_t;

/**
  * Produces initialization string for channel structure.
  */
#define mf_channel_init(Tag, NameString) {Tag, 0, {0, 0}, {0, 0}, {0, 0}, NameString}

void channel_open (channel_t *ch, int direction, int *swap_me);
int  channel_done (const channel_t *ch, int direction);
void channel_close (channel_t *ch);
void channel_activateSide (const channel_t *ch, int direction);

socket_t* channel_findSocket (const channel_t *ch, int cpu, int direction);
socket_t* socket_seekTest (const channel_t *ch, int direction);
socket_t* socket_seekWait (const channel_t *ch, int direction);
int  socket_unlock (socket_t *s);
void socket_transfer (socket_t *s);

#define mf_socket_locked(socketPointer)		((socketPointer)->locked)	/**< Returns socket \b locked state (true for locked).	*/

#include "misc_MPItags.h"

#endif
