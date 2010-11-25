/** \file parr_meshes.h
  * High level library to exchange the meshe data between cpus. Includes
  * tools to partition regions and to handle periodic wrap of the content.
  */

#ifndef MC_PARR_MESHES_HEADER
#define MC_PARR_MESHES_HEADER				///< Multiple inculdes guard.

#include "type_reg.h"
#include "misc_socket.h"

/**
  * Structure to group information for mesh syncronisation routines (channel and list of regions to fit to the socket buffers).
  */
typedef struct
{
  int 	      allocated;				///< State of the descriptor (\b open != 0 <=> in use).
  int         bytesPerNode;				///< Number of bytes per node for storage allocation.
  channel_t   channel;					///< Channel.
  regList_t  *local;					///< Local regions to process.  \todo compress.
  regList_t  *remote;					///< Remote regions to process. \todo compress.
  regList_t   map;					///< Distribution map (used for partitioning of the regions in syncMesh_addReg).
} connection_t;

/**
  * Produces initialization string for connection structure.
  */
#define mf_connection_init(Tag, NameString) {0, 0, mf_channel_init(Tag, NameString), 0, 0, mc_regList_init}

// Creation of the entire exchange structure routines.
void syncMesh_createConnection (connection_t *slot, int bytesPerNode, const reg_t *ghost);
void syncMesh_addReg (const connection_t *slot, const reg_t *reg);
void syncMesh_syncronize (connection_t *slot);
void syncMesh_verify (const connection_t *slot);
void syncMesh_dump (const connection_t *slot);
void syncMesh_closeConnection (connection_t *slot);

// Parallel exchange routines.
void syncMesh_irecv (const connection_t *slot);
void syncMesh_waitSend (const connection_t *slot);
int  syncMesh_roundIsCompleted (const connection_t *slot);

/**
  * \brief Returns pointers required for packing pass.
  */
#define mf_syncMesh_put(slot, s, N)							\
{											\
  (N) = (slot).channel.socketsN[mc_channel_outcome];					\
  (s) = (slot).channel.sockets[mc_channel_outcome];					\
}

#endif
