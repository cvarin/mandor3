/** \file parr_ghostCells.c
  * Syncronises content of the ghost-cells (they overlap due to domain decomposition).
  *
  * Few notes:
  * - Two channels are used because of regions for E and H fields are different and may engage different set of cpus.
  *   Alternative is to use barcode to distingvish regions for E and H but that means another conditional branch in the
  *   loop which is already not optimal (ghost cells are very thin).
  *
  * - It is \b syncronization, not <b>full update</b>: content of the local ghost cells must be updated by
  *   using em-module's ebc_<> group of functions. Here only parallel aspect is adressed.
  *
  * - During initialization entire extended region is given to the syncMesh module (see parr_meshes.c). It ignores
  *   sending to themself so after the syncronization only external (ghost) cells remain.
  *
  * \warning <b>THERE ARE NO SYMMETRY IN THE EXCHANGE!</b> That means that if you need data from one cpu it
  *       does not mean that this cpu needs data from you - points on the boundary are assigned to
  *       one cpu only while in reality it is OK for all of them with this boundary. It makes this
  *       assigment point->cpu asymmetric.
  */

#include <assert.h>

#include "log.h"
#include "parr_meshes.h"

#include "core/em_caps.h"

static connection_t syncGhostE = mf_connection_init(TAG_GHOST_SYNC_E, "syncGhost:E");	///< Connection for ghost cell E updates.
static connection_t syncGhostH = mf_connection_init(TAG_GHOST_SYNC_H, "syncGhost:H");	///< Connection for ghost cell H updates.

// ---------------------------------------------------------------------------
/// Opens connections for E and H fields.
// ---------------------------------------------------------------------------
static void
ghostSync_setConnection (void)
{
  // Tuners of the real domain of defined fields.
  static const reg_t ghostE = {{0, 0, 0}, {-1, -1, -1}},
                     ghostH = {{1, 1, 1}, { 0,  0,  0}};

  syncMesh_closeConnection (&syncGhostE);										// Closes connections.
  syncMesh_closeConnection (&syncGhostH);

  syncMesh_createConnection (&syncGhostE, sizeof (vec3D_t), &ghostE);							// Opens connections.
  syncMesh_createConnection (&syncGhostH, sizeof (vec3D_t), &ghostH);

  const reg_t reg = {{cpu_min[0] - mc_have_x, cpu_min[1] - mc_have_y, cpu_min[2] - mc_have_z},
                     {cpu_max[0] + mc_have_x, cpu_max[1] + mc_have_y, cpu_max[2] + mc_have_z}};

  syncMesh_addReg (&syncGhostE, &reg);											// Gets ghost cells assigned.
  syncMesh_addReg (&syncGhostH, &reg);

  syncMesh_syncronize (&syncGhostE);											// Creates frameworks.
  syncMesh_syncronize (&syncGhostH);
/*
  syncMesh_dump (&syncGhostE);												// Dumps frameworks.
  syncMesh_dump (&syncGhostH);*/
}

// ---------------------------------------------------------------------------
/// Sends ghost cell regions.
// ---------------------------------------------------------------------------
static void
ghostSync_send (meshVec_p F, connection_t *connection)
{
//   SAY_DEBUG ("ghostSync_send: sending %s mesh", F->name);
  int N;
  socket_t *s;
  mf_syncMesh_put(*connection, s, N);											// Gets packing lists.
  for (const socket_t * const end = s + N ; s < end ; ++s)								// Packs and sends F-data.
  {
    const regList_t *list = (regList_t*) s->boss;
    reg_t *reg = list->list;
    reg_t * const regEnd = list->list + list->N;
    vec3D_t *pos = (vec3D_t *) s->buffer;
    for ( ; reg < regEnd ; ++reg)
    {
//       SAY_DEBUG ("  o %s region to cpu %d", reg_printRanges (reg), s->cpu);
      for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)								/// \todo Unified reg-cycles for 2D/1D.
        for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
          for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos)
          {
            pos->x = mv_fx(F, i, j, k);
            pos->y = mv_fy(F, i, j, k);
            pos->z = mv_fz(F, i, j, k);
          }
    }
    socket_transfer (s);
  }
//   SAY_DEBUG ("  done.");
}

// ---------------------------------------------------------------------------
/// Receives ghost cell regions.
// ---------------------------------------------------------------------------
static void
ghostSync_recv (meshVec_p F, connection_t *connection)
{
//   SAY_DEBUG ("ghostSync_recv: receiving %s mesh", F->name);
  socket_t *s;
  while ((s = socket_seekWait (&connection->channel, mc_channel_income)))						// Receives and updates F-field.
  {
    const regList_t *regs = (regList_t *) s->boss;
    vec3D_t *pos = (vec3D_t *) s->buffer;
    reg_t *reg = regs->list;
    const reg_t *end = regs->list + regs->N;
    for ( ; reg < end ; ++reg)
    {
//       SAY_DEBUG ("  o %s region from cpu %d", reg_printRanges (reg), s->cpu);
      for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
        for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
          for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos)
          {
            mv_fx(F, i, j, k) = pos->x;
            mv_fy(F, i, j, k) = pos->y;
            mv_fz(F, i, j, k) = pos->z;
          }
    }
    socket_unlock (s);
  }
//   SAY_DEBUG ("  done.");
}

// ---------------------------------------------------------------------------
/// Syncronizes overlapped ghost cells.
// ---------------------------------------------------------------------------
void
ghostSync_sync (meshVec_p E, meshVec_p H)
{
  assert (!syncGhostE.channel.open && !syncGhostH.channel.open);							// Everything should be empty and closed.

  ghostSync_setConnection ();												// Sets connection.

  syncMesh_irecv (&syncGhostH);												// Starts non-blocking receives.
  syncMesh_irecv (&syncGhostE);

  ghostSync_send (H, &syncGhostH);											// Sends fields.
  ghostSync_send (E, &syncGhostE);

  ghostSync_recv (H, &syncGhostH);											// Receives fields.
  ghostSync_recv (E, &syncGhostE);

  syncMesh_waitSend (&syncGhostH);											// Finishes all sending.
  syncMesh_waitSend (&syncGhostE);

  ENSURE (syncMesh_roundIsCompleted (&syncGhostE)
       && syncMesh_roundIsCompleted (&syncGhostH),
         "all data expected be sended and received at this point");

  SAY_DEBUG ("ghost cells are syncronised\ncontent of the caps is ignored!");

  syncMesh_closeConnection (&syncGhostE);										// Closes connections.
  syncMesh_closeConnection (&syncGhostH);
}
