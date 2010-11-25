/** \file parr_meshes.c
  * \brief High-level parallel exchange by mesh-variables (syncronisation of scalar/vector fields).
  *
  * Notes:
  * - <b>Explicitly ignores sending messages to themself for any node.</b> It simplifies logic of all ghost-cell
  *   partitionings (see ghostSync_setConnection() vs old variants (before <b>/331</b> record in arsenal).
  */

#include <assert.h>

#include "type_mesh.h"

#include "log.h"
#include "parr_meshes.h"
#include "parr_regLists.h"

// ---------------------------------------------------------------------------
/// Allocates connection structure and reserves tag for it.
// ---------------------------------------------------------------------------
void
syncMesh_createConnection (connection_t *slot, int bytesPerNode,
                           const reg_t *ghost)
{
  assert (bytesPerNode > 0);

  if (slot->remote)								// Closes previously opened connection.
    syncMesh_closeConnection (slot);

  slot->allocated    = 1;								// Claims the slot.
  slot->bytesPerNode = bytesPerNode;						// Saves scale for the data buffer.
  slot->local        = syncReg_createArray ();						// Creates array for all passed local regions.

  regList_t baseMap = partition_getMapping ();					// Gets base map.

  int wrap[3] = {cpu_bc_min[0] == BC_SPLITTER && dmn_bc_min[0] == BC_PERIODIC,	// Wrapping request.
                 cpu_bc_min[1] == BC_SPLITTER && dmn_bc_min[1] == BC_PERIODIC,
                 cpu_bc_min[2] == BC_SPLITTER && dmn_bc_min[2] == BC_PERIODIC};

  // Builds partitioning map.
  reg_t domain = { {dmn_min[0], dmn_min[1], dmn_min[2] },
                   {dmn_max[0], dmn_max[1], dmn_max[2] }};
  reg_buildMap (&domain, &baseMap, &slot->map, ghost, wrap);

  regList_clean (&baseMap);							// Removes temporary.
}

// ---------------------------------------------------------------------------
/// Destroys connection.
// ---------------------------------------------------------------------------
void
syncMesh_closeConnection (connection_t *slot)
{
  if (!slot->remote)
  {
//     msg_warning (mf_here, "connection is not opened.");
    return;
  }

  channel_close (&slot->channel);						// Closes channel.

  if (slot->local)			  					// Cleans all lists of regions.
    for (int cpu = 0 ; cpu < cpu_total ; ++cpu)
      regList_clean (slot->local + cpu);

  if (slot->remote)			  					// Cleans all lists of regions.
    for (int cpu = 0 ; cpu < cpu_total ; ++cpu)
      regList_clean (slot->remote + cpu);

  slot->allocated = slot->bytesPerNode = 0;					// Cleans connection structure.
  slot->local = slot->remote = 0;
  regList_clean (&slot->map);
}

// ---------------------------------------------------------------------------
/// \brief Takes region, partitions it over CPUs and prepares for syncronisation using syncReg_addClaims(). For meaning of the \b ghost see reg_distribute().
// ---------------------------------------------------------------------------
void
syncMesh_addReg (const connection_t *slot, const reg_t *reg)
{
  ENSURE (!slot->remote, "cannot add region after parallel syncronization");

  regList_t partitioned;
  reg_distributeOnMap (reg, &slot->map, &partitioned);				// Splits domain between nodes.
  syncReg_addClaims (slot->local, &partitioned);				// Adds region to the list.
  regList_clean (&partitioned);							// Cleans temporary list.
}

// ---------------------------------------------------------------------------
/// Establishes connection (all receives now will be met with proper sends).
// ---------------------------------------------------------------------------
void
syncMesh_syncronize (connection_t *slot)
{
  ENSURE (!slot->remote, "multiple syncronization is not allowed");

/*  SAY_DEBUG ("List of claims for syncronization of the %s connection.", slot->channel.name);
  for (int cpu = 0 ; cpu < cpu_total ; ++cpu)
  {
    SAY_DEBUG ("  cpu %d:", cpu);
    regList_dump (slot->local + cpu, "  - ");
  }
  */
  slot->remote = syncReg_sync (slot->local);					// Syncronizes reg_t lists.

  for (int cpu = 0 ; cpu < cpu_total ; ++cpu)					// Wraps ghost regions to align properly.
    for (int r = 0 ; r < slot->remote[cpu].N ; ++r)
      reg_unwrap (slot->remote[cpu].list + r);

  int *sizes = (int*) calloc (cpu_total, sizeof (int));				// Allocates empty sizes array.
  ENSURE (sizes, "cannot allocate memory for 'sizes'");
  for (int cpu = 0 ; cpu < cpu_total ; ++cpu)					// Calculates exchange buffer sizes.
  {
    reg_t *reg = slot->local[cpu].list, *end = slot->local[cpu].list + slot->local[cpu].N;
    for ( ; reg < end ; ++reg)
      sizes[cpu] += reg_volume (reg)*slot->bytesPerNode;
  }

  sizes[cpu_here] = 0;								// Ignores sending to myself.
  channel_open (&slot->channel, mc_channel_income, sizes);			// Opens channel.
  free (sizes);

  socket_t *s = slot->channel.sockets[mc_channel_income];			// Associates socket buffer with local reg-list.
  for (const socket_t * const end = s + slot->channel.socketsN[mc_channel_income] ; s < end ; ++s)
    s->boss = slot->local + s->cpu;

  s = slot->channel.sockets[mc_channel_outcome];				// Associates socket buffer with remote reg-list.
  for (const socket_t * const end = s + slot->channel.socketsN[mc_channel_outcome] ; s < end ; ++s)
    s->boss = slot->remote + s->cpu;

  regList_clean (&slot->map);							// Removes map (no longer needed).
}

// ---------------------------------------------------------------------------
/// \brief Checks that each cpu will have received only one entry for each node (important for current density, because of data send are not
/// assigned but added to the total current density).
// ---------------------------------------------------------------------------
void
syncMesh_verify (const connection_t *slot)
{
  for (int cpu = 0 ; cpu < cpu_total ; ++cpu)					// Fills socket look-up table.
  {
    if (cpu == cpu_here)
      continue;

    int errcode = regList_verify (slot->local + cpu);
    ENSURE (!errcode,
            "bad local regions (%s)", (errcode == 1) ? "bad vertices"
                                                     : "overlapping");
    errcode = regList_verify (slot->remote + cpu);
    ENSURE (!errcode,
            "bad remote regions (%s)", (errcode == 1) ? "bad vertices"
                                                      : "overlapping");
  }
}

// ---------------------------------------------------------------------------
/// Dumps structure of the connection created.
// ---------------------------------------------------------------------------
void
syncMesh_dump (const connection_t *slot)
{
  SAY_DEBUG ("Connection '%s' structure is:", slot->channel.name);

  for (int dir = 0 ; dir < 2 ; ++dir)
  {
    static const char action[2][20] = {"receive", "send"};
    static const char way[2][20] = {"from", "to"};
    int N = slot->channel.socketsN[dir];					// Gets packing lists.
    socket_t *s = slot->channel.sockets[dir];
    SAY_DEBUG (" %d sockets to %s data:", N, action[dir]);
    for (const socket_t * const end = s + N ; s < end ; ++s)			// Dumps info.
    {
      const regList_t *list = (regList_t*) s->boss;
      reg_t *reg = list->list;
      SAY_DEBUG ("  - %s cpu %d:", way[dir], s->cpu);
      for (const reg_t * const regEnd = list->list + list->N  ; reg < regEnd ; ++reg)
        SAY_DEBUG ("    o %s, cap %d.", reg_printRanges(reg), reg->barcode);
    }
  }
}

// ---------------------------------------------------------------------------
/// Starts non-blocking receiving of the data.
// ---------------------------------------------------------------------------
void
syncMesh_irecv (const connection_t *slot)
{
  channel_activateSide (&slot->channel, mc_channel_income);			// Initiates all receiving.
}

// ---------------------------------------------------------------------------
/// Waits for send to complete.
// ---------------------------------------------------------------------------
void
syncMesh_waitSend (const connection_t *slot)
{
  socket_t *s;
  while ((s = socket_seekWait (&slot->channel, mc_channel_outcome)))
    socket_unlock (s);
}

// ---------------------------------------------------------------------------
/// Verifies that exchange round is done - all data sended and received successfully.
// ---------------------------------------------------------------------------
int
syncMesh_roundIsCompleted (const connection_t *slot)
{
  return (channel_done (&slot->channel, mc_channel_income) && channel_done (&slot->channel, mc_channel_outcome));
}
