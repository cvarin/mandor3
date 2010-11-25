/** \file emCap_decomposition.c
  * Parallel exchange of EM field using cached data from boundary caps (see
  * em_caps.h).
  *
  * This module supports syncronized time advancement of the EM-field in the
  * ghost cell regions and doesn't syncronize data explicitly. For function
  * to directly perform a parallel syncronization see parr_ghostCells.c.
  *
  * Data logic:
  * - Cpu domain is extended and before time step field is syncronized already
  *   (see em_caps.h).
  * - Field is cached in region big enough to use it in PIC step without any
  *   parallel exchanges (see em_caps.h).
  * - Mesh node owner is found using rule <i>"new field must depend only on
  *   the values in the <b>[0, I] x [0, J] x [0, K]</b> region.</i>
  *   So local cpu exclusively holds and send to other cpus:
  *   - E-field in <b>[0, I-1] x [0, J-1] x [0, K-1]</b> region
  *   - H-field in <b>[1, I] x [1, J] x [1, K]</b> region
  *
  * Data structures:
  * - Regions to syncronise are evaluated automatically by providing extended
  *   domain to the syncMesh_addReg() (local part is cutted off and external
  *   nodes remain).
  * - Barcode of the regions contain <b>cap id</b> to point to the cap from
  *   which prepared data must be taken.
  *
  * \note \b NaNs are used to initialize all meshes - it easily detects the
  *          usage of the uninitialized values.
  *
  * \warning <b>THERE ARE NO SYMMETRY IN THE EXCHANGE!</b> That means that if
  *          you need data from one cpu it does not mean that this cpu needs
  *          data from you - points on the boundary are assigned to one cpu
  *          only while in reality it is OK for all of them with this boundary.
  *          It makes this assigment point->cpu asymmetric.
  */

#ifndef MC_EMCAP_DECOMPOSITION_C
#define MC_EMCAP_DECOMPOSITION_C

static connection_t capConnectionH = mf_connection_init(TAG_EM_SPLITTER_H, "cap/H");	///< Connection to send/recv magnetic field.
static connection_t capConnectionE = mf_connection_init(TAG_EM_SPLITTER_E, "cap/E");	///< Connection to send/recv electric field.

// ---------------------------------------------------------------------------
/// Removes all allocated structures.
// ---------------------------------------------------------------------------
static void
capBC_parr_exit (void)
{
  syncMesh_closeConnection (&capConnectionE);										// Closes all connections.
  syncMesh_closeConnection (&capConnectionH);

  SAY_DEBUG ("capBC_parr_exit: all cleaned.");
}

// ---------------------------------------------------------------------------
/// Sets sizes of the volume caps around the boundary to use it in speculative time step and opens connections.
///
/// Algorithm: neighbour's extended domain is overlapped with region where field is defined to find out region to cache. Cached magnetic field
/// is used to advance E so regions are extended to provide all data necessary. Few point at this stage should be mentioned:
/// -# Due to creation order caps number for E and H is the same.
/// -# Extension of the magnetic field is done upward: E = f (E, H, H_{+1}) => only upper boundary conditions for H is required.
/// -# Periodic BC may refer to the lower ghost cells as well => for magnetic field mesh region is always like for main domain in transversal
///    direction and cache region is defined by E-dependence.
// ---------------------------------------------------------------------------
static void
capBC_parr_connect (connection_t *connection, const int capClaim[4][2], const reg_t *ghost)
{
   reg_t domain = {{cpu_min[0], cpu_min[1], cpu_min[2]},
                   {cpu_max[0], cpu_max[1], cpu_max[2]}};

   for (int ea = 0 ; ea < 3 ; ++ea) {
      domain.min[ea] += capClaim[cpu_bc_min[ea]] [0];
      domain.max[ea] += capClaim[cpu_bc_max[ea]] [1];
   }

   mf_reg_collapse (&domain);												// 1D/2D normalization.

   syncMesh_createConnection (connection, sizeof (vec3D_t), ghost);							// Opens connection.
   syncMesh_addReg           (connection, &domain);										// Adds region to receive.
   syncMesh_syncronize       (connection);											// Syncronizes connection.
}

// ---------------------------------------------------------------------------
/// Updates all barcodes in the array of the lists so cached field can be located easily.
// ---------------------------------------------------------------------------
static void
capBC_parr_attachCaps (cap_t *caps, const regList_t *remote)
{
   for (int cpu = 0 ; cpu < cpu_total ; ++cpu) {
      // Skips local list.
      if (cpu == cpu_here)
         continue;

      reg_t *reg = remote[cpu].list;
      for (reg_t * const end = reg + remote[cpu].N ; reg < end ; ++reg) {
         reg->barcode = -1;
         for (int b = 0 ; b < 6 ; ++b) {
            if (reg_isInside (reg, &caps[b].toFlush))
               reg->barcode = b;
         }
         ENSURE (reg->barcode >= 0,
                 "cap for reg %s is not created", reg_printRanges (reg));
      }
   }
}

// ---------------------------------------------------------------------------
/// Creates all data structures to support data exchange between CPUs.
// ---------------------------------------------------------------------------
static void
capBC_parr_init (void)
{
  capBC_parr_exit ();													// Cleans previously allocated frame.

  static const reg_t ghostE = {{0, 0, 0}, {-mc_have_x, -mc_have_y, -mc_have_z}}, ghostH = {{mc_have_x, mc_have_y, mc_have_z}, {0, 0, 0}};
  capBC_parr_connect (&capConnectionE, capClaimE, &ghostE);								// Opens connection for E.
  capBC_parr_connect (&capConnectionH, capClaimH, &ghostH);								// Opens connection for H.

  capBC_parr_attachCaps (capsE, capConnectionE.remote);									// Links send regions with cached field.
  capBC_parr_attachCaps (capsH, capConnectionH.remote);

  syncMesh_dump (&capConnectionE);											// Dumps connections.
  syncMesh_dump (&capConnectionH);
  SAY_DEBUG ("capBC_parr_init: all done.");
}

// ---------------------------------------------------------------------------
/// Calculates \f$ H^{n+1/2} \f$ for given set of nodes and sends it.
// ---------------------------------------------------------------------------
static void
capBC_parr_throw (meshVec_RO_p E, meshVec_p H)
{
   // Starts waiting for the data.
   syncMesh_irecv (&capConnectionH);
   syncMesh_irecv (&capConnectionE);

   // Gets packing lists, packs and sends H-data.
   int N;
   socket_t *s;
   mf_syncMesh_put (capConnectionH, s, N);
   for (const socket_t * const end = s + N ; s < end ; ++s) {
      const regList_t *list = (regList_t*) s->boss;
      reg_t   *reg =             list->list;
      vec3D_t *pos = (vec3D_t *) s->buffer;
      for (const reg_t * const regEnd = list->list + list->N  ; reg < regEnd ; ++reg) {
         meshVec_t *cH = &(capsH[reg->barcode].mesh);
/*      SAY_DEBUG ("Packing H-field for cpu %d: %s, cap %d.", s->cpu, reg_printRanges(reg), reg->barcode);*/
         for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
         for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
         for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos) {
            *pos = mv_f(cH, i, j, k);
/*            if (cap_isNaN (mv_fx (cH, i, j, k)) || cap_isNaN (mv_fy (cH, i, j, k)) || cap_isNaN (mv_fz (cH, i, j, k)))
              SAY_DEBUG ("packing NaN at (%d, %d, %d) / caps %d.", i, j, k, reg->barcode);*/
         }
      }
      socket_transfer (s);
   }

   // Gets packing lists, packs and sends E-data.
   mf_syncMesh_put (capConnectionE, s, N);
   for (const socket_t * const end = s + N ; s < end ; ++s) {
      const regList_t *list = (regList_t*) s->boss;
      reg_t           *reg  =              list->list;
      vec3D_t         *pos  = (vec3D_t *)  s->buffer;
      for (const reg_t * const regEnd = list->list + list->N  ; reg < regEnd ; ++reg) {
         meshVec_t *cE = &(capsE[reg->barcode].mesh);
/*      SAY_DEBUG ("Packing E-field for cpu %d: %s, cap %d.", s->cpu, reg_printRanges(reg), reg->barcode);*/
         for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
         for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
         for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos) {
            *pos = mv_f(cE, i, j, k);
/*            if (cap_isNaN(mv_fx (cE, i, j, k)) || cap_isNaN(mv_fy (cE, i, j, k)) || cap_isNaN(mv_fz (cE, i, j, k)))
              SAY_DEBUG ("packing NaN at (%d, %d, %d) / caps %d.", i, j, k, reg->barcode);*/
         }
      }
      socket_transfer (s);
   }
/*
   SAY_DEBUG ("Done\n\n");*/
}

// ---------------------------------------------------------------------------
/// Receives \f$ H^{n+1/2} \f$ for given set of nodes and installs it.
// ---------------------------------------------------------------------------
static void
capBC_parr_catchH (meshVec_p H)
{
   // Receives and updates magnetic field.
   socket_t *s;
   while ((s = socket_seekWait (& capConnectionH.channel, mc_channel_income))) {
      const regList_t *regs = (regList_t *) s->boss;
      vec3D_t         *pos  = (vec3D_t *)   s->buffer;
      reg_t           *reg  =               regs->list;
      for (const reg_t * const end = regs->list + regs->N ; reg < end ; ++reg) {
/*      SAY_DEBUG ("Getting H-field from cpu %d: %s.", s->cpu, reg_printRanges(reg));*/
         for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
         for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
         for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos) {
            mv_f(H, i, j, k) = *pos;
/*            if (cap_isNaN(pos->x) || cap_isNaN(pos->y) || cap_isNaN(pos->z))
              SAY_DEBUG ("Unpacking H-NaN at (%d, %d, %d).", i, j, k);*/
         }
      }
      socket_unlock (s);
   }
}

// ---------------------------------------------------------------------------
/// \brief Receives \f$ E^{n+1} \f$ for given set of nodes and installs it.
// ---------------------------------------------------------------------------
static void
capBC_parr_catchE (meshVec_p E)
{
  socket_t *s;
  while ((s = socket_seekWait (& capConnectionE.channel, mc_channel_income)))						// Receives and updates magnetic field.
  {
    const regList_t *regs = (regList_t *) s->boss;
    vec3D_t *pos = (vec3D_t *) s->buffer;
    reg_t *reg = regs->list;
    for (const reg_t * const end = regs->list + regs->N ; reg < end ; ++reg)
    {
/*      SAY_DEBUG ("Getting E-field from cpu %d: %s.", s->cpu, reg_printRanges(reg));*/
      for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
        for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
          for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos)
          {
            mv_f(E, i, j, k) = *pos;
/*            if (cap_isNaN(pos->x) || cap_isNaN(pos->y) || cap_isNaN(pos->z))
              SAY_DEBUG ("Unpacking E-NaN at (%d, %d, %d).", i, j, k);*/
          }
    }
    socket_unlock (s);
  }
}

// ---------------------------------------------------------------------------
/// Finalizes exchange round started in the capBC_parr_throw().
// ---------------------------------------------------------------------------
static void
capBC_parr_finishRound (void)
{
  syncMesh_waitSend (&capConnectionH);											// Finishes all sending.
  syncMesh_waitSend (&capConnectionE);

  ENSURE (syncMesh_roundIsCompleted (&capConnectionH)
       && syncMesh_roundIsCompleted (&capConnectionE),
          "all data must be sended and received at this point already");
}

#endif
