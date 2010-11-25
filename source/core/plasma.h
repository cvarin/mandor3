/** \file plasma.h
  * Boris timestep with early advance and parallel nonblocking exchange.
  *
  * All markers are stored in huge array "plasma". Array has unused space
  * in the middle, as shown below:
  *
  *  <--- shell --->     unused     <- Out -><---- core ---->
  * [xxxxxxxxxxxxxxx----------------+++++++++xxxxxxxxxxxxxxxx].
  *                 ^               ^        ^                ^
  *   countShell ---+  countSend ---+        +-- countCore    +-- countAll
  *
  * Example: total number of particles is 'countShell + countAll - countCore'
  *
  * Division details
  * ================
  *   Core   - particles far from the boundary, never contribute to other cpus;
  *          - timestep is parallel exchange independent;
  *          - no boundary conditions necessary;
  *          - outgoing particles move into the shell.
  *
  *   Shell  - boundary aligned markers, can contribute to other cpus (j, ρ);
  *          - parallel exchange: incoming particles arrive to the shell only;
  *          - outgoing particles go out to other nodes' shells or into core.
  *
  *   Out    - temporary space for outgoing particles.
  *
  *   Unused - incoming markers are added to shell in previously unused region.
  *
  *
  * Parallel timestep stages
  * ========================
  * Step 0 (initial plasma array)
  *     small character means old time,
  *     capital characters - new time,
  *     's' - shell,      'c' - core,
  *     'o' - outgoing,   'i' - incoming.
  * >>> [sssss---cccc]     // That is plasma array.
  *
  * Step 1 (shell timestep):
  *     shell particles are advanced
  *     markers going into the core are moved immediately
  * >>> [sssss---cccc] -> [SSSS---Ccccc]
  *
  * Step 2 (boundary conditions + parallel exchange)
  *     all local boundary conditions are applied,
  *     outgoing particles are moved 'out'
  *     isend for outgoing particles
  *     irecv for incoming particles directly into the shell
  * >>> [SSSS---Ccccc] -> [SS---OOCcccc] -> [SSI--OOCcccc]
  *
  * Step 3 (core timestep)
  *     core particles are advanced (ignoring fresh particles from shell)
  *     outgoing particles are placed to shell
  * >>> [SSI--OOCcccc] -> [SSIS-OO-CCCC]
  *
  * Step 4 (parallel exchange completion)
  * >>> [SSIS-OO-CCCC] -> [ssss----cccc]
  *
  * Done.
  */

#ifndef MC_PLASMA_HEADER
#define MC_PLASMA_HEADER

#include "type_mesh.h"
#include "type_marker.h"


// Particle position sometimes is compared against domain boundary to associate
// the particle with a node. To make sure that round-off effects will not turn
// a particle located exactly on the boundary into an orphan, small margin
// should be used.
//
// It should be used consistently everywhere; for example, offset in parallel
// exchange module should be smaller than the one in plasma loader.
#define MC_EPS (1e-5)


marker_t *plasma;
long int  countAll, countCore, countShell, countSend;

void plasma_init        (void);
void plasma_temperature (double *WTx, double *WTy, double *WTz);
void plasma_move        (meshVec_RO_p E, meshVec_RO_p H, meshVec_p J);

void jbc_finish (meshVec_p J);
void plasma_rho (meshDouble_p rho);

#endif
