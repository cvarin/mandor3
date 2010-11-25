#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <assert.h>

#include "type_marker.h"

#include "misc_PIC.h"
#include "log.h"

#include "profiler.h"
#include "parr_meshes.h"

#include "plasma_VSP.h"

/// Connection for charge density parallel exchange.
static connection_t connectionRho = mf_connection_init (TAG_PLASMA_RHO, "plasma:rho");

// ---------------------------------------------------------------------------
/// Applies boundary conditions on the charge density (before parallel exchange).
// ---------------------------------------------------------------------------
static void
plasmaRho_applyBoundary (meshDouble_p rho)
{
  // Packed parameters for BC sweeps.
  const int sizes[3] = {cpu_max[0] - cpu_min[0], cpu_max[1] - cpu_min[1], cpu_max[2] - cpu_min[2]};
  const int strides[3] = {mc_have_x*rho->width_yz, mc_have_y*rho->width_z, mc_have_z};
  const int BC[6] = {cpu_bc_min[0], cpu_bc_max[0], cpu_bc_min[1], cpu_bc_max[1], cpu_bc_min[2], cpu_bc_max[2]};

  static const int axis[6][3] = { {mc_x, mc_y, mc_z}, {mc_x, mc_y, mc_z}, {mc_y, mc_x, mc_z}, {mc_y, mc_x, mc_z}, {mc_z, mc_x, mc_y}, {mc_z, mc_x, mc_y} };
  static const int top[6] = {0, 1, 0, 1, 0, 1};

  for (int b = 0 ; b < 6 ; ++b)							// Goes through all boundaries.
  {
    const int ea = axis[b][0], ep = axis[b][1], eq = axis[b][2], dp = sizes[ep] + 4*ACTIVATOR[ep], dq = sizes[eq] + 4*ACTIVATOR[eq];
    const int stepOut = strides[ea]*(2*top[b] - 1);

    if (!ACTIVATOR[ea])								// Skips nonactive axises.
      continue;

    // Gets the offset of the point in the origin if the separation plane.
    int pos0 = mf_offset(rho, cpu_min[0], cpu_min[1], cpu_min[2]) - 2*ACTIVATOR[ep]*strides[ep] - 2*ACTIVATOR[eq]*strides[eq] + top[b]*strides[ea]*sizes[ea];

    switch (BC[b])								// Chooses BC method.
    {
      case BC_PERIODIC:
      {
        const int shift = strides[ea]*sizes[ea];

        if (!top[b])								// Periodic BC are processed @ bottom boundary.
        {
          for (int p = 0 ; p <= dp ; p++, pos0 += strides[ep])
          {
            int pos = pos0;
            for (int q = 0 ; q <= dq ; q++, pos += strides[eq])
            {
              int pos2 = pos - 2*strides[ea];
              for (int r = 0 ; r < 5 ; r++, pos2 += strides[ea])
              {
                mv_unrl (rho, pos2 + shift) += mv_unrl (rho, pos2);
                mv_unrl (rho, pos2) = mv_unrl (rho, pos2 + shift);
              }
            }
          }
        }
      }
      break;

      case BC_MIRROR:
      case BC_OPEN:
      {
        for (int p = 0 ; p <= dp ; p++, pos0 += strides[ep])
        {
          int pos = pos0;
          for (int q = 0 ; q <= dq ; q++, pos += strides[eq])
          {
            mv_unrl (rho, pos - stepOut) += mv_unrl (rho, pos + stepOut);
            mv_unrl (rho, pos - 2*stepOut) += mv_unrl (rho, pos + 2*stepOut);

            // Adds mirror charges.
            mv_unrl (rho, pos) = 0;
            mv_unrl (rho, pos + stepOut) = - mv_unrl (rho, pos - stepOut);
            mv_unrl (rho, pos + 2*stepOut) = - mv_unrl (rho, pos - 2*stepOut);
          }
        }
      }
      break;

      default:
        DIE ("unsupported boundary condition");
      case BC_SPLITTER:
        break;
    }
  }
}

// ---------------------------------------------------------------------------
/// Starts all non-blocking exchanges of charge density.
// ---------------------------------------------------------------------------
void
plasmaRho_startExchange (meshDouble_p rho)
{
  syncMesh_irecv (&connectionRho);						// Starts non-blocking receive.

  int N;
  socket_t *s;
  mf_syncMesh_put(connectionRho, s, N);						// Gets packing lists.
  for (const socket_t * const end = s + N ; s < end ; ++s)			// Packs and sends charge density.
  {
    double *pos = (double *) s->buffer;
    const regList_t *list = (regList_t*) s->boss;
    reg_t *reg = list->list;
    for (const reg_t * const regEnd = list->list + list->N ; reg < regEnd ; ++reg)
    {
//       SAY_DEBUG ("Sending rho to cpu %d: %s.", s->cpu, reg_printRanges(reg));
      for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
        for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
          for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos)
            *pos = mv_f (rho, i, j, k);
    }
    socket_transfer (s);
  }
}

// ---------------------------------------------------------------------------
/// Finishes all non-blocking exchanges of the charge density.
// ---------------------------------------------------------------------------
void
plasmaRho_finishExchange (meshDouble_p rho)
{
  socket_t *s;
  while ((s = socket_seekWait (&connectionRho.channel, mc_channel_income)))	// Receives and adds charge density.
  {
    double *pos = (double *) s->buffer;
    const regList_t *regs = (regList_t *) s->boss;
    reg_t *reg = regs->list;
    for (const reg_t * const end = regs->list + regs->N ; reg < end ; ++reg)
    {
//       SAY_DEBUG ("Getting rho from cpu %d: %s.", s->cpu, reg_printRanges(reg));
      for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
        for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
          for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos)
            mv_f(rho, i, j, k) += *pos;
    }
    socket_unlock (s);
  }

  syncMesh_waitSend (&connectionRho);						// Finishes all sending.
  assert (syncMesh_roundIsCompleted (&connectionRho));				// All data must be sended/received.

  plasmaRho_applyBoundary (rho);

  VSP_scalarSweep (rho);							// VSP convolution.
}

// ---------------------------------------------------------------------------
/// Accumulates charge density.
// ---------------------------------------------------------------------------
void
plasmaRho_add (meshDouble_p rho, marker_t *plasma, long int N)
{
  for (marker_t *end = plasma + N ; plasma < end ; ++plasma)
  {
    int i, j, k;
    double sigmaX, sigmaY, sigmaZ;

    MF_PIC_SIGMA (plasma->x, h1, i, sigmaX);
    MF_PIC_SIGMA (plasma->y, h2, j, sigmaY);
    MF_PIC_SIGMA (plasma->z, h3, k, sigmaZ);
    mf_PIC_distribute (rho, i, sigmaX, j, sigmaY, k, sigmaZ, plasma->rho);
  }
}

// ---------------------------------------------------------------------------
/// Initializes parallel exchange infrastructure to do parallel accumulation of charge density in ghost cells.
// ---------------------------------------------------------------------------
void
plasmaRho_setupConnections (void)
{
  // Domain of influence of particles.
  const reg_t reg = {.min = {cpu_min[0] - 2*mc_have_x, cpu_min[1] - 2*mc_have_y, cpu_min[2] - 2*mc_have_z},
                     .max = {cpu_max[0] + 2*mc_have_x, cpu_max[1] + 2*mc_have_y, cpu_max[2] + 2*mc_have_z}},
              ghostJ = {.min = {-1, -1, -1}, .max = {1, 1, 1}};			// Extends to grab all ghost cells.

  syncMesh_createConnection (&connectionRho, sizeof (double), &ghostJ);
  syncMesh_addReg (&connectionRho, &reg);
  syncMesh_syncronize (&connectionRho);						// Builds framework.
  syncMesh_verify (&connectionRho);						// Checks that all inputs are unique.
  syncMesh_dump (&connectionRho);
}
