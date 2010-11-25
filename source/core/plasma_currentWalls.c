/** \file plasma_currentWalls.c
  *
  * THIS MODULE IS INCLUDED INTO plasma.c AND INHERITS ALL DOMAIN SIZES FROM THERE.
  *
  * This module just applies boundary conditions for current density.
  *
  * NOTES:
  *   cbc prefixs stands for Currents Boundary Conditions.
  */

#ifndef MC_PLASMA_CURRENTWALLS_C
#define MC_PLASMA_CURRENTWALLS_C

#include <assert.h>

#include "parr_meshes.h"
#include "misc_definedKeys.h"

static void jbc_start (meshVec_p J);

// ---------------------------------------------------------------------------
/// \brief Updates current density in the boundary cells with respect to proper boundary condition for
/// particles (periodic or rigid wall). It also evaluates current density at i = mc_imax + 1, ..
/// layers to use in charge conservation routines without messing up that part.
// ---------------------------------------------------------------------------
static void
plasmaCurrents_boundaryConditions (void)
{
  const int sizes[3] = {cpu_max[0] - cpu_min[0], cpu_max[1] - cpu_min[1], cpu_max[2] - cpu_min[2]};	// Packed parameters.
  const int strides[3] = {mc_have_x*plasma_J->width_yz, mc_have_y*plasma_J->width_z, mc_have_z};
  const int BC[6] = {cpu_bc_min[0], cpu_bc_max[0], cpu_bc_min[1], cpu_bc_max[1], cpu_bc_min[2], cpu_bc_max[2]};

  static const int axis[6][3] = { {mc_x, mc_y, mc_z}, {mc_x, mc_y, mc_z}, {mc_y, mc_x, mc_z}, {mc_y, mc_x, mc_z}, {mc_z, mc_x, mc_y}, {mc_z, mc_x, mc_y} };
  static const int top[6] = {0, 1, 0, 1, 0, 1};

  meshVecI_p J = mcast_meshVecI (plasma_J);

  for (int b = 0 ; b < 6 ; ++b)							// Does all boundaries.
  {
    int pos0;
    const int ea = axis[b][0], ep = axis[b][1], eq = axis[b][2],
              dp = sizes[ep] + 4*ACTIVATOR[ep], dq = sizes[eq] + 4*ACTIVATOR[eq];
    const int stepOut = strides[ea]*(2*top[b] - 1);

    if (!ACTIVATOR[ea])								// Skips nonactive axises.
      continue;

    // Gets the offset of the point in the origin if the separation plane.
    pos0 = mf_offset(plasma_J, cpu_min[0], cpu_min[1], cpu_min[2]) - 2*ACTIVATOR[ep]*strides[ep] - 2*ACTIVATOR[eq]*strides[eq] + top[b]*strides[ea]*sizes[ea];

    switch (BC[b])								// Chooses BC method.
    {
      case BC_PERIODIC:
      {
        int shift = strides[ea]*sizes[ea], p;

        if (!top[b])								// Periodic done @ bottom boundary.
        {
          for (p = 0 ; p <= dp ; p++, pos0 += strides[ep])
          {
            int pos = pos0, q;
            for (q = 0 ; q <= dq ; q++, pos += strides[eq])
            {
              int pos2 = pos - 2*strides[ea], r;
              for (r = 0 ; r < 5 ; r++, pos2 += strides[ea])
              {
                mv_unrl (J, pos2 + shift).r[0] += mv_unrl (J, pos2).r[0];	// All components are added, so
                mv_unrl (J, pos2 + shift).r[1] += mv_unrl (J, pos2).r[1];	//   {ea, ep, eq} written as {0, 1, 2}.
                mv_unrl (J, pos2 + shift).r[2] += mv_unrl (J, pos2).r[2];
                mv_unrl (J, pos2) = mv_unrl (J, pos2 + shift);
              }
            }
          }
        }
      }
      break;

      case BC_MIRROR:
      case BC_OPEN:
      {
        int internal = (1 - top[b])*strides[ea];	// J_axis is shifted by half of the node so to access it I need to remember
        int external = top[b]*strides[ea];		// that sign of shift is different for top and bottom cases.

        for (int p = 0 ; p <= dp ; ++p, pos0 += strides[ep])
        {
          int pos = pos0, q;
          for (q = 0 ; q <= dq ; ++q, pos += strides[eq])
          {
            mv_unrl (J, pos + internal).r[ea] -= mv_unrl (J, pos + external).r[ea];	// E_normal -= E_normal'
            mv_unrl (J, pos - stepOut).r[ep] += mv_unrl (J, pos + stepOut).r[ep];	// E_tangent += E_tangent'
            mv_unrl (J, pos - stepOut).r[eq] += mv_unrl (J, pos + stepOut).r[eq];	// E_tangent += E_tangent'

            // VSP extention of the same stuff.
            mv_unrl (J, pos + internal - stepOut).r[ea] -= mv_unrl (J, pos + external + stepOut).r[ea];
            mv_unrl (J, pos - 2*stepOut).r[ep] += mv_unrl (J, pos + 2*stepOut).r[ep];
            mv_unrl (J, pos - 2*stepOut).r[eq] += mv_unrl (J, pos + 2*stepOut).r[eq];

            // Forces full charge density of the mirrored charges.
            mv_unrl (J, pos).r[ep] = mv_unrl (J, pos).r[eq] = 0;

            mv_unrl (J, pos + external).r[ea] = mv_unrl (J, pos + internal).r[ea];	// E_normal -= E_normal'
            mv_unrl (J, pos + stepOut).r[ep] = - mv_unrl (J, pos - stepOut).r[ep];	// E_tangent += E_tangent'
            mv_unrl (J, pos + stepOut).r[eq] = - mv_unrl (J, pos - stepOut).r[eq];	// E_tangent += E_tangent'

            mv_unrl (J, pos + external + stepOut).r[ea] = mv_unrl (J, pos + internal - stepOut).r[ea];
            mv_unrl (J, pos + 2*stepOut).r[ep] = - mv_unrl (J, pos - 2*stepOut).r[ep];
            mv_unrl (J, pos + 2*stepOut).r[eq] = - mv_unrl (J, pos - 2*stepOut).r[eq];
          }
        }
      }
      break;

      default:
        DIE ("unsupported boundary condition");
      break;

      case BC_SPLITTER:
      break;
    }
  }
}

// ---------------------------------------------------------------------------
/// Initializes all connection infrastructure to do parallel accumulation of the \f$ \vec j,\ \rho \f$ in the ghost cells.
// ---------------------------------------------------------------------------
static void
plasma_setupConnections (void)
{
  const reg_t ghostJ = {{-1, -1, -1}, {1, 1, 1}, 0, 0, {0, 0, 0}};		// Extends to grab all ghost cells.
  const reg_t reg = {{cpu_min[0] - 2*mc_have_x, cpu_min[1] - 2*mc_have_y, cpu_min[2] - 2*mc_have_z},// Domain of influence of particles.
                     {cpu_max[0] + 2*mc_have_x, cpu_max[1] + 2*mc_have_y, cpu_max[2] + 2*mc_have_z}};

  syncMesh_createConnection (&connectionJ, sizeof (vec3D_t), &ghostJ);
  syncMesh_addReg (&connectionJ, &reg);
  syncMesh_syncronize (&connectionJ);						// Builds frameworks.
  syncMesh_verify (&connectionJ);						// Checks that all inputs are unique.
  syncMesh_dump (&connectionJ);
}

// ---------------------------------------------------------------------------
/// Starts all non-blocking exchanges.
/// \todo Post syncMesh_irecv() in the beginning of the plasma_step or even of
///       the main simulation loop.
// ---------------------------------------------------------------------------
static void
jbc_start (meshVec_p J)
{
  syncMesh_irecv (&connectionJ);						// Starts non-blocking receive.

  int N;
  socket_t *s;
  mf_syncMesh_put(connectionJ, s, N);
  for (const socket_t * const end = s + N ; s < end ; ++s)			// Packs and sends H-data.
  {
    const regList_t *list = (regList_t*) s->boss;
    reg_t *reg = list->list;
    reg_t * const regEnd = list->list + list->N;
    vec3D_t *pos = (vec3D_t *) s->buffer;
    for ( ; reg < regEnd ; ++reg)
      for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
        for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
          for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos)
            *pos = mv_f(J, i, j, k);
    socket_transfer (s);
  }
}

// ---------------------------------------------------------------------------
/// Finishes all non-blocking exchanges.
// ---------------------------------------------------------------------------
void
jbc_finish (meshVec_p J)
{
  profiler_begin (mc_prof_plasma_jbc_mesh);
  socket_t *s;
  while ((s = socket_seekWait (&connectionJ.channel, mc_channel_income)))	// Receives and adds currents.
  {
    const regList_t *regs = (regList_t *) s->boss;
    vec3D_t *pos = (vec3D_t *) s->buffer;
    reg_t *reg = regs->list;
    const reg_t *end = regs->list + regs->N;
    for ( ; reg < end ; ++reg)
    {
//       SAY_DEBUG ("Getting J-field contribution from cpu %d: %s.", reg->cpu, reg_printRanges(reg));
      for (int i = reg->min[0] ; i <= reg->max[0] ; ++i)
        for (int j = reg->min[1] ; j <= reg->max[1] ; ++j)
          for (int k = reg->min[2] ; k <= reg->max[2] ; ++k, ++pos)
          {
            mv_fx(J, i, j, k) += pos->x;
            mv_fy(J, i, j, k) += pos->y;
            mv_fz(J, i, j, k) += pos->z;
          }
    }
    socket_unlock (s);
  }

  syncMesh_waitSend (&connectionJ);						// Finishes all sending.
  assert (syncMesh_roundIsCompleted (&connectionJ));				// All must be sended / received at this point.

  profiler_endBegin (mc_prof_plasma_jbc_localBC);
  plasmaCurrents_boundaryConditions ();						// Boundary condition for current density.

  profiler_endBegin (mc_prof_plasma_jbc_sweep);
  VSP_vectorSweep (mcast_meshVec (plasma_J));					// Convolution PIC->VSP.

  profiler_endBegin (mc_prof_plasma_jbc_particles);
  // XXX pbc_finalize ();							// Finalizes || exchange of particles.
  profiler_end ();
}

#endif
