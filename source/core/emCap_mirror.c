/** \file emCap_mirror.c
  * Electromagnetic module: mirror boundary condition functions used in the em_cap.c.
  *
  * This boundary condition uses (anti) symmetry of the solution with respect to the boundary.
  * Symmetric components are \f$ E_n, \vec H_\tau \f$, antisymmetric components are \f$ H_n, \vec E_\tau \f$ ( \b n means component normal
  * to boundary and \f$ \tau \f$ - tangential component). Antisymmetric components on the boundary itself are equal to zero. Please keep in
  * mind that default shift is -1/2 (see em.h).
  */

#ifndef mc_emCap_mirror_includer
#define mc_emCap_mirror_includer						///< Multiple including guard.

// ---------------------------------------------------------------------------
/// Applies \b mirror BC for magnetic field on the boundary pointed by \b ea and \b top.
// ---------------------------------------------------------------------------
static void
capBC_mirror_H (meshVecI_p H, int ea, int top)
{
  const int ep = capAxisFrameFinder[ea+1], eq = capAxisFrameFinder[ea+2];
  const int stride[3] = {mc_have_x*H->width_yz, mc_have_y*H->width_z, mc_have_z};					/// \todo Strides in unrolled array.
  const int plane[2][3] = {{cpu_min[0], cpu_min[1], cpu_min[2]}, {cpu_max[0], cpu_max[1], cpu_max[2]}};				/// \todo Replace by 2 pointers to cpu.min/max.
  const int min[3] = {H->imin, H->jmin, H->kmin}, max[3] = {H->imax, H->jmax, H->kmax};					/// \todo Replace by pointers to reg_t H->min.
  const int dp = max[ep] - min[ep], dq = max[eq] - min[eq];								/// \todo Move to the cpu.wrap/minH.wrap array.
  const int intl = (1 - top)*stride[ea], extnl = top*stride[ea];							// Shifts to different sides of the boundary.
/*
  say ("Applying mirror BC for mesh '%s' along axis %d / %d.", H->name, ea, top);*/

  int pos0 = stride[ea]*plane[top][ea] + stride[ep]*min[ep] + stride[eq]*min[eq];					// Inits cursor.

  for (int p = 0 ; p <= dp ; ++p, pos0 += stride[ep])
  {
    int pos = pos0;
    for (int q = 0 ; q <= dq ; ++q, pos += stride[eq])
    {
      mv_unrl (H, pos).r[ea] = 0;											// Erases H_normal.
      mv_unrl (H, pos + extnl).r[ep] = mv_unrl (H, pos + intl).r[ep];							// Copies H_tangent.
      mv_unrl (H, pos + extnl).r[eq] = mv_unrl (H, pos + intl).r[eq];							// Copies H_tangent.
    }
  }
}

// ---------------------------------------------------------------------------
/// Applies \b mirror BC for electric field on the boundary pointed by \b ea and \b top.
// ---------------------------------------------------------------------------
static void
capBC_mirror_E (meshVecI_p E, int ea, int top)
{
   const int ep = capAxisFrameFinder[ea+1],
             eq = capAxisFrameFinder[ea+2];
   const int stride[3] = {mc_have_x*E->width_yz, mc_have_y*E->width_z, mc_have_z};					/// \todo Strides in unrolled array.
   const int plane[2][3] = {{cpu_min[0], cpu_min[1], cpu_min[2]},
                            {cpu_max[0], cpu_max[1], cpu_max[2]}};				/// \todo Replace by 2 pointers to cpu.min/max.
   const int min[3] = {E->imin, E->jmin, E->kmin},
             max[3] = {E->imax, E->jmax, E->kmax};					/// \todo Replace by pointers to reg_t E->min.
   const int dp = max[ep] - min[ep],
             dq = max[eq] - min[eq];								/// \todo Move to the cpu.wrap/minE.wrap array.
   const int intl  = (1 - top)*stride[ea],
             extnl = top*stride[ea];							// Shifts to different sides of the boundary.
   /*
   say ("Applying mirror BC for mesh '%s' along axis %d / %d.", E->name, ea, top);*/

   int pos0 = stride[ea]*plane[top][ea] + stride[ep]*min[ep] + stride[eq]*min[eq];					// Inits cursor.

   for (int p = 0 ; p <= dp ; ++p, pos0 += stride[ep]) {
      int pos = pos0;
      for (int q = 0 ; q <= dq ; ++q, pos += stride[eq]) {
         mv_unrl (E, pos).r[ep] = 0;	// Erases tangent component of E.
         mv_unrl (E, pos).r[eq] = 0;
         mv_unrl (E, pos + extnl).r[ea] = mv_unrl (E, pos + intl).r[ea];							// Copies E_normal.
      }
   }
}

#endif
