/** \file emCap_Mur.c
  * Electromagnetic module: first order Mur's absorbing boundary condition.
  *
  * <b>First order Mur's absorbing boundary condition (short derivation).</b>
  *
  * -# LHS of the 1D Maxwell equation is written in form of superposition of two
  *    differential cancellation operators (for left/right propagating waves):
  *    \f[ (c^2\cdot \partial^2_{zz} - \partial^2_{tt}) E =
  *         (c\cdot \partial_{z} - \partial_{t})
  *         (c\cdot \partial_{z} + \partial_{t}) E \f]
  *
  * -# Lets consider left wall, where equation below should be satisfied to
  *    ensure full absorption (no wave is propagating to the right):
  *    \f[ (c\cdot \partial_{z} - \partial_{t}) E = 0 \f]
  *    Lets use index \b 0 for boundary and \b 1 for the next node inside of the
  *    domain (in other words, \b 1 points to the closest to the boundary
  *    \b inner node of the domain).
  *
  * -# Using mesh node placement from Yee solver (see em.h) this equation may be
  *    approximated in finite difference form using central interpolation:
  *    \f[ \frac c2\cdot \left(\frac {E^{n+1}_1 - E^{n+1}_0}{h} +
  *                            \frac {E^{n}_1 - E^{n}_0}{h}\right) -
  *        \frac 12\cdot \left(\frac {E^{n+1}_1 - E^{n}_1}{\tau} +
  *                            \frac {E^{n+1}_0 - E^{n}_0}{\tau}\right) = 0 \f]
  *
  * -# Straighforward math delivers final expression:
  *    \f[ E^{n+1}_0  = E^{n}_1 + \frac {c\tau - h}{c\tau + h}\cdot (E^{n+1}_1 - E^{n}_0) \f]
  *
  * This expression is the boundary condition used in code.
  *
  * Full theory is described in <b>Gerrit Mur, "Absorbing boundary condition for
  * the finite-difference approximation of the time-domain electromagnetic-field
  * equation", IEEE transaction on electromagnetic compability, 1981,
  * vol.EMC-23, No.4, November, pp.377-382.</b> or see a book of Taflove
  * (referenced in em_TFSF.h).
  *
  * \note In boundary condition old value on the boundary is used. This may
  * create some problems if node is under the influence of few boundaries. To
  * avoid any conflicts of this type old values are stored in the buffers
  * (allocated here and released by caller). This way all updates are
  * order-independent.
  *
  * \note Some component of the field are copied to external ghost cells so
  * particles may safely interpolate field in the vicinity of the boundary.
  */

#ifndef MC_EMCAP_MUR_HEADER
#define MC_EMCAP_MUR_HEADER

// ---------------------------------------------------------------------------
/// Applies \b Mur BC for magnetic field on the boundary pointed by \b ea and \b top.
/// \todo Add one of the shifts to the \b pos: "mv_unrl (H, pos + extnl).r[ep] = mv_unrl (H, pos + intl).r[ep]" ->
/// "mv_unrl (H, pos + shiftOut).r[ep] = mv_unrl (H, pos).r[ep]"
// ---------------------------------------------------------------------------
static void
capBC_Mur_H (meshVecI_p H, int ea, int top)
{
  const int ep = capAxisFrameFinder[ea+1], eq = capAxisFrameFinder[ea+2];
  const int stride[3] = {mc_have_x*H->width_yz, mc_have_y*H->width_z, mc_have_z};					/// \todo Strides in unrolled array.
  const int plane[2][3] = {{cpu_min[0], cpu_min[1], cpu_min[2]}, {cpu_max[0], cpu_max[1], cpu_max[2]}};				/// \todo Replace by 2 pointers to cpu.min/max.
  const int min[3] = {H->imin, H->jmin, H->kmin}, max[3] = {H->imax, H->jmax, H->kmax};					/// \todo Replace by pointers to reg_t H->min.
  const int dp = max[ep] - min[ep], dq = max[eq] - min[eq];								/// \todo Move to the cpu.wrap/minH.wrap array.
  const int intl = (1 - top)*stride[ea], extnl = top*stride[ea];							// Shifts to different sides of the boundary.
/*
  say ("Applying Mur BC for mesh '%s' along axis %d / %d.", H->name, ea, top);*/

  int pos0 = stride[ea]*plane[top][ea] + stride[ep]*min[ep] + stride[eq]*min[eq];					// Inits cursor.

  for (int p = 0 ; p <= dp ; ++p, pos0 += stride[ep])
  {
    int pos = pos0;
    for (int q = 0 ; q <= dq ; ++q, pos += stride[eq])
    {
      mv_unrl (H, pos + extnl).r[ep] = mv_unrl (H, pos + intl).r[ep];							// Copies H_tangent outside.
      mv_unrl (H, pos + extnl).r[eq] = mv_unrl (H, pos + intl).r[eq];							// Copies H_tangent outside.
    }
  }
}

/*
 * Unrolled offsets in the storage array.
 */
#define mc_Ep0		0	///< Ep on boundary.
#define mc_Eq0		1	///< Eq on boundary.
#define mc_Ep1		2	///< Ep on the first innermost point.
#define mc_Eq1		3	///< Eq on the first innermost point.
#define mc_chunk	4	///< Size of the total record.

// ---------------------------------------------------------------------------
/// \brief This function is used to prepare Mur boundary condition buffers and reset external field. Allocated memory should be released by caller.
// ---------------------------------------------------------------------------
static void
capBC_Mur_init (meshVecI_p E, meshVec_p capE, int ea, int top, capMur_t *pack, double tau, double h)
{
  const int ep = capAxisFrameFinder[ea+1], eq = capAxisFrameFinder[ea+2];
  const int stride[3] = {mc_have_x*E->width_yz, mc_have_y*E->width_z, mc_have_z};					/// \todo Strides in unrolled array.
  const int plane[2][3] = {{cpu_min[0], cpu_min[1], cpu_min[2]}, {cpu_max[0], cpu_max[1], cpu_max[2]}};				/// \todo Replace by 2 pointers to cpu.min/max.
  const int min[3] = {capE->imin, capE->jmin, capE->kmin}, max[3] = {capE->imax, capE->jmax, capE->kmax};		/// \todo Replace by pointers to reg_t E->min.
  const int dp = max[ep] - min[ep], dq = max[eq] - min[eq];								/// \todo Move to the cpu.wrap/minE.wrap array.
  const int intl = (1 - top)*stride[ea], extnl = top*stride[ea], shiftIn = (1 - 2*top)*stride[ea];			// Shifts to different sides of the boundary.
/*
  say ("Preparing Mur BC storages for mesh '%s' along axis %d / %d.", E->name, ea, top);*/

  ENSURE (!pack->buffer, "double allocation of memory");

  pack->buffer = (double*) malloc (mc_chunk*(dp + 1)*(dq + 1)*sizeof (double));
  ENSURE (pack->buffer, "out of RAM (%.2f Kb)", dp*dq*sizeof (double)/1024.0);

  pack->alpha = (tau - h)/(tau + h);											// Inits packed coefficients.
  pack->ea = ea;
  pack->top = top;
  pack->E = mcast_meshVecI (capE);

  double *buffer = pack->buffer;

  int pos0 = stride[ea]*plane[top][ea] + stride[ep]*min[ep] + stride[eq]*min[eq];					// Inits cursor.
  int old = 0;

  for (int p = 0 ; p <= dp ; ++p, pos0 += stride[ep])
  {
    int pos = pos0;
    for (int q = 0 ; q <= dq ; ++q, pos += stride[eq], old += mc_chunk)
    {
      buffer[old+mc_Ep0] = mv_unrl (E, pos).r[ep];									// Updates buffered old E values.
      buffer[old+mc_Eq0] = mv_unrl (E, pos).r[eq];
      buffer[old+mc_Ep1] = mv_unrl (E, pos + shiftIn).r[ep];
      buffer[old+mc_Eq1] = mv_unrl (E, pos + shiftIn).r[eq];

      mv_unrl (E, pos + extnl).r[ea] = mv_unrl (E, pos + intl).r[ea];							// Copies E_normal.
    }
  }
}

// ---------------------------------------------------------------------------
/// Applies \b Mur BC for electric field on the boundary pointed by \b ea and \b top.
// ---------------------------------------------------------------------------
static void
capBC_Mur_E (capMur_t *pack)
{
  meshVecI_p E = pack->E;												// Uses packed parameters passed by pointer.
  const int ea = pack->ea, top = pack->top;

  const int ep = capAxisFrameFinder[ea+1], eq = capAxisFrameFinder[ea+2];						// Usual construction of frame.
  const int stride[3] = {mc_have_x*E->width_yz, mc_have_y*E->width_z, mc_have_z};					/// \todo Strides in unrolled array.
  const int plane[2][3] = {{cpu_min[0], cpu_min[1], cpu_min[2]}, {cpu_max[0], cpu_max[1], cpu_max[2]}};				/// \todo Replace by 2 pointers to cpu.min/max.
  const int min[3] = {E->imin, E->jmin, E->kmin}, max[3] = {E->imax, E->jmax, E->kmax};					/// \todo Replace by pointers to reg_t E->min.
  const int dp = max[ep] - min[ep], dq = max[eq] - min[eq];								/// \todo Move to the cpu.wrap/minE.wrap array.
  const int intl = (1 - top)*stride[ea], extnl = top*stride[ea], shiftIn = (1 - 2*top)*stride[ea];			// Shifts to different sides of the boundary.
/*
  say ("Applying Mur BC for mesh '%s' along axis %d / %d.", E->name, ea, top);*/

  int pos0 = stride[ea]*plane[top][ea] + stride[ep]*min[ep] + stride[eq]*min[eq];					// Inits cursor.
  int old = 0;
  double * const buffer = pack->buffer;											// Local copies to avoid indirect references.
  const double alpha = pack->alpha;

  for (int p = 0 ; p <= dp ; ++p, pos0 += stride[ep])
  {
    int pos = pos0;
    for (int q = 0 ; q <= dq ; ++q, pos += stride[eq], old += mc_chunk)
    {
      // Copy for particles.
      mv_unrl (E, pos + extnl).r[ea] = mv_unrl (E, pos + intl).r[ea];

      // BC for tangent component of E.
      mv_unrl (E, pos).r[ep] = buffer[old+mc_Ep1] + alpha*(mv_unrl (E, pos + shiftIn).r[ep] - buffer[old+mc_Ep0]);
      mv_unrl (E, pos).r[eq] = buffer[old+mc_Eq1] + alpha*(mv_unrl (E, pos + shiftIn).r[eq] - buffer[old+mc_Eq0]);

      buffer[old+mc_Ep0] = mv_unrl (E, pos).r[ep];									// Updates buffered old E values.
      buffer[old+mc_Eq0] = mv_unrl (E, pos).r[eq];
      buffer[old+mc_Ep1] = mv_unrl (E, pos + shiftIn).r[ep];
      buffer[old+mc_Eq1] = mv_unrl (E, pos + shiftIn).r[eq];
    }
  }
}

#undef mc_Ep0
#undef mc_Eq0
#undef mc_Ep1
#undef mc_Eq1
#undef mc_chunk

#endif
