/** \file emCap_periodic.c
  * Electromagnetic module: periodic boundary condition functions used in the em_cap.c.
  */

#ifndef mc_emCap_periodic_includer
#define mc_emCap_periodic_includer						///< Multiple including guard.

// ---------------------------------------------------------------------------
/// Applies \b periodic BC for magnetic field.
// ---------------------------------------------------------------------------
static void
capBC_periodic_H (meshVec_p minH, meshVec_p maxH, int ea)
{
   int ep = capAxisFrameFinder[ea+1],
       eq = capAxisFrameFinder[ea+2];
   int strideBtm[3] = {mc_have_x*minH->width_yz, mc_have_y*minH->width_z, mc_have_z};					/// \todo Strides in unrolled array.
   int strideTop[3] = {mc_have_x*maxH->width_yz, mc_have_y*maxH->width_z, mc_have_z};				/// \todo Keep in reg.wrap field of the mesh.
   int min[3] = {minH->imin, minH->jmin, minH->kmin},
       max[3] = {minH->imax, minH->jmax, minH->kmax};	/// \todo Replace by pointers to reg_t minH->min.
   int dp     = max[ep] - min[ep],
       dq     = max[eq] - min[eq];								/// \todo Move to the cpu.wrap/minH.wrap array.
/*
  say ("Applying periodic BC for mesh '%s/%s' along axis %d.", minH->name, maxH->name, ea);*/

   // Inits cursors.
   int posBtm0 = strideBtm[ea]*cpu_min[ea] + strideBtm[ep]*min[ep] + strideBtm[eq]*min[eq];
   int posTop0 = strideTop[ea]*cpu_max[ea] + strideTop[ep]*min[ep] + strideTop[eq]*min[eq];

   for (int p = 0 ; p <= dp ; ++p, posBtm0 += strideBtm[ep], posTop0 += strideTop[ep]) {
      int posTop = posTop0,
          posBtm = posBtm0;
      for (int q = 0 ; q <= dq ; ++q, posTop += strideTop[eq], posBtm += strideBtm[eq]) {
         mv_unrl (minH, posBtm) = mv_unrl (maxH, posTop);									// Sets values on bottom boundary.
         mv_unrl (maxH, posTop + strideTop[ea]) = mv_unrl (minH, posBtm + strideBtm[ea]);					// Sets values on top boundary.
      }
   }
}

// ---------------------------------------------------------------------------
/// Applies \b periodic BC for electric field.
// ---------------------------------------------------------------------------
static void
capBC_periodic_E (meshVec_p minE, meshVec_p maxE, int ea)
{
   int ep = capAxisFrameFinder[ea+1],
       eq = capAxisFrameFinder[ea+2];
   int strideBtm[3] = {mc_have_x*minE->width_yz, mc_have_y*minE->width_z, mc_have_z};				/// \todo Strides in unrolled array.
   int strideTop[3] = {mc_have_x*maxE->width_yz, mc_have_y*maxE->width_z, mc_have_z};				/// \todo Keep in reg.wrap field of the mesh.
   int min[3] = {minE->imin, minE->jmin, minE->kmin},
       max[3] = {minE->imax, minE->jmax, minE->kmax};		/// \todo Replace by pointers to reg_t minE->min.
   int dp     = max[ep] - min[ep],
       dq     = max[eq] - min[eq];								/// \todo Move to the cpu.wrap/minH.wrap array.
/*
  SAY_DEBUG ("Applying periodic BC for mesh '%s/%s' along axis %d.", minE->name, maxE->name, ea);*/

  // Inits cursors.
   int posBtm0 = strideBtm[ea]*cpu_min[ea] + strideBtm[ep]*min[ep] + strideBtm[eq]*min[eq];
   int posTop0 = strideTop[ea]*cpu_max[ea] + strideTop[ep]*min[ep] + strideTop[eq]*min[eq];

   for (int p = 0 ; p <= dp ; ++p, posBtm0 += strideBtm[ep], posTop0 += strideTop[ep]) {
      int posTop = posTop0,
          posBtm = posBtm0;
      for (int q = 0 ; q <= dq ; ++q, posTop += strideTop[eq], posBtm += strideBtm[eq]) {
         mv_unrl (maxE, posTop) = mv_unrl (minE, posBtm);									// Sets values on bottom boundary.
         mv_unrl (maxE, posTop + strideTop[ea]) = mv_unrl (minE, posBtm + strideBtm[ea]);					// Sets values on top boundary.
      }
   }
}

#endif
