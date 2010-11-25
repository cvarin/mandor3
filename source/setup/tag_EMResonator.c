/** \file tag_EMResonator.c
  * Interface to the standing wave test.
  *
  * Excites standing EM wave in a metal box. Polarisation: nonzero components are Ex, Hy, Hz.
  */

#include <math.h>
#include <stdlib.h>

#include "type_mesh.h"

#include "mf_vector.h"

#include "main.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"

static int    m[3]  = {0, 0, 0};	///< Wave numbers.
static double k[3]  = {0, 0, 0};	///< Wave vector (can be initialized numerically / physically, see tag_EMWave.h).
static double omega = 0;		///< Cyclic frequency (see tag_EMWave.h).
static double Dt    = 0;		///< Cyclic frequency modified by dispersion of Yee scheme (see tag_EMWave.h).
static double Dr[3] = {0, 0, 0};	///< Wave vector modified by dispersion of Yee scheme (see tag_EMWave.h).

// ----------------------------------------------------------------------------
/// Initializes ω and k using exact physical dispersion equation and sets
/// Dr, Dt to physical values.
// ----------------------------------------------------------------------------
static void
EMResonator_setPhysFrame (void)
{
   k[0] = mc_have_x*mc_pi*m[0]/(h1*(1 - mc_have_x + dmn_max[0] - dmn_min[0]));	// Turns wave numbers to wave vector.
   k[1] = mc_have_y*mc_pi*m[1]/(h2*(1 - mc_have_y + dmn_max[1] - dmn_min[1]));
   k[2] = mc_have_z*mc_pi*m[2]/(h3*(1 - mc_have_z + dmn_max[2] - dmn_min[2]));

   omega = Dt = sqrt (mf_dotProduct (k, k));										// Gets omega from normal dispersion equation.
   mf_vectorCopy (k, Dr);												// Updates numerical wave vector.
}

// ----------------------------------------------------------------------------
/// Initializes ω and k using exact numerical dispersion equation, sets Dr, Dt.
// ----------------------------------------------------------------------------
static void
EMResonator_setYee2Frame (void)
{
   EMResonator_setPhysFrame ();
   mf_vectorSet(2*sin (0.5*k[0]*h1)/h1,
                2*sin (0.5*k[1]*h2)/h2,
                2*sin (0.5*k[2]*h3)/h3, Dr);
   Dt    = sqrt (mf_dotProduct (Dr, Dr));
   omega = asin (0.5*tau*Dt)*2/tau;
}

// ----------------------------------------------------------------------------
/// Interface to the tag initializator.
// ----------------------------------------------------------------------------
void
tag_EMResonator (FILE *fp, meshVec_p E)
{
   const char *typeWord = cfg_readWord (fp);											// Reads type of the initialization.
   int type = cfg_identifyWord (typeWord, "physical", 0, "Yee_2nd_order", 1, mc_cfgTermGuesses);
   say ("tag_EMResonator:");
   say ("  - %s dispersion equation used,", typeWord);

   m[0] = cfg_readInt (fp)*mc_have_x;											// Reads wave numbers.
   m[1] = cfg_readInt (fp)*mc_have_y;
   m[2] = cfg_readInt (fp)*mc_have_z;

   typeWord = cfg_readWord (fp);
   const int ea = cfg_identifyWord (typeWord, "Ex", 0, "Ey", 1, "Ez", 2, mc_cfgTermGuesses);
   ENSURE (ea != -1, "unknown field component '%s'", typeWord);
   say ("  - polarization '%s'", typeWord);

   m[ea] = 0;														// Accounts for polarization.
   ENSURE (mf_dotProduct(m, m),
           "bad wave numbers (%le, %le, %le)", m[0], m[1], m[2]);

   // Solves dispertion equation for omega.
   switch (type) {
   case 0:
      EMResonator_setPhysFrame ();
      break;

   case 1:
      EMResonator_setYee2Frame ();
      break;

   default:
      DIE ("unknown initialization type '%s'", typeWord);
      break;
   }
   say ("  - k  = (%e, %e, %e), omega = %e,", k[0], k[1], k[2], omega);
   say ("  - Dr = (%e, %e, %e), Dt    = %e,", Dr[0], Dr[1], Dr[2], Dt);

   double E0 = cfg_readDouble (fp);											// Reads field amplitude.
   say ("  - amplitude %e", E0);

   if (memEstimateOnly)
      return;

   k[0] *= h1;														// Includes mesh steps into k.
   k[1] *= h2;
   k[2] *= h3;

   const int ep = (ea + 1)%3,
             eq = (ea + 2)%3;
   int pos[3];
   for (pos[0] = E->imin ; pos[0] <= E->imax ; ++pos[0])
   for (pos[1] = E->jmin ; pos[1] <= E->jmax ; ++pos[1])
   for (pos[2] = E->kmin ; pos[2] <= E->kmax ; ++pos[2]) {
      mv_fi(mcast_meshVecI(E), pos[0], pos[1], pos[2], ea) +=
                                                      E0*sin (k[ep]*pos[ep])
                                                        *sin (k[eq]*pos[eq]);
   }
}
