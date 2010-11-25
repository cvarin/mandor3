/** \file tag_EMWave.c
  * \brief Excites uniform plane EM-wave in the domain (theory is in
  *        'tag_EMWave.h').
  */

#include <math.h>
#include <stdlib.h>

#include "type_mesh.h"

#include "mf_vector.h"

#include "log.h"
#include "main.h"
#include "misc_units.h"
#include "misc_cfgReader.h"

static int    m[3]  = {0, 0, 0};	///< Wave numbers.
static double k[3]  = {0, 0, 0};	///< Wave vector (can be initialized
                                        ///<   numerically / physically, see
                                        ///<   'tag_EMWave.h').
static double omega = 0;		///< Cyclic frequency.
static double Dt    = 0;		///< Cyclic frequency modified by dispersion of Yee scheme (see tag_EMWave.h).
static double Dr[3] = {0, 0, 0};	///< Wave vector modified by dispersion of Yee scheme (see tag_EMWave.h).

// ---------------------------------------------------------------------------
/// Initializes \f$ \omega \f$ and \f$ \vec k \f$ using exact physical
/// dispersion equation and sets \f$ \vec D_r, D_t \f$ to physical values.
// ---------------------------------------------------------------------------
static void
EMWave_setPhysFrame (void)
{
   // Turns wave numbers to wave vector.
   k[0] = mc_have_x*2*mc_pi*m[0]/(h1*(1 - mc_have_x + dmn_max[0] - dmn_min[0]));
   k[1] = mc_have_y*2*mc_pi*m[1]/(h2*(1 - mc_have_y + dmn_max[1] - dmn_min[1]));
   k[2] = mc_have_z*2*mc_pi*m[2]/(h3*(1 - mc_have_z + dmn_max[2] - dmn_min[2]));

   // Gets omega from normal dispersion equation, updates numerical wave vector.
   omega = Dt = sqrt (mf_dotProduct (k, k));
   mf_vectorCopy (k, Dr);
}

// ---------------------------------------------------------------------------
/// Initializes \f$ \omega \f$ and \f$ \vec k \f$ using exact numerical
/// dispersion equation and updates \f$ \vec D_r, D_t \f$.
// ---------------------------------------------------------------------------
static void
EMWave_setYee2Frame (void)
{
   EMWave_setPhysFrame ();			// Sets k and c vectors.
   mf_vectorSet(2*sin (0.5*k[0]*h1)/h1,		// Applies dispersion factors
                2*sin (0.5*k[1]*h2)/h2,		//   to wave vector.
                       2*sin (0.5*k[2]*h3)/h3, Dr);
   Dt    = sqrt (mf_dotProduct (Dr, Dr));	// Applies dispersion factors to frequency.
   omega = asin (0.5*tau*Dt)*2/tau;		// Gets real omega.
}

// ---------------------------------------------------------------------------
/// Excites uniform plane EM-wave (for periodic boundary conditions only).
// ---------------------------------------------------------------------------
void
tag_EMWave (FILE *fp, meshVec_p E, meshVec_p H)
{
   const char *typeWord = cfg_readWord (fp);
   int         type     = cfg_identifyWord (typeWord, "physical", 0,
                                                      "Yee_2nd_order", 1,
                                                      mc_cfgTermGuesses);
   m[0] = cfg_readInt (fp)*mc_have_x;
   m[1] = cfg_readInt (fp)*mc_have_y;
   m[2] = cfg_readInt (fp)*mc_have_z;
   ENSURE (mf_dotProduct(m, m),
           "zero wave numbers (%d, %d, %d)\n", m[0], m[1], m[2]);

   // Solves dispertion equation for omega.
   switch (type) {
   case 0:
      EMWave_setPhysFrame ();
      break;

   case 1:
      EMWave_setYee2Frame ();
      break;

   default:
      DIE ("unknown initialization type '%s'", typeWord);
      break;
   }

   double E0 = cfg_readDouble (fp),
          e[3];
   e[0] = cfg_readDouble (fp);
   e[1] = cfg_readDouble (fp);
   e[2] = cfg_readDouble (fp);

   // Corrects direction and amplitude of the polarization vector:
   //   e = e - k*(e,k)/(k,k),
   //   e = e*E0/|e|.
   double ek  = mf_dotProduct (Dr, e);
   double Dr2 = mf_dotProduct (Dr, Dr),
          tmp[3];
   mf_vectorCopy     (Dr,     tmp);
   mf_vectorScale    (ek/Dr2, tmp, tmp);
   mf_vectorSubtract (e,      tmp, e);

   Dr2 = sqrt (mf_dotProduct (e, e));											// Normalizes length of 'e'.
   ENSURE (Dr2 >= 0.1,
           "polarisation vector is too small: n_perp = %le\n", Dr2);
   mf_vectorScale (E0/Dr2, e, e);											// Sets amplitude of e to E.

  // Sets polarization vector for H: h = - [e, Dr]/Dt.
  double h[3];
  mf_crossProduct (Dr,     e, h);												// Gets direction of H.
  mf_vectorScale  (1.0/Dt, h, h);											// Gets direction of H.

  say ("tag_EMWave: ");
  say ("  - electromagnetic wave is excited using %s dispertion equation,",
       (type) ? "numerical" : "physical");
  say ("  - wave numbers are (%d, %d, %d),", m[0], m[1], m[2]);
  say ("  - E = (%e, %e, %e), E0 = %e", e[0], e[1], e[2], sqrt (mf_dotProduct (e, e)));
  say ("  - H = (%e, %e, %e), H0 = %e", h[0], h[1], h[2], sqrt (mf_dotProduct (h, h)));
  say ("  - k  = (%e, %e, %e), w  = %e", k[0], k[1], k[2], omega);
  say ("  - Dr = (%e, %e, %e), Dt = %e", Dr[0], Dr[1], Dr[2], Dt);

  typeWord = cfg_readWord (fp);												// Reads type of the initialization.
  type     = cfg_identifyWord (typeWord, "right", 0,
                                         "plane", 1,
                                         "left",  2,
                                         mc_cfgTermGuesses);
   ENSURE (type != -1, "unknown polarization '%s'", typeWord);

   --type;
   say ("  - polarization of the wave is '%s'", typeWord);

   if (dmn_bc_min[0] != BC_PERIODIC
   ||  dmn_bc_min[1] != BC_PERIODIC
   ||  dmn_bc_min[2] != BC_PERIODIC) {
      SAY_WARNING ("boundary conditions are not periodic");
   }

   if (memEstimateOnly)
      return;

   // Includes boundaries to ensure that initialization is done everywhere.
   if (mf_mesh_pointIsOutside (H, E->imin, E->jmin, E->kmin)
   ||  mf_mesh_pointIsOutside (H, E->imax, E->jmax, E->kmax)) {
      DIE ("bad limits of the initialization loop");
   }

   k[0] *= h1;														// Scales k to work with indices directly.
   k[1] *= h2;
   k[2] *= h3;
   for (int i = E->imin ; i <= E->imax ; ++i)
   for (int j = E->jmin ; j <= E->jmax ; ++j)
   for (int l = E->kmin ; l <= E->kmax ; ++l) {
      double phase = k[0]*i + k[1]*j + k[2]*l;

      mv_fx(E, i, j, l) += e[0]*sin (phase - 0.5*k[0]) + type*h[0]*cos (phase - 0.5*k[0]);
      mv_fy(E, i, j, l) += e[1]*sin (phase - 0.5*k[1]) + type*h[1]*cos (phase - 0.5*k[1]);
      mv_fz(E, i, j, l) += e[2]*sin (phase - 0.5*k[2]) + type*h[2]*cos (phase - 0.5*k[2]);

      phase += 0.5*omega*tau;
      mv_fx(H, i, j, l) += h[0]*sin (phase - 0.5*k[1] - 0.5*k[2]) - type*e[0]*cos (phase - 0.5*k[1] - 0.5*k[2]);
      mv_fy(H, i, j, l) += h[1]*sin (phase - 0.5*k[0] - 0.5*k[2]) - type*e[1]*cos (phase - 0.5*k[0] - 0.5*k[2]);
      mv_fz(H, i, j, l) += h[2]*sin (phase - 0.5*k[0] - 0.5*k[1]) - type*e[2]*cos (phase - 0.5*k[0] - 0.5*k[1]);
   }
}
