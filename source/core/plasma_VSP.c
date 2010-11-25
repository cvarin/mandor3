/** \file plasma_VSP.c
  * Routines for transformation of the PIC computed current/charge density to the VSP one.
  */

#include <math.h>
#include <assert.h>

#include "plasma_VSP.h"

#include "log.h"

/// Convolution parameter.
static double VSPalpha = 1.0;

// ---------------------------------------------------------------------------
/// Configures convolution kernel.
// ---------------------------------------------------------------------------
void
VSP_configure (double Alpha)
{
  ENSURE (Alpha > - 0.1 && Alpha <= 1.1,
          "too weird parameter; test it and report to maintainer to include");

  VSPalpha = Alpha;
  say ("VSP_configure:\n  o alpha = %f\n  o kernel [%f %f %f]", Alpha, 0.5*(1 - Alpha), Alpha, 0.5*(1 - Alpha));
}

// ---------------------------------------------------------------------------
/// Convolutes charge density with \f$(\beta, \alpha, \beta)\f$ function where \f$2*\beta + \alpha = 1\f$.
// ---------------------------------------------------------------------------
void
VSP_scalarSweep (meshDouble_t *rho)
{
  if (fabs (VSPalpha - 1.0) < 1e-10)
    return;

  assert ( (!mf_mesh_pointIsOutside (rho, cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2)) && 				// Checks sizes before accessing mesh.
           (!mf_mesh_pointIsOutside (rho, cpu_min[0] + 2, cpu_min[1] + 2, cpu_min[2] + 2)) );

  const double alpha = VSPalpha, beta = 0.5*(1.0 - alpha);

#if mc_have_z
  for (int i = rho->imin ; i <= rho->imax ; ++i)
    for (int j = rho->jmin ; j <= rho->jmax ; ++j)
    {
      double untaintedRho = mv_f (rho, i, j, cpu_min[2] - 2);
      mv_f (rho, i, j, cpu_min[2] - 2) = beta*mv_f (rho, i, j, cpu_min[2] - 1) + alpha*mv_f (rho, i, j, cpu_min[2] - 2);
      for (int k = cpu_min[2] - 1 ; k <= cpu_max[2] + 1; ++k)
      {
        double VSPrho = alpha*(mv_f (rho, i, j, k)) + beta*(untaintedRho + mv_f (rho, i, j, k + 1));
        untaintedRho = mv_f (rho, i, j, k);
        mv_f (rho, i, j, k) = VSPrho;
      }
      mv_f (rho, i, j, cpu_max[2] + 2) = beta*untaintedRho + alpha*mv_f (rho, i, j, cpu_max[2] + 2);
    }
#endif

#if mc_have_y
  for (int i = rho->imin ; i <= rho->imax ; ++i)
    for (int k = rho->kmin ; k <= rho->kmax ; ++k)
    {
      double untaintedRho = mv_f (rho, i, cpu_min[1] - 2, k);
      mv_f (rho, i, cpu_min[1] - 2, k) = beta*mv_f (rho, i, cpu_min[1] - 1, k) + alpha*mv_f (rho, i, cpu_min[1] - 2, k);
      for (int j = cpu_min[1] - 1 ; j <= cpu_max[1] + 1; ++j)
      {
        double VSPrho = alpha*(mv_f (rho, i, j, k)) + beta*(untaintedRho + mv_f (rho, i, j + 1, k));
        untaintedRho = mv_f (rho, i, j, k);
        mv_f (rho, i, j, k) = VSPrho;
      }
      mv_f (rho, i, cpu_max[1] + 2, k) = beta*untaintedRho + alpha*mv_f (rho, i, cpu_max[1] + 2, k);
    }
#endif

#if mc_have_x
  for (int j = rho->jmin ; j <= rho->jmax ; ++j)
    for (int k = rho->kmin ; k <= rho->kmax ; ++k)
    {
      double untaintedRho = mv_f (rho, cpu_min[0] - 2, j, k);
      mv_f (rho, cpu_min[0] - 2, j, k) = beta*mv_f (rho, cpu_min[0] - 1, j, k) + alpha*mv_f (rho, cpu_min[0] - 2, j, k);
      for (int i = cpu_min[0] - 1 ; i <= cpu_max[0] + 1; ++i)
      {
        double VSPrho = alpha*(mv_f (rho, i, j, k)) + beta*(untaintedRho + mv_f (rho, i + 1, j, k));
        untaintedRho = mv_f (rho, i, j, k);
        mv_f (rho, i, j, k) = VSPrho;
      }
      mv_f (rho, cpu_max[0] + 2, j, k) = beta*untaintedRho + alpha*mv_f (rho, cpu_max[0] + 2, j, k);
    }
#endif
}

// ---------------------------------------------------------------------------
/// Convolutes current density with \f$(\beta, \alpha, \beta)\f$ function where \f$2*\beta + \alpha = 1\f$.
// ---------------------------------------------------------------------------
void
VSP_vectorSweep (meshVec_t *J)
{
  if (fabs (VSPalpha - 1.0) < 1e-10)
    return;

  assert ( (!mf_mesh_pointIsOutside (J, cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2)) && 					// Checks sizes before accessing mesh.
           (!mf_mesh_pointIsOutside (J, cpu_min[0] + 2, cpu_min[1] + 2, cpu_min[2] + 2)) );

  const double alpha = VSPalpha, beta = 0.5*(1.0 - alpha);

#if mc_have_z
  for (int i = J->imin ; i <= J->imax ; ++i)
    for (int j = J->jmin ; j <= J->jmax ; ++j)
    {
      vec3D_t untaintedJ = mv_v (J, i, j, cpu_min[2] - 2);
      mv_fx (J, i, j, cpu_min[2] - 2) = beta*mv_fx (J, i, j, cpu_min[2] - 1) + alpha*mv_fx (J, i, j, cpu_min[2] - 2);
      mv_fy (J, i, j, cpu_min[2] - 2) = beta*mv_fy (J, i, j, cpu_min[2] - 1) + alpha*mv_fy (J, i, j, cpu_min[2] - 2);
      mv_fz (J, i, j, cpu_min[2] - 2) = beta*mv_fz (J, i, j, cpu_min[2] - 1) + alpha*mv_fz (J, i, j, cpu_min[2] - 2);
      for (int k = cpu_min[2] - 1 ; k <= cpu_max[2] + 1; ++k)
      {
        double VSPJx = alpha*(mv_fx (J, i, j, k)) + beta*(untaintedJ.x + mv_fx (J, i, j, k + 1));
        double VSPJy = alpha*(mv_fy (J, i, j, k)) + beta*(untaintedJ.y + mv_fy (J, i, j, k + 1));
        double VSPJz = alpha*(mv_fz (J, i, j, k)) + beta*(untaintedJ.z + mv_fz (J, i, j, k + 1));
        untaintedJ = mv_v (J, i, j, k);
        mv_fx (J, i, j, k) = VSPJx;
        mv_fy (J, i, j, k) = VSPJy;
        mv_fz (J, i, j, k) = VSPJz;
      }
      mv_fx (J, i, j, cpu_max[2] + 2) = beta*untaintedJ.x + alpha*mv_fx (J, i, j, cpu_max[2] + 2);
      mv_fy (J, i, j, cpu_max[2] + 2) = beta*untaintedJ.y + alpha*mv_fy (J, i, j, cpu_max[2] + 2);
      mv_fz (J, i, j, cpu_max[2] + 2) = beta*untaintedJ.z + alpha*mv_fz (J, i, j, cpu_max[2] + 2);
    }
#endif

#if mc_have_y
  for (int i = J->imin ; i <= J->imax ; ++i)
    for (int k = J->kmin ; k <= J->kmax ; ++k)
    {
      vec3D_t untaintedJ = mv_v (J, i, cpu_min[1] - 2, k);
      mv_fx (J, i, cpu_min[1] - 2, k) = beta*mv_fx (J, i, cpu_min[1] - 1, k) + alpha*mv_fx (J, i, cpu_min[1] - 2, k);
      mv_fy (J, i, cpu_min[1] - 2, k) = beta*mv_fy (J, i, cpu_min[1] - 1, k) + alpha*mv_fy (J, i, cpu_min[1] - 2, k);
      mv_fz (J, i, cpu_min[1] - 2, k) = beta*mv_fz (J, i, cpu_min[1] - 1, k) + alpha*mv_fz (J, i, cpu_min[1] - 2, k);
      for (int j = cpu_min[1] - 1 ; j <= cpu_max[1] + 1; ++j)
      {
        double VSPJx = alpha*(mv_fx (J, i, j, k)) + beta*(untaintedJ.x + mv_fx (J, i, j + 1, k));
        double VSPJy = alpha*(mv_fy (J, i, j, k)) + beta*(untaintedJ.y + mv_fy (J, i, j + 1, k));
        double VSPJz = alpha*(mv_fz (J, i, j, k)) + beta*(untaintedJ.z + mv_fz (J, i, j + 1, k));
        untaintedJ = mv_f (J, i, j, k);
        mv_fx (J, i, j, k) = VSPJx;
        mv_fy (J, i, j, k) = VSPJy;
        mv_fz (J, i, j, k) = VSPJz;
      }
      mv_fx (J, i, cpu_max[1] + 2, k) = beta*untaintedJ.x + alpha*mv_fx (J, i, cpu_max[1] + 2, k);
      mv_fy (J, i, cpu_max[1] + 2, k) = beta*untaintedJ.y + alpha*mv_fy (J, i, cpu_max[1] + 2, k);
      mv_fz (J, i, cpu_max[1] + 2, k) = beta*untaintedJ.z + alpha*mv_fz (J, i, cpu_max[1] + 2, k);
    }
#endif

#if mc_have_x
  for (int j = J->jmin ; j <= J->jmax ; ++j)
    for (int k = J->kmin ; k <= J->kmax ; ++k)
    {
      vec3D_t untaintedJ = mv_v (J, cpu_min[0] - 2, j, k);
      mv_fx (J, cpu_min[0] - 2, j, k) = beta*mv_fx (J, cpu_min[0] - 1, j, k) + alpha*mv_fx (J, cpu_min[0] - 2, j, k);
      mv_fy (J, cpu_min[0] - 2, j, k) = beta*mv_fy (J, cpu_min[0] - 1, j, k) + alpha*mv_fy (J, cpu_min[0] - 2, j, k);
      mv_fz (J, cpu_min[0] - 2, j, k) = beta*mv_fz (J, cpu_min[0] - 1, j, k) + alpha*mv_fz (J, cpu_min[0] - 2, j, k);
      for (int i = cpu_min[0] - 1 ; i <= cpu_max[0] + 1; ++i)
      {
        double VSPJx = alpha*(mv_fx (J, i, j, k)) + beta*(untaintedJ.x + mv_fx (J, i + 1, j, k));
        double VSPJy = alpha*(mv_fy (J, i, j, k)) + beta*(untaintedJ.y + mv_fy (J, i + 1, j, k));
        double VSPJz = alpha*(mv_fz (J, i, j, k)) + beta*(untaintedJ.z + mv_fz (J, i + 1, j, k));
        untaintedJ = mv_f (J, i, j, k);
        mv_fx (J, i, j, k) = VSPJx;
        mv_fy (J, i, j, k) = VSPJy;
        mv_fz (J, i, j, k) = VSPJz;
      }
      mv_fx (J, cpu_max[0] + 2, j, k) = beta*untaintedJ.x + alpha*mv_fx (J, cpu_max[0] + 2, j, k);
      mv_fy (J, cpu_max[0] + 2, j, k) = beta*untaintedJ.y + alpha*mv_fy (J, cpu_max[0] + 2, j, k);
      mv_fz (J, cpu_max[0] + 2, j, k) = beta*untaintedJ.z + alpha*mv_fz (J, cpu_max[0] + 2, j, k);
    }
#endif
}
