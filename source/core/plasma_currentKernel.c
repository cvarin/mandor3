#ifndef MC_PLASMA_CURRENTKERNEL_HEADER
#define MC_PLASMA_CURRENTKERNEL_HEADER		///< \internal Guard.

/*
 * THIS MODULE IS INCLUDED INTO plasma.c AND INHERITS ALL DOMAIN SIZES FROM THERE.
 *
 * This module is a plasma_currentKernel: input into current density by single particle.
 * Particle does not cross boundaries of cell. All other cases are reduced to multiple
 * call to this subroutine. Current density is found doing intergration of transferred
 * charge and dividing by time step.
 */


// ---------------------------------------------------------------------------
/// Multi-stage-defined approach to PIC distribution of a current component is used.
/// Goal of it is to remove unnecessary computations explicitly for 1D/2D runs (in this
/// case we no longer need few variables (tmpX, tmpY, tmpZ) and few interpolations along
/// unreferenced dimension(s)).
// ---------------------------------------------------------------------------
static inline void
plasmaCurrents_kernel (const int i, const int j, const int k,
                       const double x1, const double x2,
                       const double y1, const double y2,
                       const double z1, const double z2,
                       const double rho_part)
{
  const double dX = x2 - x1;
  const double dY = y2 - y1;
  const double dZ = z2 - z1;

#if mc_have_x != 0
  const double tmpX = 0.5*(x1 + x2) - i*h1;
#endif

#if mc_have_y != 0
  const double tmpY = 0.5*(y1 + y2) - j*h2;
#endif

#if mc_have_z != 0
  const double tmpZ = 0.5*(z1 + z2) - k*h3;
#endif

#if defined (MC_CURRENT_KERNEL_RANGE_TEST) && MC_CURRENT_KERNEL_RANGE_TEST
  // Checks violation of the plasma_J array.
  ENSURE (!mf_mesh_pointIsOutside (plasma_J, i, j, k) && !mf_mesh_pointIsOutside (plasma_J, i+1, j+1, k+1),
              "Array '%s' violation: %d %d %d / %d %d %d.", plasma_J->name, i, j, k, i + 1, j + 1, k + 1);
#endif

/*
 *     oooooo             |  Original Jx distribution routine is:
 *        M    ooo   ooo  |
 *        M     "Mo oM"   |  mv_fx(J, i+1, j,   k)   += rho_part*dX/(tau*h2*h3)*( (h2 - tmpY)*(h3 - tmpZ) + dY*dZ/12 );
 *  M     M       M"M     |  mv_fx(J, i+1, j+1, k)   += rho_part*dX/(tau*h2*h3)*( tmpY*       (h3 - tmpZ) - dY*dZ/12 );
 *  Mo   oM     o"   "o   |  mv_fx(J, i+1, j,   k+1) += rho_part*dX/(tau*h2*h3)*( (h2 - tmpY)*tmpZ - dY*dZ/12 );
 *    """"    """"   """" |  mv_fx(J, i+1, j+1, k+1) += rho_part*dX/(tau*h2*h3)*( tmpY*       tmpZ + dY*dZ/12 );
 */

#if mc_have_y
  #define mf_Jx_distributeY(k, h3, expr_tmpZ, expr_dYdZ)                                                	\
    mv_fx (plasma_J, i+1, j,   (k)) += rho_part*dX/(tau*h2*(h3))*( (h2 - tmpY)*(expr_tmpZ) + (expr_dYdZ) );    	\
    mv_fx (plasma_J, i+1, j+1, (k)) += rho_part*dX/(tau*h2*(h3))*( tmpY*(expr_tmpZ) - (expr_dYdZ) );
#else
  #define mf_Jx_distributeY(k, h3, expr_tmpZ, expr_dYdZ)                                                	\
    mv_fx (plasma_J, i+1, j,   (k)) += rho_part*dX/(tau*(h3))*(expr_tmpZ);
#endif

#if mc_have_z
  mf_Jx_distributeY (k, h3, h3 - tmpZ, dY*dZ/12);
  mf_Jx_distributeY (k+1, h3, tmpZ, - dY*dZ/12);
#else
  mf_Jx_distributeY (k, 1, 1, 0);
#endif

/*
 *     oooooo             |  Original routine is:
 *        M    ooo   ooo  |
 *        M     "o   o"   |  mv_fy(J, i,   j+1, k)   += rho_part*dY/(tau*h1*h3)*( (h1 - tmpX)*(h3 - tmpZ) + dX*dZ/12 );
 *  M     M      M  o"    |  mv_fy(J, i+1, j+1, k)   += rho_part*dY/(tau*h1*h3)*( tmpX*(h3 - tmpZ) - dX*dZ/12 );
 *  Mo   oM       Mo"     |  mv_fy(J, i,   j+1, k+1) += rho_part*dY/(tau*h1*h3)*( (h1 - tmpX)*tmpZ - dX*dZ/12 );
 *    """"        o"      |  mv_fy(J, i+1, j+1, k+1) += rho_part*dY/(tau*h1*h3)*( tmpX*tmpZ + dX*dZ/12 );
 *             oooMo
 */

#if mc_have_x
  #define mf_Jy_distributeX(k, h3, expr_tmpZ, expr_dXdZ)                                                	\
    mv_fy (plasma_J, i,   j+1, (k)) += rho_part*dY/(tau*h1*(h3))*( (h1 - tmpX)*(expr_tmpZ) + (expr_dXdZ) );    	\
    mv_fy (plasma_J, i+1, j+1, (k)) += rho_part*dY/(tau*h1*(h3))*( tmpX*(expr_tmpZ) - (expr_dXdZ) );
#else
  #define mf_Jy_distributeX(k, h3, expr_tmpZ, expr_dXdZ)                                                	\
    mv_fy (plasma_J, i, j+1,   (k)) += rho_part*dY/(tau*(h3))*(expr_tmpZ);
#endif

#if mc_have_z
  mf_Jy_distributeX (k, h3, h3 - tmpZ, dX*dZ/12);
  mf_Jy_distributeX (k+1, h3, tmpZ, - dX*dZ/12);
#else
  mf_Jy_distributeX (k, 1, 1, 0);
#endif

/*
 *     oooooo             |  Original Jz distribution routine is:
 *        M     ooooooo   |
 *        M     M   o"    |  mv_fz(J, i,   j,   k+1) += rho_part*dZ/(tau*h1*h2)*( (h1 - tmpX)*(h2 - tmpY) + dX*dY/12 );
 *  M     M       oM"     |  mv_fz(J, i+1, j,   k+1) += rho_part*dZ/(tau*h1*h2)*( tmpX*(h2 - tmpY) - dX*dY/12 );
 *  Mo   oM     oM"   M   |  mv_fz(J, i,   j+1, k+1) += rho_part*dZ/(tau*h1*h2)*( (h1 - tmpX)*tmpY - dX*dY/12 );
 *    """"      """""""   |  mv_fz(J, i+1, j+1, k+1) += rho_part*dZ/(tau*h1*h2)*( tmpX*tmpY + dX*dY/12 );
 */

#if mc_have_x
  #define mf_Jz_distributeX(j, h2, expr_tmpY, expr_dXdY)                                              		\
    mv_fz (plasma_J, i,   (j), k+1) += rho_part*dZ/(tau*h1*(h2))*( (h1 - tmpX)*(expr_tmpY) + (expr_dXdY) );  	\
    mv_fz (plasma_J, i+1, (j), k+1) += rho_part*dZ/(tau*h1*(h2))*( tmpX*(expr_tmpY) - (expr_dXdY) );
#else
  #define mf_Jz_distributeX(j, h2, expr_tmpY, expr_dXdY)                                                	\
    mv_fz (plasma_J, i, (j),   k+1) += rho_part*dZ/(tau*(h2))*(expr_tmpY);
#endif

#if mc_have_y
  mf_Jz_distributeX (j, h2, h2 - tmpY, dX*dY/12);
  mf_Jz_distributeX (j+1, h2, tmpY, - dX*dY/12);
#else
  mf_Jz_distributeX (j, 1, 1, 0);
#endif
}

#undef mf_Jx_distributeY
#undef mf_Jy_distributeX
#undef mf_Jz_distributeX


// ---------------------------------------------------------------------------
/// Swaps two integer numbers.
// ---------------------------------------------------------------------------
#define mf_swap_int(x, y)  			\
{                       			\
  const int tmp = x;    			\
  x = y;                			\
  y = tmp;              			\
}

// ---------------------------------------------------------------------------
/// Swaps two double numbers.
// ---------------------------------------------------------------------------
#define mf_swap_double(x, y)  			\
{                          			\
  const double tmp = x;    			\
  x = y;                   			\
  y = tmp;                 			\
}

// ---------------------------------------------------------------------------
/// \brief Takes particle's trajectory, splits it by cell boundaries and evaluates input of each piece into charge transferred through cells' boundary
/// (or current density). Index \b C stands for \b Cross, \b T for \b Tmp. \b Warning: particle is supposed to move to neighbour cell only - it means
/// di, dj, dk are only 0 or 1!
// ---------------------------------------------------------------------------
static void
plasmaCurrents_multiCellStep (const int i1, const int j1, const int k1,
                              const int i2, const int j2, const int k2,
                              const double x1, const double y1, const double z1,
                              const double x2, const double y2, const double z2,
                              const double rho_part)
{
  const int di = abs (i2 - i1);
  const int dj = abs (j2 - j1);
  const int dk = abs (k2 - k1);

#if defined (MC_CURRENT_KERNEL_RANGE_TEST) && MC_CURRENT_KERNEL_RANGE_TEST
  // Checks violation of the plasma_J array.
  ENSURE (di < 2 && dj < 2 && dk < 2,
             "Too big step for current density evaluator: %d %d %d -> %d %d %d   |   %.3e %.3e %.3e -> %.3e %.3e %.3e",
             i1, j1, k1, i2, j2, k2, x1, y1, z1, x2, y2, z2);
#endif

  if (di + dj + dk == 2)
  {
    int    i = i1, j = j1, k = k1;
    double xA, xB, yA, yB, zA, zB;
    const double sigmaX = di*((i1 + i2 + 1)*0.5*h1 - x1)/(x2 - x1 + ((1 - di)<<16));
    const double sigmaY = dj*((j1 + j2 + 1)*0.5*h2 - y1)/(y2 - y1 + ((1 - dj)<<16));
    const double sigmaZ = dk*((k1 + k2 + 1)*0.5*h3 - z1)/(z2 - z1 + ((1 - dk)<<16));
    double kappa1 = sigmaX + (1 - di)*sigmaY;
    double kappa2 = sigmaZ + (1 - dk)*sigmaY;
    int deltaIndex[2][3] = {{0, 0, 0}, {0, 0, 0}};
    int step1 = 0;

    deltaIndex[0][0] = di*(i2 - i1);
    deltaIndex[0][1] = (1 - di)*(j2 - j1);
    deltaIndex[1][1] = (1 - dk)*(j2 - j1);
    deltaIndex[1][2] = dk*(k2 - k1);

    if (kappa1 > kappa2)
    {
      mf_swap_double (kappa1, kappa2);
      step1 = 1;
    }

    xA = x1 + kappa1*(x2 - x1);
    yA = y1 + kappa1*(y2 - y1);
    zA = z1 + kappa1*(z2 - z1);
    xB = x1 + kappa2*(x2 - x1);
    yB = y1 + kappa2*(y2 - y1);
    zB = z1 + kappa2*(z2 - z1);

    // 3 stage pass (2 boundary crossing).
    plasmaCurrents_kernel (i, j, k, x1, xA, y1, yA, z1, zA, rho_part);

    i += deltaIndex[step1][0];
    j += deltaIndex[step1][1];
    k += deltaIndex[step1][2];
    plasmaCurrents_kernel (i, j, k, xA, xB, yA, yB, zA, zB, rho_part);

    i += deltaIndex[1 - step1][0];
    j += deltaIndex[1 - step1][1];
    k += deltaIndex[1 - step1][2];
    plasmaCurrents_kernel (i, j, k, xB, x2, yB, y2, zB, z2, rho_part);

    return;
  }
  else
  {
    /* XYZ */
    const double xC = (i1 + i2 + 1)*0.5*h1;
    const double yC = (j1 + j2 + 1)*0.5*h2;
    const double zC = (k1 + k2 + 1)*0.5*h3;
    double sigma1 = (xC - x1)/(x2 - x1);
    double sigma2 = (yC - y1)/(y2 - y1);
    double sigma3 = (zC - z1)/(z2 - z1);
    double rCross[3][3];
    int dIndex[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    int    step1 = 0, step2 = 1, step3 = 2, i = i1, j = j1, k = k1;

    rCross[0][0] = xC;													// X-plane cross point.
    rCross[0][1] = y1 + sigma1*(y2 - y1);
    rCross[0][2] = z1 + sigma1*(z2 - z1);
    dIndex[0][0] = i2 - i1;

    rCross[1][0] = x1 + sigma2*(x2 - x1);										// Y-plane cross point.
    rCross[1][1] = yC;
    rCross[1][2] = z1 + sigma2*(z2 - z1);
    dIndex[1][1] = j2 - j1;

    rCross[2][0] = x1 + sigma3*(x2 - x1);										// Z-plane cross point.
    rCross[2][1] = y1 + sigma3*(y2 - y1);
    rCross[2][2] = zC;
    dIndex[2][2] = k2 - k1;

    if (sigma1 > sigma2)
    {
      mf_swap_int (step1, step2);
      mf_swap_double (sigma1, sigma2);
    }

    if (sigma1 > sigma3)
    {
      mf_swap_int (step1, step3);
      mf_swap_double (sigma1, sigma3);
    }

    if (sigma2 > sigma3)
    {
      mf_swap_int (step2, step3);
      mf_swap_double (sigma2, sigma3);
    }

    // 4-stage pass (crosses all boundaries).
    plasmaCurrents_kernel (i, j, k, x1, rCross[step1][0], y1, rCross[step1][1], z1, rCross[step1][2], rho_part);
    i += dIndex[step1][0];
    j += dIndex[step1][1];
    k += dIndex[step1][2];

    plasmaCurrents_kernel (i, j, k, rCross[step1][0], rCross[step2][0], rCross[step1][1], rCross[step2][1], rCross[step1][2], rCross[step2][2], rho_part);
    i += dIndex[step2][0];
    j += dIndex[step2][1];
    k += dIndex[step2][2];

    plasmaCurrents_kernel (i, j, k, rCross[step2][0], rCross[step3][0], rCross[step2][1], rCross[step3][1], rCross[step2][2], rCross[step3][2], rho_part);
    plasmaCurrents_kernel (i2, j2, k2, rCross[step3][0], x2, rCross[step3][1], y2, rCross[step3][2], z2, rho_part);
  }
}

#endif
