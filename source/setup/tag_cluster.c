/** \file tag_cluster.c
  * Cold uniform charge-neutral cluster creation (see tag_cluster.h).
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "log.h"
#include "type_mesh.h"
#include "type_marker.h"

#include "misc_units.h"
#include "misc_cfgReader.h"

#include "main.h"
#include "tag_cluster.h"

enum {eLinear, eGaussian} profiles;	///< Types of the profile.

// ---------------------------------------------------------------------------
/// Packed parameters of cluster.
// ---------------------------------------------------------------------------
static struct
{
  double qDivM;				///< Charge to mass ratio.
  double x, y, z;			///< Center of the cluster.
  double R, rho;			///< External radius and density in the center of the cluster.
  int    nx, ny, nz;			///< Markers sampling frequency (particles per axis per mesh cell).

  int    profile;			///< Type of the density profile (one of the enumerated above).
  double rho2, L;			///< Either density on the edge or gaussian scale in exp (-r^2/scale^2).
} cl;

// ---------------------------------------------------------------------------
/// Returns \b approximate memory consuming by cluster (nx, ny and nz are normalized on 1 (for 1D/2D) by caller).
// ---------------------------------------------------------------------------
// XXX static
double
tag_clusterEstimate (double R, int nx, int ny, int nz)
{
  double l = 0, volume = ((h1 - 1)*mc_have_x + 1)*((h2 - 1)*mc_have_y + 1)*((h3 - 1)*mc_have_z + 1);
  switch (mc_have_x + mc_have_y + mc_have_z)
  {
    case 1:	l = 2*R/volume*nx*ny*nz;			break;
    case 2:	l = mc_pi*R*R/volume*nx*ny*nz;			break;
    case 3:	l = 4.0/3.0*mc_pi*R*R*R/volume*nx*ny*nz;	break;

    default:
      DIE ("bad flags mc_have_x = %d, mc_have_y = %d, mc_have_z = %d",
           mc_have_x, mc_have_y, mc_have_z);
    break;
  }
  return l*sizeof (marker_t);
}

// ---------------------------------------------------------------------------
/// Reads parameters for tag [cluster] and creates cluster itself.
// ---------------------------------------------------------------------------
// XXX static
void
tag_cluster_read (FILE *fp)
{
  const double micron = units (mc_micron);										// Micron in dimensionless coordinates.

  cl.x = cfg_readDouble (fp)*micron*mc_have_x;										// Center position.
  cl.y = cfg_readDouble (fp)*micron*mc_have_y;
  cl.z = cfg_readDouble (fp)*micron*mc_have_z;

  cl.R = cfg_readDouble (fp)*micron;											// Radius.
  ENSURE (cl.R > 4*(h1*mc_have_x + h2*mc_have_y + h3*mc_have_z),
          "badly resolved R (%f): h1 = %f, h2 = %f, h3 = %f",
          cl.R, h1*mc_have_x, h2*mc_have_y, h3*mc_have_z);

  double n = cfg_readDouble (fp), n2 = 10;										// Concentration of particles.

  const char *type = cfg_readWord (fp);											// Reads profile type.

  cl.profile = cfg_identifyWord (type, "linear", eLinear, "gaussian", eGaussian, mc_cfgTermGuesses);
  switch (cl.profile)
  {
    case eLinear:
      n2 = cfg_readOptDouble (fp);
      cl.L = cl.R;													// To turn r -> alpha \in [0, 1].
    break;

    case eGaussian:
      cl.L = cfg_readOptDouble (fp)*micron;
    break;

    default:
      DIE ("bad profile (%d / %s)", cl.profile, type);
  }

  const double q = cfg_readDouble (fp);											// Markers type - charge in |e|.
  const double M = cfg_readDouble (fp);											// Markers type - mass in |me|.
  ENSURE (M >= 0 && fabs (M) >= 1, "bad M = %e", M);
  cl.qDivM = q/M;													// Sets q/M after all checks.

  cl.nx = mc_have_x*cfg_readInt (fp) + (1 - mc_have_x);									// Resolution.
  cl.ny = mc_have_y*cfg_readInt (fp) + (1 - mc_have_y);
  cl.nz = mc_have_z*cfg_readInt (fp) + (1 - mc_have_z);
  ENSURE (cl.nx > 0 && cl.ny > 0 && cl.nz > 0,
          "bad 'particles per cell' parameters (%d, %d, %d)",
          cl.nx, cl.ny, cl.nz);

  cl.rho = (n > 0) ? n*units (mc_ne_critical) : - n;									// Evaluates n_\alpha in [cm^-3].
  cl.rho *= q*mc_CGS_e/units (mc_rho0);											// Evaluates rho in dimensionless units.
  ENSURE (q*cl.rho > 0.01, "q = %e, rho_center = %e", q, cl.rho);

  cl.rho2 = (n2 > 0) ? n2*units (mc_ne_critical) : - n2;								// Evaluates n_\alpha in [cm^-3].
  cl.rho2 *= q*mc_CGS_e/units (mc_rho0);										// Evaluates rho in dimensionless units.
  ENSURE (q*cl.rho2 > 0.01, "q = %e, rho_edge = %e", q, cl.rho2);

  say ("tag_cluster: ");
  say ("  - q = %f [e], M = %f [me]", q, M);
  say ("  - center at (%.3f, %.3f, %.3f) [r0] = (%.3f, %.3f, %.3f) [micron]", cl.x, cl.y, cl.z, cl.x/micron, cl.y/micron, cl.z/micron);
  say ("  - external radius = %.3f [r0] = %.3f [micron]", cl.R, cl.R/micron);
  n = cl.rho*units (mc_rho0)/(q*mc_CGS_e);
  say ("    o nCenter = %e [cm^{-3}] = %e [n_e_cr] = %e [n_alpha_cr]", n, n/units (mc_ne_critical), n*q*q/(M*units (mc_ne_critical)));
  say ("    o rhoCenter = %e [rho0]", cl.rho);
  switch (cl.profile)
  {
    case eLinear:
      say ("  - linear density profile:");
      n = cl.rho2*units (mc_rho0)/(q*mc_CGS_e);
      say ("    o nEdge = %e [cm^{-3}] = %e [n_e_cr] = %e [n_alpha_cr]", n, n/units (mc_ne_critical), n*q*q/(M*units (mc_ne_critical)));
      say ("    o rhoEdge = %e [rho0]", cl.rho2);
    break;

    case eGaussian:
      say ("  - gaussian density profile (exp (-r^2/L^2):");
      say ("    o L = %e [r0] = %e [micron]", cl.L, cl.L/micron);
      n = cl.rho*exp (- cl.R*cl.R/(cl.L*cl.L))*units (mc_rho0)/(q*mc_CGS_e);
      say ("    o nEdge = %e [cm^{-3}] = %e [n_e_cr] = %e [n_alpha_cr]", n, n/units (mc_ne_critical), n*q*q/(M*units (mc_ne_critical)));
      say ("    o rhoEdge = %e [rho0]", cl.rho2);
    break;
  }
}

// ---------------------------------------------------------------------------
/// Reads parameters for tag [cluster] and creates cluster itself.
// ---------------------------------------------------------------------------
double
tag_cluster (FILE *fp)
{
#if 0
  tag_cluster_read (fp);

  if (memEstimateOnly)
    return tag_clusterEstimate (cl.R, cl.nx, cl.ny, cl.nz);

  const double dx = h1/(double) cl.nx;											// Markers placing intervals.
  const double dy = h2/(double) cl.ny;
  const double dz = h3/(double) cl.nz;

  long long int particles = 0;

  reg_t clstr = {{(cl.x - cl.R)/h1 - 2, (cl.y - cl.R)/h2 - 2, (cl.z - cl.R)/h3 - 2},					// Region aroung cluster.
                 {(cl.x + cl.R)/h1 + 2, (cl.y + cl.R)/h2 + 2, (cl.z + cl.R)/h3 + 2}};
  reg_t cp = {{cpu_min[0], cpu_min[1], cpu_min[2]}, {cpu_max[0] - mc_have_x, cpu_max[1] - mc_have_y, cpu_max[2] - mc_have_z}};	// Node domain.

  reg_overlap (&clstr, &cp, NULL);											// Computes overalpped part.
  mf_reg_collapse (&clstr);												// Adjustes for 1D/2D.

  cache_startChapter ();												// Starts new chapter.
  cache_setSpecie (cl.qDivM);

  for (int i = clstr.min[0] ; i <= clstr.max[0] ; ++i)
    for (int j = clstr.min[1] ; j <= clstr.max[1] ; ++j)
      for (int k = clstr.min[2] ; k <= clstr.max[2] ; ++k)
        for (int di = 0 ; di < cl.nx ; ++di)
          for (int dj = 0 ; dj < cl.ny ; ++dj)
            for (int dk = 0 ; dk < cl.nz ; ++dk)
            {
              double x = i*h1 + (di + 0.5)*dx;
              double y = j*h2 + (dj + 0.5)*dy;
              double z = k*h3 + (dk + 0.5)*dz;
              double r2 = mc_have_x*(x - cl.x)*(x - cl.x) + mc_have_y*(y - cl.y)*(y - cl.y) + mc_have_z*(z - cl.z)*(z - cl.z);

              if (r2 <= cl.R*cl.R)
              {
                marker_t *marker;
                marker = cache_getParticle ();										// Claims pointer on cached marker.
                marker->x = x*mc_have_x;
                marker->y = y*mc_have_y;
                marker->z = z*mc_have_z;
                marker->vx = 0;
                marker->vy = 0;
                marker->vz = 0;
                marker->rho = sqrt (r2)/cl.L;										// Computes interpolation coefficient.
                ++particles;
              }
            }

  cache_flush ();													// Finalizes chapter.

  const double weight = 1.0/((double) cl.nx*cl.ny*cl.nz);								// Charge density -> marker density factor.
  const int ID = cache_lastChapter ();
  markerIterator_t page;
  for (markerPage_first (ID, &page) ; page.df ; markerPage_next (&page))						// Computes charge density.
  {
    switch (cl.profile)
    {
      case eLinear:
        for (int p = 0 ; p < page.N ; ++p)
          page.df[p].rho = weight*(cl.rho + (cl.rho2 - cl.rho)*page.df[p].rho);
      break;

      case eGaussian:
        for (int p = 0 ; p < page.N ; ++p)
          page.df[p].rho = weight*cl.rho*exp (- page.df[p].rho*page.df[p].rho);
      break;

      default:
        DIE ("unknow profile");
      break;
    }
  }

  return ((double) particles)*sizeof (marker_t);
#endif
  DIE ("not re-implemented");
}
