/** \file tag_photoelectrons.c
  * Fills domain with uniform photo-ionized plasma (see 'tag_photoelectrons.h'
  * for theoretical details, formulas, etc).
  *
  * \warning Staggered allocation (to increase resolution in velocity space by
  * loosing resolution in coordinate space) is experimental. I observe high
  * frequency noise excitation when I use it; possible solution is to remove
  * mirror/spatial coupling (otherwise I generate microstreams directed in opposite
  * direction with steps about \f$ (N_x\cdot h_1, N_y\cdot h_2, N_z\cdot h_3) \f$.
  *
  * \sa tag_photoDF_parameters(), tag_photoelectrons.h, tag_seed_PITS() and setup_distrMapper.h.
  */

#include <math.h>
#include <stdlib.h>

#include <mpi.h>

#include "misc_units.h"
#include "log.h"
#include "misc_cfgReader.h"

#include "setup/main.h"
#include "setup/plasma.h"
#include "tag_photoelectrons.h"
#include "setup_denavit.h"
#include "setup_distrMapper.h"

#define mc_nSamples 10000				///< Number of samples in theta-mapper.

/// Parameters of the distribution function.
static int    nx, ny, nz, N;				// Fine sub-cell meshing parameters.
static int    Nx, Ny, Nz;				// Staggered mesh steps.
static int    mirror, nRotations;			// Symmetry of the pattern flags/parameters.

/// \b set confirms the plasma allocation and \b ID holds plasma chapter reference number.
static int    set = 0, ID = -1;

/// Parameters of the plasma. Density is referenced through \f$ \omega^2_{pe} \f$. Avaliable through interface tag_photoDF_parameters().
static double V0 = 0, omega2_pe = 0;

static double tag_photoInit (const double rho, const int uniformWeight, double qDivM);

// ---------------------------------------------------------------------------
/// Returns parameter of the photo-DF to the modification routine (like perturbation creation, etc). Values of the param are predefined in tag_photoDF.h.
// ---------------------------------------------------------------------------
void
tag_photoDF_parameters (int param, void *pntr)
{
  ENSURE (set, "There are no photoDF plasma allocated.");

  switch (param)
  {
    case mc_photoDF_V0:
      *((double*) pntr) = V0;
    break;

    case mc_photoDF_omega2_pe:
      *((double*) pntr) = omega2_pe;
    break;

    case mc_photoDF_ID:
      *((int*) pntr) = ID;
    break;

    default:
      DIE ("unknown id; please update interface function");
    break;
  }
}

// ---------------------------------------------------------------------------
/// Creates quiet start photo-ionization produced plasma DF (for paper with Valery).
// ---------------------------------------------------------------------------
double
tag_photoelectrons (FILE *fp)
{
  int uniformWeight;													// Flag to choose sampling of the DF.
  double qDivM, chargeDensity, eta;

  nx = cfg_readInt (fp);												// Reads parameters.
  ny = cfg_readInt (fp);
  nz = cfg_readInt (fp);
  N  = cfg_readInt (fp);
  mirror = cfg_readInt (fp);
  nRotations = cfg_readInt (fp);

  eta = cfg_readDouble (fp);
  qDivM = cfg_readDouble (fp);
  V0 = cfg_readDouble (fp);

  Nx = cfg_readInt (fp);							// Reads staggered mesh steps.
  Ny = cfg_readInt (fp);
  Nz = cfg_readInt (fp);

  uniformWeight = cfg_readInt (fp) == 1;

  V0 = (V0 >= 0) ? V0 : sqrt (-2.0*V0*mc_CGS_eV/mc_CGS_m)/units (mc_v0);	// Converts velocity to our units.

  chargeDensity = (eta > 0) ? eta*units (mc_ne_critical)/qDivM : - eta;		// Evaluates n_\alpha[cm^-3].
  chargeDensity *= mc_CGS_e/units (mc_rho0);					// Evaluates rho in dimensionless units.
  if (qDivM*chargeDensity <= 0)
    chargeDensity *= - 1;

  // Updates omega_pe after updating rho.
  omega2_pe = 4*mc_pi * chargeDensity*units (mc_rho0) * qDivM*mc_CGS_e/mc_CGS_m * units (mc_t0)*units (mc_t0);

  if (nx <= 0 || ny <= 0 || nz <= 0 || Nx <= 0 || Ny <= 0 || Nz <= 0)
    DIE ("nx, ny, nz, Nx, Ny, Nz must be positive");

  if (dmn_max[0] % Nx || dmn_max[1] % Ny || dmn_max[2] % Nz)
    DIE ("staggered steps size doesn't fit domain size");

  ENSURE (nRotations >= 0, "nRotations (%d) must be >= 0", nRotations);

  mirror = (mirror != 0);							// Clamps mirror flag to \in {0, 1}.

  double rDebay = sqrt (V0*V0/(5.0*omega2_pe));					// Estimates debay scale from disp. eq..

  say ("tag_photoDF:");
  say ("  - photo-DF plasma component is added,");
  say ("  - V0 = %.3e (energy = %.3e eV),", V0, mc_CGS_m*pow (V0*units (mc_v0), 2)/2/mc_CGS_eV);
  say ("  - (omega_pe/omega_0)^2 = %.3e, q/M = %.3e, n_0 = %.3e [cm^-3],", omega2_pe/(4*mc_pi*mc_pi), qDivM, chargeDensity*units (mc_rho0)/mc_CGS_e);
  say ("  - resolution of debay scale (nodes on rD) is");
  say ("    %e (x), %e (y), %e (z),", rDebay/h1, rDebay/h2, rDebay/h3);
  say ("  - particles placement:");
  say ("    o uniform %s of markers,", (uniformWeight) ? "weight" : "spacing");
  say ("    o mirror: %d, rotations: %d,", mirror, nRotations);
  say ("    o %d x %d x %d sub-cell mesh,", nx, ny, nz);
  say ("    o %d x %d x %d staggered mesh steps,", Nx, Ny, Nz);
  say ("    o %d particles per cell,", nx*ny*nz*N*(1 + mirror)*(1 + nRotations));
  say ("    o %d particles per pattern,", nx*ny*nz*Nx*Ny*Nz*N*(1 + mirror)*(1 + nRotations));
  double N_part = nx*ny*nz*N*(dmn_max[0] - dmn_min[0] + 1 - mc_have_x)*
                             (dmn_max[1] - dmn_min[1] + 1 - mc_have_y)*
                             (dmn_max[2] - dmn_min[2] + 1 - mc_have_z)*(1.0 + mirror)*(1.0 + nRotations);
  say ("    o %.3e particles total,", N_part);

  if (memEstimateOnly)
  {
    set = 1;									// To permit seed parameters checks.
    return (N_part + nx*ny*nz*Nx*Ny*Nz*N*(1.0 + mirror)*(1.0 + nRotations))*sizeof (marker_t);
  }

  plasma_newObject ();								// Starts new chapter.

  double Ntotal = 0, Nlocal;							// Allocates particles and gets number.
  Nlocal = tag_photoInit (chargeDensity/(nx*ny*nz*N*(1.0 + mirror)*(1.0 + nRotations)), uniformWeight, qDivM);
  say ("    o %.3e particles added / cpu %d.", Nlocal, cpu_here);

  MPI_Allreduce (&Nlocal, &Ntotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);	// Checks total number of particles created.
  ENSURE (Ntotal == N_part, "lost %.1f particles", N_part - Ntotal);

  set = 1;
  // Size (pattern size included as storage).
  return (Nlocal + nx*ny*nz*Nx*Ny*Nz*N*(1.0 + mirror)*(1.0 + nRotations))*sizeof (marker_t);
}

// ---------------------------------------------------------------------------
/// \brief Angular distribution of the photo-electrons \f$ f(\theta) = \sin (\theta)\cdot\cos^2 (theta) \f$.
// ---------------------------------------------------------------------------
double
tag_photoDF_DF (double theta)
{
  return sin (theta)*cos (theta)*cos (theta);
}

// ---------------------------------------------------------------------------
/// \brief Uniform angular distribution of the markers (not photo-electrons) \f$ f(\theta) = \sin (\theta) \f$.
// ---------------------------------------------------------------------------
double
tag_photoDF_markers (double theta)
{
  return sin (theta);
}

// ---------------------------------------------------------------------------
/// Creates the DF itself used values precalculated at tag_photoDF(FILE *fp).
// ---------------------------------------------------------------------------
static double
tag_photoInit (const double rho, const int uniformWeight, double qDivM)
{
  int nMax = 0;
  int patternSize = nx*ny*nz*Nx*Ny*Nz*N*(1 + mirror)*(1 + nRotations);
  marker_t *pattern = (marker_t*) malloc (patternSize*sizeof (marker_t));
  FILE *fp = NULL;

  mapper_t *mapper;								// Creates mapper.
  if (uniformWeight)
    mapper = mapper_create (tag_photoDF_DF, 0, mc_pi, 1000);
  else
    mapper = mapper_create (tag_photoDF_markers, 0, mc_pi, 1000);

  if (!cpu_here)
  {
    fp = cfg_open ("output/photoDF_pattern.dat", "wt", "tag_photoInit");	// Saves pattern for examination.
    fprintf (fp, "variables = <greek>j/p</greek>, <greek>q/p</greek>, x, y, z, v<sub>x</sub>, v<sub>y</sub>, v<sub>z</sub>\nzone t=\"Basic pattern\", f = point\n");
  }

  int p = 0;									// Makes basic pattern.
  for (int i = 0 ; i < nx*Nx ; ++i)
    for (int j = 0 ; j < ny*Ny ; ++j)
      for (int k = 0 ; k < nz*Nz ; ++k)
        for (int l = 0 ; l < N ; ++l)
        {
          double phi, theta, tmp;
          denavit_createQuartet (p, nx*ny*nz*Nx*Ny*Nz*N, &tmp, &theta, &phi, &tmp);	// Gets Denavit-Wallsh coverage sample.

          theta /= 1.0 + mirror;						// Fills vz > 0 part (mirror will fill the rest).
          phi /= 1.0 + nRotations;						// Fills one slice of pie in V-space.

          pattern[p].x =  h1*(i + 0.5)/(double) nx;				// Uniform mesh in space.
          pattern[p].y =  h2*(j + 0.5)/(double) ny;
          pattern[p].z =  h3*(k + 0.5)/(double) nz;

          theta = mapper_invoke (theta, mapper);				// Morphes angle to get desired distribution.

          pattern[p].vx = V0*sin (theta)*cos(phi*2*mc_pi);
          pattern[p].vy = V0*sin (theta)*sin(phi*2*mc_pi);
          pattern[p].vz = V0*cos (theta);

          pattern[p].rho = rho;
          pattern[p].rho = qDivM;

          if (!cpu_here)
            fprintf (fp, "%e %e %e %e %e %e %e %e\n", atan2 (pattern[p].vx, pattern[p].vy)/mc_pi, acos (pattern[p].vz/V0)/mc_pi,
                     pattern[p].x, pattern[p].y, pattern[p].z, pattern[p].vx, pattern[p].vy, pattern[p].vz);

          ++p;
        }
  nMax = p;

  mapper_destroy (mapper);							// Destroys mapper to free memory.

  if (!uniformWeight)								// Correction of the weights of the markers.
  {
    double correction = 0;
    for (int i = 0 ; i < nMax ; ++i)						// Sets angular distribution of charge density.
    {
      pattern[i].rho = pattern[i].vz*pattern[i].vz;
      correction += pattern[i].rho;
    }

    correction = rho*nMax/correction;						// Sets proper magnitudes.
    for (int i = 0 ; i < nMax ; ++i)
      pattern[i].rho *= correction;
  }

  if (!cpu_here)
    fprintf (fp, "\nzone t=\"Mirrored\", f = point\n");
  for (int i = 0 ; i < nMax && mirror ; ++i)					// Adds mirrored basic pattern to the group.
  {
    pattern[p] = pattern[p-nMax];						// Full copy of all fields.

    /// \todo Remix mirror pattern (line allocate not in form pppp/p'p'p'p'p'/p_rotp_rot, but interleave) to avoid explicit excitation of the modes with \f$ lambda = N_x h_x \f$.
    pattern[p].z =  h3*Nz - pattern[p].z;					// Mirror: Z -> Z'.
    pattern[p].vz = - pattern[p].vz;						// Mirror: Vz - Vz'.

    if (!cpu_here)
      fprintf (fp, "%e %e %e %e %e %e %e %e\n", atan2 (pattern[p].vx, pattern[p].vy)/mc_pi, acos (pattern[p].vz/V0)/mc_pi,
               pattern[p].x, pattern[p].y, pattern[p].z, pattern[p].vx, pattern[p].vy, pattern[p].vz);

    ++p;
  }
  nMax = p;

  // Pattern += (basic pattern & mirror) rotated by 2*\pi*j/(nRotations + 1), j = 1, .., nRot.
  for (int l = 1 ; l <= nRotations ; ++l)
  {
    if (!cpu_here)
      fprintf (fp, "\nzone t=\"Rotated by %d<sup>o</sup>\", f = point\n", (int) l*360/(nRotations + 1));
    double phi = 2*mc_pi*l/(double)(nRotations + 1);
    for (int k = 0 ; k < nMax ; ++k)
    {
      pattern[p] =  pattern[k];							// Full copy of all fields.

      pattern[p].vx = pattern[k].vx*cos (phi) - pattern[k].vy*sin (phi);	// Rotation of the V vector.
      pattern[p].vy = pattern[k].vx*sin (phi) + pattern[k].vy*cos (phi);

      if (!cpu_here)
        fprintf (fp, "%e %e %e %e %e %e %e %e\n", atan2 (pattern[p].vx, pattern[p].vy)/mc_pi, acos (pattern[p].vz/V0)/mc_pi,
                 pattern[p].x, pattern[p].y, pattern[p].z, pattern[p].vx, pattern[p].vy, pattern[p].vz);

      ++p;
    }
  }
  if (!cpu_here)
    fclose (fp);

  ENSURE (p == patternSize,
          "number of initialized markers differs from number of allocated");

  double xMin = +1e100, yMin = +1e100, zMin = +1e100;				// Bounding box of the pattern.
  double xMax = -1e100, yMax = -1e100, zMax = -1e100;
  for (p = 0 ; p < patternSize ; ++p)
  {
    xMin = (pattern[p].x < xMin) ? pattern[p].x : xMin;
    yMin = (pattern[p].y < yMin) ? pattern[p].y : yMin;
    zMin = (pattern[p].z < zMin) ? pattern[p].z : zMin;
    xMax = (pattern[p].x > xMax) ? pattern[p].x : xMax;
    yMax = (pattern[p].y > yMax) ? pattern[p].y : yMax;
    zMax = (pattern[p].z > zMax) ? pattern[p].z : zMax;
  }
  xMin -= 1e-4*h1;								// Extends bbox to ensure overlap and real check.
  yMin -= 1e-4*h2;
  zMin -= 1e-4*h3;
  xMax += 1e-4*h1;
  yMax += 1e-4*h2;
  zMax += 1e-4*h3;

  double nPart = 0;
  for (int i = dmn_min[0] ; i < dmn_max[0] + 1 - mc_have_x ; i += Nx)	// Fills domain by pattern.
  {
    for (int j = dmn_min[1] ; j < dmn_max[1] + 1 - mc_have_y ; j += Ny)
    {
      for (int k = dmn_min[2] ; k < dmn_max[2] + 1 - mc_have_z ; k += Nz)
      {
        if ( ((i*h1 + xMax < cpu_min[0]*h1 || i*h1 + xMin > cpu_max[0]*h1) && mc_have_x) ||		// Comparison by bounding box.
             ((j*h2 + yMax < cpu_min[1]*h2 || j*h2 + yMin > cpu_max[1]*h2) && mc_have_y) ||
             ((k*h3 + zMax < cpu_min[2]*h3 || k*h3 + zMin > cpu_max[2]*h3) && mc_have_z) )
          continue;

        for (int l = 0 ; l < patternSize ; ++l)
        {
          double x, y, z;
          x = pattern[l].x + i*h1;
          y = pattern[l].y + j*h2;
          z = pattern[l].z + k*h3;

          if (((x >= cpu_min[0]*h1 && x < cpu_max[0]*h1) || !mc_have_x) && 		// Fine filter for || setup.
              ((y >= cpu_min[1]*h2 && y < cpu_max[1]*h2) || !mc_have_y) &&
              ((z >= cpu_min[2]*h3 && z < cpu_max[2]*h3) || !mc_have_z))
          {
            marker_t *marker = plasma_marker ();				// Gets pointer on uninitialized marker.
            marker->x = x;
            marker->y = y;
            marker->z = z;
            marker->vx = pattern[l].vx;
            marker->vy = pattern[l].vy;
            marker->vz = pattern[l].vz;
            marker->rho = pattern[l].rho;
            nPart++;
          }
        }
      }
    }
  }

  free (pattern);								// Releases memory used to store stencil.

  return nPart;
}
