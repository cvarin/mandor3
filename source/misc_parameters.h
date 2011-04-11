/** \file misc_parameters.h
  * See misc_parameters.c for full description.
  *
  * Few parameters are always 'must have': domain size, boundary conditions,
  * number of cpu the code operates on, etc. This module holds global
  * parameters (cpu_here, cpu_total, etc) which are protected by 'const'
  * qualifier. Initialization is done using simple interface to pass data
  * behind the const-protection barrier.
  *
  * \attention One may use many ways to scan all boundariers and sometimes it
  * is desirable to have all 6 parameters (6 boundary conditions, for example)
  * packed into single array. It is done already - feel free to use dmn_min
  * as an array with 6 elements {imin, jmin, kmin, imax, jmax, kmax} safely –
  * it is allocated properly to ensure it. The same is true for dmn_bc_min.
  *
  * \warning External reference for tau, h1, ..., are 'extern const', so if you
  * update parameter, be careful - sometimes compiler can overoptimize it. I
  * did set 'volatile' qualifier on time but I want to keep free hands to
  * optimize the core of the engine: global data is set at start-up and is
  * never changed so I do NOT use 'volatile' qualifier on mesh size or stuff
  * like that.
  *
  * UPDATE: gcc with aggressive optimization didn't respect 'volatile' flag so
  * declaration for Time is changed back to a traditional one.
  */

#ifndef MC_MISC_PARAMETERS_HEADER
#define MC_MISC_PARAMETERS_HEADER

/// Maximal theoretical number of the node (simular to INT_MAX).
/// Used to scan files created by many nodes without boring tests, etc.
#define CPU_MAX (10000)

#include "frame.h"

/// Enumerated boundary conditions (last element is used to specify array
/// mappings 'BC -> something').
enum { BC_PERIODIC,
       BC_MIRROR,
       BC_OPEN,
       BC_SPLITTER,
       BC_ENUM_LENGHT };

void parameter_dump (void);
void parameter_save (void);
void parameter_load (void);

void parameter_enterMPI    (int argc, char *argv[], int continue_log);
void parameter_setupMesh   (int i, int j, int k,
                            double X, double Y, double Z, double TAU);
void parameter_setupBounds (int xMin, int xMax, int yMin, int yMax, int zMin, int zMax);

void parameter_setTime     (double time);

#ifndef MC_MISC_PARAMETER_SOURCE
   extern const int    cpu_here, cpu_total;
   extern const double tau, h1, h2, h3, Lx, Ly, Lz;

   // Global parameters of the computational domain - size of the mesh (no
   // ghost cells included!) and boundary conditions.
   extern const int *dmn_min,
                    *dmn_max,
                    *dmn_bc_min,
                    *dmn_bc_max;
#endif

double Time;

/// Packed state of the axises: presented or removed at compile-time (I keep it
/// local to help in optimization).
static const int ACTIVATOR[6] = { mc_have_x, mc_have_y, mc_have_z,
                                  mc_have_x, mc_have_y, mc_have_z };

#endif
