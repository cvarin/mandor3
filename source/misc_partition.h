/** \file misc_partition.h
  * \brief Prototypes and predefined constants for the partitioning unit.
  *
  * <h3>Intro.</h3>
  *
  * This module maintains multilevel tree of domain subdivisions. If domain is
  * splitted on N = p1*p2*..*pL cpus one can subdivide domain on p1 groups,
  * than divide each group on p2 subgroups and so on. Independecy of subdividions
  * introduces enough free parameters to perfectly balance the tree for any
  * distribution of work load over cluster.
  *
  * <h3>Command line arguments.</h3>
  *
  * Module can read options in form of <b>-p:argString</b> or <b>--partition:argString</b>.
  * \b argString contains no spaces and can be sequence of this commands:
  * - \b -x or \b -X: disables splitting of X-axis
  * - \b -y or \b -Y: disables splitting of Y-axis
  * - \b -z or \b -Z: disables splitting of Z-axis
  * - <b>xnum</b> or <b>Xnum</b>: direct input of partitioning sequence (\b num slices across X axis)
  * - <b>ynum</b> or <b>Ynum</b>: direct input of partitioning sequence (\b num slices across Y axis)
  * - <b>znum</b> or <b>Znum</b>: direct input of partitioning sequence (\b num slices across Z axis)
  * - <b>~x[num]</b> or <b>~X[num]</b>: displaces partitioning planes randomly with
  *   amplitude equal to \b num (if no \b num is given default value misc_partition.c::MC_SKEW_DEFAULT is used)
  * - <b>~y[num]</b> or <b>~Y[num]</b>: jitter for Y-axis
  * - <b>~z[num]</b> or <b>~Z[num]</b>: jitter for Z-axis
  *
  * Purpose of this command line interface is to provide easy testing facility and
  * full control over partitioning to an independent caller. User can even directly
  * tell partitioning he wants.
  *
  * Examples:
  * - partitioning along Z axis only:
  *   - mpirun -np 3 ./core.out -p:-y-z
  *   - mpirun -np 3 ./core.out -p:-Y --partition:-X
  * - partitioning along XY axises only with random shift in [-5,5] for X and [-10, 10] for Y:
  *   - mpirun -np 20 ./core.out -p:-z~x5~Y10
  *   - mpirun -np 20 ./core.out --partition:-Z --partition:~x5 -p:~Y10
  * - 3 slices across X axises and 2 slices across Z:
  *   - mpirun -np 6 ./core.out -p:x3z2
  *   - mpirun -np 6 ./core.out --partition:X3 --partition:Z2
  */

#ifndef MC_MISC_PARTITION_HEADER
#define MC_MISC_PARTITION_HEADER					///< Include guard.

#include "type_reg.h"

#include "misc_parameters.h"

#ifndef MC_MISC_PARALLEL_SOURCE
  /// Protected declaration of local (current node) domain sizes and local boundary conditions.
  extern const int *cpu_min, *cpu_max;
  extern const int *cpu_bc_min, *cpu_bc_max;
#endif

void   partition_init (void);
void   partition_show (void);
void   partition_save (const char *name);
void   partition_load (const char *name);
regList_t partition_getMapping (void);

#endif
