/** \file misc_parameters.c
  * Parameters manager.
  */

#define MC_MISC_PARAMETER_SOURCE	///< Removes 'const' protection from h1, .., dmn_bc_max (see 'misc_parameters.h').

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "log.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

/// Global parameters of the simpulation (initialized as full crap to ensure explicit initialization).
int    cpu_here  = -1,
       cpu_total = -1;

double h1   = NAN,
       h2   = NAN,
       h3   = NAN,
       Lx   = NAN,
       Ly   = NAN,
       Lz   = NAN,
       tau  = NAN,
       Time = 0;

static int dmn_bc[6]   = {-1, -1, -1, -1, -1, -1};
static int dmn_mesh[6] = { 0,  0,  0, -1, -1, -1};

int *dmn_bc_min = dmn_bc,
    *dmn_bc_max = dmn_bc + 3;

int *dmn_min    = dmn_mesh,
    *dmn_max    = dmn_mesh + 3;

// ----------------------------------------------------------------------------
/// Broadcasts parameters (master node updates all other nodes).
// ----------------------------------------------------------------------------
static void
parameter_bcast (void)
{
   struct tmp_s {
      double h1, h2, h3, tau;
      double Lx, Ly, Lz;
   } package = { h1, h2, h3, tau, Lx, Ly, Lz };

   MPI_Bcast (&package,   sizeof (package), MPI_BYTE, 0, MPI_COMM_WORLD);
   MPI_Bcast (dmn_bc_min, 6,                MPI_INT,  0, MPI_COMM_WORLD);
   MPI_Bcast (dmn_min,    6,                MPI_INT,  0, MPI_COMM_WORLD);

   h1  = package.h1;
   h2  = package.h2;
   h3  = package.h3;
   Lx  = package.Lx;
   Ly  = package.Ly;
   Lz  = package.Lz;
   tau = package.tau;
}

// ----------------------------------------------------------------------------
/// Updates time (used for main engine timestep only).
// ----------------------------------------------------------------------------
void
parameter_setTime (double time)
{
   Time = time;
}

// ----------------------------------------------------------------------------
/// Reports global envoronment.
// ----------------------------------------------------------------------------
void
parameter_dump (void)
{
   char states[2][20] = {"deactivated", "active"};
   char BC    [3][20] = {"periodic", "mirror", "open"};

   say ("global parameters:");
   say ("  - time step: %e", tau);
   say ("  - mesh steps: %e, %e, %e", h1, h2, h3);
   // XXX Use reg-pack + reg_print.
   say ("  - mesh size: [%d, %d] x [%d, %d] x [%d, %d]", dmn_min[0], dmn_max[0],
                                                         dmn_min[1], dmn_max[1],
                                                         dmn_min[2], dmn_max[2]);
   say ("  - domain size: %f x %f x %f", Lx, Ly, Lz);
   say ("  - boundary conditions:");
   say ("    o X[%s, %s], %s", BC[dmn_bc_min[0]], BC[dmn_bc_max[0]], states[mc_have_x]);
   say ("    o Y[%s, %s], %s", BC[dmn_bc_min[1]], BC[dmn_bc_min[1]], states[mc_have_y]);
   say ("    o Z[%s, %s], %s", BC[dmn_bc_min[2]], BC[dmn_bc_min[2]], states[mc_have_z]);
}

// ----------------------------------------------------------------------------
/// Saves parameters.
// ----------------------------------------------------------------------------
void
parameter_save (void)
{
   // Only master saves the data.
   if (cpu_here)
      return;

   FILE *fp = cfg_open ("binData/domain.bin", "wb", __func__);
   double doubles[7] = { h1, h2, h3, Lx, Ly, Lz, tau };
   fwrite (doubles,    sizeof (double), 7, fp);
   fwrite (dmn_min,    sizeof (int),    6, fp);
   fwrite (dmn_bc_min, sizeof (int),    6, fp);
   fclose (fp);
}

// ----------------------------------------------------------------------------
/// Loads parameters.
// ----------------------------------------------------------------------------
void
parameter_load (void)
{
   FILE *fp = cfg_open ("binData/domain.bin", "rb", __func__);
   fread (&h1,  sizeof (double), 1, fp);
   fread (&h2,  sizeof (double), 1, fp);
   fread (&h3,  sizeof (double), 1, fp);
   fread (&Lx,  sizeof (double), 1, fp);
   fread (&Ly,  sizeof (double), 1, fp);
   fread (&Lz,  sizeof (double), 1, fp);
   fread (&tau, sizeof (double), 1, fp);
   fread (dmn_min,    sizeof (int), 6, fp);
   fread (dmn_bc_min, sizeof (int), 6, fp);
   fclose (fp);

   parameter_bcast ();
}


// ----------------------------------------------------------------------------
/// Wrapped 'MPI_Finalize' (for 'atexit').
// ----------------------------------------------------------------------------
static void
parameter_exitMPI (void)
{
   ENSURE (MPI_Finalize () == MPI_SUCCESS, "cannot finish MPI session");
   SAY_DEBUG ("MPI (cpu #%d): session closed", cpu_here);
}

// ----------------------------------------------------------------------------
/// Opens MPI session, initializes 'cpu_here', 'cpu_total', logger, etc.
// ----------------------------------------------------------------------------
void
parameter_enterMPI (int argc, char *argv[], int continue_log)
{
   static int MPI_is_not_set = 1;
   ENSURE (MPI_is_not_set, "double initializaton");

   MPI_Init      (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &cpu_total);
   MPI_Comm_rank (MPI_COMM_WORLD, &cpu_here);

   ENSURE (cpu_total < CPU_MAX, "CPU_MAX is too small, fix 'misc_parameters.h'"
                                " (we need at least %d instead of %d)",
				cpu_total + 1, CPU_MAX);

   int  name_lenght;
   char processorName[MPI_MAX_PROCESSOR_NAME];
   MPI_Get_processor_name (processorName, &name_lenght);

   // Finds name of the executable (strips all path info).
   char *exec_name = argv[0];
   while (strstr (exec_name, "/")) {
      exec_name = strstr (exec_name, "/") + 1;
   }

   // Opens log file with exec-file/cpu encoded in.
   log_open (cpu_here, continue_log, _("output/logs/%s_%03d.log", exec_name,
                                                                  cpu_here));
   say ("MPI: %d cpu(s) total, node number is %d (%s)",
        cpu_total, cpu_here, processorName);
   ENSURE (!atexit (parameter_exitMPI), "cannot submit 'parameter_exitMPI'");

   MPI_is_not_set = 0;
}

// ----------------------------------------------------------------------------
/// Initialization of domain size, mesh size, spatial and temporal steps.
// ----------------------------------------------------------------------------
void
parameter_setupMesh (int i, int j, int k,
                     double X, double Y, double Z, double TAU)
{
   ENSURE (mc_have_x == (i > 1) && mc_have_y == (j > 1) && mc_have_z == (k > 1),
           "badly compiled code; axises are %d%d%d, domain is %d x %d x %d",
           mc_have_x, mc_have_y, mc_have_z, i, j, k);

   ENSURE (i >= 0 && j >= 0 && k >= 0 && X > 0 && Y > 0 && Z > 0,
           "bad domain size (%.3e x %.3e x %.3e) or mesh size (%d x %d x %d)",
           X , Y, Z, i, j, k);

   // Updates mesh/domain sizes.
   int tmp[6] = { 0, 0, 0, i, j, k };
   memcpy (dmn_min, tmp, 6*sizeof (int));
   Lx = X;
   Ly = Y;
   Lz = Z;

   // Updates mesh steps.
   h1 = X / (double) (i + (i == 0));
   h2 = Y / (double) (j + (j == 0));
   h3 = Z / (double) (k + (k == 0));
   tau = TAU;

   parameter_bcast ();
}

// ----------------------------------------------------------------------------
/// Sets boundary conditions.
// ----------------------------------------------------------------------------
void
parameter_setupBounds (int xMin, int xMax, int yMin, int yMax, int zMin, int zMax)
{
   // Makes sure degenerated axises have periodic BC.
   xMin *= mc_have_x;
   xMax *= mc_have_x;
   yMin *= mc_have_y;
   yMax *= mc_have_y;
   zMin *= mc_have_z;
   zMax *= mc_have_z;

   // Tests that global BC are 'physical'.
   ENSURE (xMin >= 0 && xMin <= BC_OPEN && xMax >= 0 && xMax <= BC_OPEN &&
           yMin >= 0 && yMin <= BC_OPEN && yMax >= 0 && yMax <= BC_OPEN &&
           zMin >= 0 && zMin <= BC_OPEN && zMax >= 0 && zMax <= BC_OPEN,
           "bad boundary condition: X(%d/%d), Y(%d/%d), Z(%d/%d)",
           xMin, xMax, yMin, yMax, zMin, zMax);

   ENSURE ((xMin == BC_PERIODIC) == (xMax == BC_PERIODIC),
           "X axis: periodic BC should be set on BOTH boundaries");
   ENSURE ((yMin == BC_PERIODIC) == (yMax == BC_PERIODIC),
           "Y axis: periodic BC should be set on BOTH boundaries");
   ENSURE ((zMin == BC_PERIODIC) == (zMax == BC_PERIODIC),
           "Z axis: periodic BC should be set on BOTH boundaries");

   // Updates boundary conditions.
   int tmp[6] = { xMin, yMin, zMin, xMax, yMax, zMax };
   memcpy (dmn_bc_min, tmp, 6*sizeof (int));

   parameter_bcast ();
}
