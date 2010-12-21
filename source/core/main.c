/** \file main.c
  * Numeric core, runs the sumulation and saves binary data for analysis.
  *
  * Engine has few setting you can set in Makefile (nicely commented macro-keys) and
  * few parameters are read from config files at startup ('run_mandor.cfg' and
  * 'run_probes.cfg') and from periodically during a run ('tmp/stop.flag').
  *
  * You can stop the code at run-time by changing file 'tmp/stop.flag' - see
  * main_throwStopFlag ().
  */

#define mc_activateProfiler 	1	///< Activates build-in profiler.

#include <stdlib.h>
#include <mpi.h>

#include "timer.h"
#include "type_reg.h"
#include "type_marker.h"

#include "log.h"

#include "test_mesh.h"

#include "misc_units.h"
#include "misc_MPItags.h"
#include "misc_cfgReader.h"
#include "misc_definedKeys.h"

#include "IO_sys.h"
#include "IO_tecplot.h"
#include "IO_fileMap.h"

#include "profiler.h"
#include "spectr_dump.h"
#include "commandLine.h"

#include "core/em.h"
#include "core/plasma.h"
#include "core/plasma_IO.h"
#include "core/plasma_VSP.h"
#include "core/plasma_parallel.h"
#include "core/plasma_gaussTest.h"

#include "core/diag_probe.h"
#include "core/diag_WDensity.h"

// TEMP
static double tmp_gamma[100];
const double tmp_max_gamma = 2.0;

void tmp_update_gamma(int start, int finish) {
  for (int i = start; i < finish; i++) if (plasma[i].qDivM < 0) {
    marker_t *p = &plasma[i];
    const double gamma = sqrt(1.0 + p->vx*p->vx + p->vy*p->vy + p->vz*p->vz);
    const j = (int)(100.0 * (gamma - 1.0) / (tmp_max_gamma - 1.0));
    if (j < 100) tmp_gamma[j] -= p->rho;
  }
}

void tmp_save_gamma(int record, int cpu) {
  char name[256];
  sprintf(name, "binData/gamma_%06d_%03d.bin", record, cpu);

  FILE *fp = fopen(name, "wt");
  ENSURE(fp, "Cannot open file '%s' for writing.", name);

  for (int i = 0; i < 100; i++) {
    fprintf(fp, "%lg %lg\n", 0.511 * double(i) * (tmp_max_gamma - 1.0) / 100.0, tmp_gamma[i]);
    tmp_gamma[i] = 0.0;
  }
}

// END TEMP


// Variables for non-blocking check of the stop requests/timeouts.
static timeTick_t  stopFlag_startTime;		///< Simulation start time (to check timeouts).
static double      stopFlag_saveTime = 0;	///< Simulation checkpoint saving time (to check timeouts).
static MPI_Request stopFlag_requests[2000];	///< Request for stop flag state transfer.

/// Check-pointing steps and simulation length (in timesteps).
static int totalSteps          =  0,
           chPoint_tecplot     = -1,
           chPoint_full        = 10,
           chPoint_stopRequest = 20,
           chPoint_spectr      = -1;

static meshDouble_t rho;	///< Charge density ρ.
static meshVec_t    E, H, J;	///< Vector fields \f$ (\vec E,\ \vec H,\ \vec J). \f$

// ---------------------------------------------------------------------------
/// Allocates dummy meshes to initialize all "type" fields and set the names.
// ---------------------------------------------------------------------------
static void
main_allocateMeshes (void)
{
   mesh_allocate (mcast_mesh(&E),   0, 0, 0, 0, 0, 0, "vec_E",    mc_vec3D_t);
   mesh_allocate (mcast_mesh(&H),   0, 0, 0, 0, 0, 0, "vec_H",    mc_vec3D_t);
   mesh_allocate (mcast_mesh(&J),   0, 0, 0, 0, 0, 0, "vec_J",    mc_vec3D_t);
   mesh_allocate (mcast_mesh(&rho), 0, 0, 0, 0, 0, 0, "scal_rho", mc_double);
}

// ---------------------------------------------------------------------------
/// Reconfigures mesh sizes to adjust it to the local node geometry.
// ---------------------------------------------------------------------------
static void
main_reconfigureMeshes (void)
{
   mesh_resize (mcast_mesh(&E),   cpu_min[0] - 1, cpu_min[1] - 1, cpu_min[2] - 1, cpu_max[0] + 1, cpu_max[1] + 1, cpu_max[2] + 1);
   mesh_resize (mcast_mesh(&H),   cpu_min[0] - 1, cpu_min[1] - 1, cpu_min[2] - 1, cpu_max[0] + 1, cpu_max[1] + 1, cpu_max[2] + 1);
   mesh_resize (mcast_mesh(&J),   cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2, cpu_max[0] + 2, cpu_max[1] + 2, cpu_max[2] + 2);
   mesh_resize (mcast_mesh(&rho), cpu_min[0] - 2, cpu_min[1] - 2, cpu_min[2] - 2, cpu_max[0] + 2, cpu_max[1] + 2, cpu_max[2] + 2);
}

// ---------------------------------------------------------------------------
/// Checks external termination request or time-out (decision is made <b>by master node only</b> thus everything is coherent).
// ---------------------------------------------------------------------------
static void
main_throwStopFlag (void)
{
   // Workers decide nothing.
   if (cpu_here)
      return;

   static int buffer[2];
   FILE *fp = cfg_open ("tmp/stop.flag", "rt", __func__);
   // "Fire alarm" and time limit in seconds.
   fscanf (fp, "%d %d", buffer, buffer + 1);
   fclose (fp);

   static timeTick_t now,
                     prevCheckTime;
   time_get (&now);

   // Checks if next check-point time < walltime.
   static double dT = 0;
   double estimate = time_elapsed (&stopFlag_startTime, &now)
                   + 1.2*dT + 1.5*stopFlag_saveTime;
   if (estimate >= buffer[1] || buffer[0]) {
      say ("main_throwStopFlag: time limit reached\n"
           "  o time = %f\n  o save time = %f\n"
           "  o interval between checks = %f\n"
           "  o estimate next hit time = %f\n  o time limit = %d\n",
           time_elapsed(&stopFlag_startTime, &now),
           stopFlag_saveTime,
           dT, estimate, buffer[1]);
      buffer[0] = 1;
   }

   // Updates estimate of step duration.
   static int first = 1;
   if (!first)
      dT = 0.8*time_elapsed (&prevCheckTime, &now) + 0.2*dT;
   prevCheckTime = now;

   ENSURE (cpu_total <= 1000,
           "increase 'stopFlag_requests' array at least up to %d", cpu_total);

   // Sends stop flag to everybody.
   buffer[1] = cpu_total;
   for (int cpu = 0 ; cpu < cpu_total ; cpu++) {
      MPI_Isend (buffer, 2, MPI_INT, cpu, TAG_STOP_REQUEST, MPI_COMM_WORLD,
                 stopFlag_requests + cpu);
   }
}

// ---------------------------------------------------------------------------
/// Checks external termination request (introduces soft barrier: nobody will
/// start new time step without permission from master, issued by
/// main_throwStopFlag()).
// ---------------------------------------------------------------------------
static int
main_catchStopFlag (void)
{
   // Receives stop flag.
   int buffer[2];
   MPI_Recv (buffer, 2, MPI_INT, 0, TAG_STOP_REQUEST, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
   ENSURE   (buffer[1] == cpu_total,
             "number of cpus is changed (%ld -> %d)", buffer[1], cpu_total);

   // Master waits for all Isends to complete.
   if (!cpu_here)
      MPI_Waitall (cpu_total, stopFlag_requests, MPI_STATUSES_IGNORE);

   return buffer[0];
}

// ---------------------------------------------------------------------------
/// Reads <b>run_mandor.cfg</b> file and sets all global static parameters (called from main()).
// ---------------------------------------------------------------------------
static void
main_startUp (void)
{
   // Loading time is a first estimate of a saving time.
   timeTick_t startUpStart;
   time_get (&startUpStart);

   main_allocateMeshes ();

   FILE *fp;
   int buffer[6];
   if (!cpu_here) {
      fp = cfg_open ("run_mandor.cfg", "rt", __func__);	// Reads checkpointing steps and startup flags.
      buffer[0] = cfg_readInt (fp);			//  - totalSteps.
      buffer[1] = cfg_readInt (fp);			//  - chPoint_full.
      buffer[2] = cfg_readInt (fp);			//  - chPoint_tecplot.
      buffer[3] = (1 == cfg_readInt (fp));		//  - continueRun.
      buffer[4] = cfg_readInt (fp);			//  - chPoint_stopRequest
      buffer[5] = cfg_readInt (fp);			//  - chPoint_spectr
      fclose (fp);
   }

   // Broadcasts parameters.
   MPI_Bcast (buffer, 6, MPI_INT, 0, MPI_COMM_WORLD);
   totalSteps          = buffer[0];
   chPoint_full        = buffer[1];
   chPoint_tecplot     = buffer[2];
   int continueRun     = buffer[3];
   chPoint_stopRequest = buffer[4];
   chPoint_spectr      = buffer[5];

   // Checks if sizes can be divided by.
   ENSURE (totalSteps*chPoint_full*chPoint_tecplot
                     *chPoint_stopRequest*chPoint_spectr,
           "'run_mandor.cfg' contains zero step(s)");

   // Sets time barrier.
   if (!(fp = fopen ("tmp/stop.flag", "rt"))) {
      fp = cfg_open ("tmp/stop.flag", "wt", __func__);
      // Sets timeout.
      SAY_DEBUG ("Setting timeout to 23h30m.");
      fprintf (fp, "0\n%d\n", 23*60*60 + 30*60);
      fclose  (fp);
   }

   // Record are numbered as 0, .., N - 1.
   int point = sysIO_recordsTotal () - 1;
   MPI_Bcast (&point, 1, MPI_INT, 0, MPI_COMM_WORLD);

   // Chooses start-up mode.
   ENSURE (point >= 0, "initial checkpoint is not prepared");

   if (continueRun && point != 0) {
      fileMap_init (mc_fileMap_contn);
      sysIO_setRecordNum  (-1);		// Continues check-pointing.
      tecIO_setRecordNum  (-1);		// Continues tecplot output.
      spectr_continueDump ();		// Continues spectr output.
   } else {
      fileMap_init (mc_fileMap_start);
      sysIO_setRecordNum (1);		// Starts to write just after setup record.
      tecIO_setRecordNum (1);		// Starts to write just after setup record.
      spectr_startDump   ();		// Starts new spectr output sequence.
      point = 0;
   }

   // Initializes total energy diagnostic.
   wDensity_prepare (!(continueRun && point != 0));

   /// \todo Partitioning options outside.
   SAY_WARNING("Explicit partitioning is blocked.");
   if (/*continueRun && point != 0 &&*/ 0) {
      partition_load ("output/partition_onStart.cfg");	// Loads old partitioning.
   } else {
      partition_init ();				// Initializes partitioning.
      partition_save ("output/partition_onStart.cfg");	// Saves partitioning used.
   }

   // Sets the size of meshes using the partitioning.
   main_reconfigureMeshes ();

   double time;
   const reg_t to_load = {{cpu_min[0], cpu_min[1], cpu_min[2]},
                          {cpu_max[0], cpu_max[1], cpu_max[2]}, 0};
   sysIO_loadEM (point, &time, &to_load, &E, &H);
   plasma_load  (point);

   mf_mesh_clean (&J);
   mf_mesh_clean (&rho);

   parameter_setTime (time);

   parameter_dump ();		// Prints parameters loaded.
   partition_show ();		// Tecplot visualization of partititon.

   // Allocates probes and opens associated files.
   probe_allocate (continueRun && point > 0);

   say ("main_startUp: config file is imported. Parameters of the run are:");
   say ("  - %d total steps",                totalSteps);
   say ("  - %d steps between check-points", chPoint_full);
   if (chPoint_tecplot > 0) {
      say ("  - %d steps between tecplot shots", chPoint_tecplot);
   } else {
      say ("  - tecplot shots are off");
   }
   if (chPoint_spectr > 0) {
      say ("  - %d steps between spectral data dumps", chPoint_spectr);
   } else {
      say ("  - EM spectral energy distribution diagnostic is off");
   }
   say ("  - %d steps between updating a stop request", chPoint_stopRequest);
   say ("  - simulation is %s from check-point #%03d.",
                                    (continueRun && point != 0) ? "continued"
                                                                : "started",
                                    point);

   // Loading time is a good estimate of saving time.
   timeTick_t startUpFinish;
   time_get (&startUpFinish);
   stopFlag_saveTime = time_elapsed (&startUpStart, &startUpFinish);

   // Configures convolutioner.
   double VSP_weight = 1.0;
   if (cl_findDouble ("VSP:weight", &VSP_weight))
      VSP_weight = (VSP_weight > 0.5 && VSP_weight < 1.00001) ? VSP_weight
                                                              : 1.00;
   VSP_configure (VSP_weight);
}

// ---------------------------------------------------------------------------
/// Mandor2 core.
// ---------------------------------------------------------------------------
int
main (int argc, char **argv)
{
   time_get (&stopFlag_startTime);
   cl_import ("./source/import.pl", argc, argv, 0);

   parameter_enterMPI (argc, argv, 1);
   parameter_load ();

   // Initializes profiler.
   profiler_init (_("output/prof_core_%03d.txt", cpu_here));

   units_load ();

   // Configures all modules.
   main_startUp ();
   em_init      (&E, &H);
   plasma_init  ();

   for (int tmp = 0; tmp < 100; tmp++) tmp_gamma[tmp] = 0.0; // TEMP

   say ("System is initialized successfully.");

   // Total time of the work stage - starts here.
   int    stopFlag   = 0;
   double computTime = MPI_Wtime (),
          oldTime    = -1;
   for (int t = 1 ; t <= totalSteps && !stopFlag ; t++) {
      double W_E, W_M, WTx = 0, WTy = 0, WTz = 0;
      profiler_startLoop ();

      // XXX Add it to the timing counter.
      comm_plasma_start ();

      // Reads data from external file (works as mailbox)..
      profiler_begin (mc_prof_throwStopFlag);
      if (chPoint_stopRequest > 0 && t%chPoint_stopRequest == 0)
         main_throwStopFlag ();

      // H^(n-1/2) -> H^n.
      profiler_endBegin (mc_prof_emh1);
      em_HHalfStep (mcast_meshVec_RO (&E), &H);

      // Energy density of the field at t = n*tau.
      profiler_endBegin (mc_prof_EMenergy);
      em_energy (mcast_meshVec_RO(&E), mcast_meshVec_RO(&H), &W_E, &W_M);
      profiler_endBegin (mc_prof_probes);
      probe_postData (Time, mcast_meshVec_RO(&E), mcast_meshVec_RO(&H));

      // Gauss law test beginning (profiler call is inside).
      gauss_before ();

      profiler_endBegin (mc_prof_plasma_move);
      plasma_move (mcast_meshVec_RO(&E), mcast_meshVec_RO(&H), &J);

      // Energy density: gets plasma termal energy.
      profiler_endBegin (mc_prof_plasma_Txyz);
      plasma_temperature (&WTx, &WTy, &WTz);

      // Tecplot's diagnostic check-pointing.
      if (chPoint_tecplot > 0 && t%chPoint_tecplot == 0) {
         say_doing ("writing tecplot mesh...");
         profiler_endBegin (mc_prof_tecRho);
         plasma_rho (&rho);
         profiler_endBegin (mc_prof_tecEH);
         tecIO_saveFields (Time, mcast_meshVec_RO(&E), mcast_meshVec_RO(&H));
      }

      // Spectral energy density diagnostic.
      if (chPoint_spectr > 0 && t%chPoint_spectr == 0) {
         say_doing ("writing spectral energy density dump...");
         profiler_endBegin (mc_prof_spectr);
         reg_t reg = {{cpu_min[0], cpu_min[1], cpu_min[2]},
                      {cpu_max[0] - mc_have_x, cpu_max[1] - mc_have_y, cpu_max[2] - mc_have_z}};
         reg_t map = {{E.imin, E.jmin, E.kmin},
                      {E.imax, E.jmax, E.kmax}};
         spectr_dump (Time, cpu_here, cpu_total, &reg, &map, E.storage, H.storage);
      }

      // Saves local energies to the buffer.
      profiler_endBegin (mc_prof_wDensity);
      wDensity_addPoint (W_E, W_M, WTx, WTy, WTz);

      // H^n -> H^(n+1/2).
      profiler_endBegin (mc_prof_emh2);
      em_HStep (mcast_meshVec_RO(&E), &H);

      // E^n -> E^(n+1).
      profiler_endBegin (mc_prof_eme);
      em_EStep_start (&E, mcast_meshVec_RO(&H));

      // Gets full picture of the currents.
      profiler_endBegin (mc_prof_jbcFinish);
      jbc_finish (&J);

      profiler_endBegin (mc_prof_plasma_pbcRecv);
      comm_plasma_complete (1000);

      profiler_endBegin (mc_prof_EFinish);
      em_EStep_finish (&E, &H, mcast_meshVec_RO(&J));

      // Gauss law test ending.
      gauss_after (mcast_meshVec_RO(&J), mcast_meshVec_RO(&E));

      // Tecplot's diagnostic check-pointing.
      if (chPoint_tecplot > 0 && t%chPoint_tecplot == 0) {
         say_doing ("writing tecplot mesh...");
         profiler_endBegin (mc_prof_tecRhoJ);
         tecIO_saveCurrents (mcast_meshVec_RO(&J), mcast_meshDouble_RO(&rho));
      }

      // Updates time after successful time step.
      parameter_setTime (Time + tau);

      tmp_update_gamma(countCore, countAll); // TEMP
      tmp_update_gamma(0, countShell); // TEMP

      // System check-point.
      if (t && t%chPoint_full == 0) {
        tmp_save_gamma(sysNum, cpu_here); // TEMP

         timeTick_t startWrite;
         time_get (&startWrite);
         say_doing ("writing check-point...");

         profiler_endBegin (mc_prof_sysSave);
         plasma_save (sysNum, cpu_here);
         sysIO_save  (Time, mcast_meshVec_RO(&E), mcast_meshVec_RO(&H));

         profiler_endBegin (mc_prof_wDensityFlush);
         // WARNING: barrier-type calls inside.
         wDensity_flushData ();
         oldTime = Time;

         // Updates estimate of saving time.
         timeTick_t endWrite;
         time_get (&endWrite);
         stopFlag_saveTime = time_elapsed (&startWrite, &endWrite);
      }

      profiler_endBegin (mc_prof_catchStop);

      // Gets data sent earlier.
      if (chPoint_stopRequest > 0 &&  t % chPoint_stopRequest == 0) {
         stopFlag = main_catchStopFlag ();
      }
      profiler_end ();

      say_doing ("time=%.4f (%.2f %%)   ", Time, 100.0*t/(double) totalSteps);
//       say ("step %d: time=%.4f (%.2f %%)", t, Time
//                                               , 100.0*t/(double) totalSteps);
      profiler_finishLoop ();
   }
   computTime = MPI_Wtime () - computTime;
   say ("main: main loop have taken %f sec.", computTime);

   // System check-point.
   if (oldTime != Time)	{
      say ("Addional check-point to save the final state.");
      plasma_save (sysNum, cpu_here);
      sysIO_save  (Time, mcast_meshVec_RO(&E), mcast_meshVec_RO(&H));
   }

   say ("Final barrier.");
   MPI_Barrier (MPI_COMM_WORLD);

   // Removes tmp file to signal the end of run.
   if (!cpu_here)
      remove ("tmp/stop.flag");

   return EXIT_SUCCESS;
}
