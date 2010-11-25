/** \file log_test.c
  * \brief Simple demo of the logging subsystem (see 'log.c').
  */

#include <mpi.h>
#include <unistd.h>

#include "log.h"

// ---------------------------------------------------------------------------
/// Autoclosing of the log in exit (see 'log.c' for docs).
// ---------------------------------------------------------------------------
static void
exit_MPI (void)
{
   ENSURE (MPI_Finalize () == MPI_SUCCESS, "cannot finish MPI session");
}


int
main (int argc, char **argv)
{
   int cpu_here, cpu_total;

   MPI_Init (&argc, &argv);
   MPI_Comm_size (MPI_COMM_WORLD, &cpu_total);
   MPI_Comm_rank (MPI_COMM_WORLD, &cpu_here);

   int  nameLenght;
   char processorName[MPI_MAX_PROCESSOR_NAME];
   MPI_Get_processor_name (processorName, &nameLenght);

   ENSURE (atexit (exit_MPI) == 0, "cannot submit 'MPI_Finalize' wrapper");

   log_open (cpu_here, 0, _("output/logs/logtest_%02d.log", cpu_here));

   MPI_Barrier (MPI_COMM_WORLD);

   say ("Node: %d", cpu_here);
   SAY_WARNING ("that is only a demo --- actual sources can show a real fun");
   SAY_DEBUG ("node %d out of %d", cpu_here, cpu_total);

   for (double i = 0 ; i < 2e6 ; ++i) {
      say_doing ("having fun: %.2f ...", i);
//       sleep (1);
   }

   int I_have = 2,
       I_need = 8;
   SAY_OMG ("scary error");
   ENSURE  (I_have >= I_need, "no fuel: have %d, need %d", I_need, I_have);

   return EXIT_SUCCESS;
}
