/** \file profiler.c
  * Run-time profiler (see profiler.h).
  */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "timer.h"

#define mc_blockDummyMacroDeclarations		///< To activate prototypes anyway.
#include "profiler.h"

#include "log.h"

#define mc_debugChecks	0			///< Define flag to activate array overflow/underflow run-time checks.
#define mc_eventCacheN  (100)			///< Size of the temporary array to store incoming events before flush.

#define mc_maxLevel     64			///< Maximum level allowed (see profLevelBeginnings).

// ---------------------------------------------------------------------------
/// Sampled interval (time of the start/end, level of the event (used to
/// calculated percentage compare to the higher level caller), id of caller.
// ---------------------------------------------------------------------------
typedef struct {
   timeTick_t startAt;	   ///< Time of the event's beginning.
   double     accumulated; ///< Accumulator: total time spend in the interval.
   int        id;	   ///< ID of the event.
   int        level;	   ///< Level (deepness).
} sample_t;

// ---------------------------------------------------------------------------
/// Descriptor of the profiler event (sampled interval).
// ---------------------------------------------------------------------------
typedef struct {
   int        nameLenght;		///< Lenght of the name of the event (used to speed-up formatted output).
   const char name[mc_prof_nameSize];	///< Name of the event.
} nameBind_t;

// ---------------------------------------------------------------------------
/// Registered events.
// ---------------------------------------------------------------------------
static nameBind_t profBinds[mc_prof_N] = {
   { 0, "main loop" 		},
   { 0, "HalfHStep" 		},
   { 0, "EStep" 		},
   { 0, "HStep" 		},
   { 0, "throwStopFlag" 	},
   { 0, "EM-energy" 		},
   { 0, "probes" 		},
   { 0, "divJ preparation" 	},
   { 0, "plasma/move" 		},
   { 0, "plasma/Txyz" 		},
   { 0, "tecIO/get rho" 	},
   { 0, "tecIO/write EH" 	},
   { 0, "spectr" 		},
   { 0, "wDensity cache" 	},
   { 0, "jbc finish" 		},
   { 0, "ebc finish" 		},
   { 0, "divJ testing" 		},
   { 0, "tecIO/write RhoJ" 	},
   { 0, "sysIO/save" 		},
   { 0, "wDensity/flush" 	},
   { 0, "defrag markers" 	},
   { 0, "catchStopFlag" 	},
   { 0, "TFSF/preH" 		},
   { 0, "cap/cacheEH" 		},
   { 0, "step" 			},
   { 0, "BC" 			},
   { 0, "cleanJ" 		},
   { 0, "push" 			},
   { 0, "pbc/pop" 		},
   { 0, "pbc/recv" 		},
   { 0, "jbc/send" 		},
   { 0, "recv mesh" 		},
   { 0, "recv markers" 		},
   { 0, "merge markers" 	},
   { 0, "VSP sweep" 		},
   { 0, "jbc/local bounds" 	},
   { 0, "ebc/cap"		},
   { 0, "ebc/recv"		},
   { 0, "ebc/jPass"		},
   { 0, "ebc/misc"		}
};

typedef sample_t *sample_p;

static int        profLevel       = -1;			///< Current level of inclusion.
static sample_p   profLevelBeginnings[mc_maxLevel];	///< Positions of the beginning of all opened blocks (for fast interval contraction).

static sample_t   profEvents[mc_eventCacheN];		///< Buffer to cache all event calls between main loop start/end.
static sample_p   profEvent       = profEvents;		///< Empty slot for the next event coming.
static sample_p   profAccumEvent  = 0;			///< Marker to point on the block in which accumulation is activated.

static int        profAccumLevel  = -1;
static sample_p   profAccumLevelBeginnings[mc_maxLevel];///< Positions of the beginning of all opened blocks (for fast interval contraction).

static FILE      *profFile        = NULL;		///< Output file stream.
static int        profLongestName = 10;			///< Lenght of the longest name (for evaluating the with of the columns for formatted output).
static timeTick_t profStartTime;			///< Moscow time :-))).

// ---------------------------------------------------------------------------
/// Closes output file.
// ---------------------------------------------------------------------------
static void
profiler_exit (void)
{
   if (profFile)
      fclose (profFile);
   profFile = NULL;
}

// ---------------------------------------------------------------------------
/// Gets initial reference point to measure time from (I assume that MPI_Wtime can be skewed).
// ---------------------------------------------------------------------------
void
profiler_init (const char *name)
{
   profFile = fopen (name, "wt");
   ENSURE (profFile, "cannot open output file stream");

   time_get (&profStartTime);
   double res, acc;
   res = time_resolution (&acc);
   fprintf (profFile,
            "Accuracy of the timer is %e sec +/- %e, start time = %e.\n",
            res, acc, time_seconds (&profStartTime));
   fprintf (profFile, "Name: (group %%) (execution time) (main loop %%)\n");

   // Sets initial level to "main loop".
   profLevel      = -1;
   profEvent      = profEvents;
   profAccumEvent = 0;

   ENSURE (!atexit (profiler_exit), "cannot submit deallocator");

   // Cleans all accumulated values.
   memset (profEvents, 0, mc_eventCacheN*sizeof (sample_t));

   // Updates precalcs for formatted output.
   profLongestName = 0;
   for (int bind = 0 ; bind < mc_prof_N ; ++bind) {
      profBinds[bind].nameLenght = strlen (profBinds[bind].name);
      profLongestName = (profLongestName < profBinds[bind].nameLenght) ? profBinds[bind].nameLenght
                                                                       : profLongestName;
   }
}

// ---------------------------------------------------------------------------
/// Gets initial reference point to measure time from (timers are assumed to be
/// skewed across the cluster).
// ---------------------------------------------------------------------------
void
profiler_syncClocks (void)
{
   time_get (&profStartTime);
}

// ---------------------------------------------------------------------------
/// Starts accumulation of the all events in the current sampling interval.
// ---------------------------------------------------------------------------
void
profiler_accumulate (void)
{
   // Remembers start of the summation region.
   if (!profAccumEvent) {
      profAccumEvent = profEvent - 1;
      memcpy (profAccumLevelBeginnings,
              profLevelBeginnings, sizeof (sample_t*)*(profLevel + 1));
   }
}

// ---------------------------------------------------------------------------
/// Registers the beginning of the block.
// ---------------------------------------------------------------------------
void
profiler_begin (int id)
{
#if defined (mc_debugChecks) && mc_debugChecks == 1
   ENSURE (profFile, "profiler is not initialized");

   ENSURE (profEvent < profEvents + mc_eventCacheN,
           "small cache, %d events already registered", profEvent - profEvents);
#endif

   if (profAccumEvent) {
      sample_t *event = profAccumEvent;
      while (event < profEvent && event->id != id) {
         ++event;
      }

      if (event < profEvent) {
         // Remembers start time and event type.
         time_get (&event->startAt);
         profAccumLevel = event->level;
         profAccumLevelBeginnings[profAccumLevel] = event;
         return;
      }
   }

   // Remembers start time and event type.
   time_get (&profEvent->startAt);
   profEvent->id    = id;
   profEvent->level = profAccumLevel = ++profLevel;
   // Remembers descriptor of the start of event.
   profAccumLevelBeginnings[profLevel] = profLevelBeginnings[profLevel]
                                       = profEvent;
   ++profEvent;

#if defined (mc_debugChecks) && mc_debugChecks == 1
   // Saves the "begin" to match with "end" later.
   ENSURE (profLevel >= 0 && profLevel < mc_maxLevel,
           "bad level (%d)", profLevel);
#endif
}

// ---------------------------------------------------------------------------
/// Registers the end of the block.
// ---------------------------------------------------------------------------
void
profiler_end (void)
{
#if defined (mc_debugChecks) && mc_debugChecks == 1
   ENSURE (profFile, "profiler is not initialized");

   ENSURE (profLevel >= 0 && profLevel < mc_maxLevel,
           "bad level (%d)", profLevel);
#endif

   // Remembers finalization time of event.
   timeTick_t endAt;
   time_get (&endAt);

   if (profAccumEvent) {
      sample_t *event = profAccumLevelBeginnings[profAccumLevel];
      event->accumulated += time_elapsed (&event->startAt, &endAt);
      if (event <= profAccumEvent)
         profAccumEvent = 0;

      if (profLevel == profAccumLevel
      &&  profAccumLevelBeginnings[profLevel] == profLevelBeginnings[profLevel]) {
         // New blocks should be closed as usual.
         --profLevel;
      }

    --profAccumLevel;
   } else {
      sample_t *event     = profLevelBeginnings[profLevel];
      event->accumulated += time_elapsed (&event->startAt, &endAt);
      --profLevel;
   }
}

// ---------------------------------------------------------------------------
/// Registers the end of the block and starts new block (kind of separator with
/// just one call to timer).
// ---------------------------------------------------------------------------
void
profiler_endBegin (int id)
{
#if defined (mc_debugChecks) && mc_debugChecks == 1
   ENSURE (profFile, "profiler is not initialized");
   ENSURE (profLevel >= 0 && profLevel < mc_maxLevel,
           "bad level (%d)", profLevel);
   ensure (profEvent < profEvents + mc_eventCacheN,
           "small cache, %d events already registered", profEvent - profEvents);
#endif

   // Remembers finalization time of event.
   timeTick_t endAt;
   time_get (&endAt);

   // End part.
   if (profAccumEvent) {
      sample_t *event     = profAccumLevelBeginnings[profAccumLevel];
      event->accumulated += time_elapsed (&event->startAt, &endAt);
      if (event <= profAccumEvent)
         profAccumEvent = 0;

      // New blocks should be closed as usual.
      if (profLevel == profAccumLevel
      &&  profAccumLevelBeginnings[profLevel] == profLevelBeginnings[profLevel]) {
         --profLevel;
      }

      --profAccumLevel;
   } else {
      sample_t *event     = profLevelBeginnings[profLevel];
      event->accumulated += time_elapsed (&event->startAt, &endAt);
      --profLevel;
   }

   // Begin part.
   if (profAccumEvent) {
      sample_t *event = profAccumEvent;
      while (event < profEvent && event->id != id) {
         ++event;
      }

      if (event < profEvent) {
         // Remembers start time and event type.
         event->startAt = endAt;
         profAccumLevel = event->level;
         profAccumLevelBeginnings[profAccumLevel] = event;
         return;
      }
   }

   // Remembers start time and event type.
   profEvent->startAt = endAt;
   profEvent->id = id;
   profEvent->level = profAccumLevel = ++profLevel;
   // Remembers descriptor of the start of event.
   profAccumLevelBeginnings[profLevel] = profLevelBeginnings[profLevel]
                                       = profEvent;
   ++profEvent;
}

// ---------------------------------------------------------------------------
/// Signals to profiler that new simulation loop is started.
// ---------------------------------------------------------------------------
void
profiler_startLoop (void)
{
#if defined (mc_debugChecks) && mc_debugChecks == 1
   ENSURE (profFile, "profiler is not initialized");
#endif

   // Resets cache.
   profEvent = profEvents;
   profLevel = -1;
   profAccumEvent = 0;

   // Opens main profiler interval.
   profiler_begin (mc_prof_main);
}

// ---------------------------------------------------------------------------
/// Signals to profiler that simulation loop is finished.
// ---------------------------------------------------------------------------
void
profiler_finishLoop (void)
{
   // Closes main profiler interval.
   profiler_end ();

   // Empty strings for output formatting.
   char   empty1[60],
          empty2[60];
   double invBlockTime[mc_maxLevel] = {1};

   // Prepares empty string to print indents.
   memset (empty1, ' ', sizeof (char)*60);
   memset (empty2, ' ', sizeof (char)*60);
   fprintf (profFile, "--------------- %d intervals ---------------\n",
                        (int) (profEvent - profEvents));
   sample_t *event = profEvents;
   for ( ; event < profEvent ; ++event) {
      const int level = event->level;
      const nameBind_t * const name = profBinds + event->id;

      empty1[level*2] = 0;							// Sets length of the offset.
      empty2[profLongestName - name->nameLenght + 8 - 2*level] = 0;		// Sets span between to reach the next column.

      const double dt     = event->accumulated;					// Alias for time interval.
      invBlockTime[level] = 100.0/(dt + 1e-10);					// 100 to turn everything into percents.
      const int levelDn   = (level > 0) ? level - 1 : 0;
      fprintf (profFile, "%s%s: %s%8.4f%% %.4e sec %8.4f%%\n",
               empty1, name->name, empty2, dt*invBlockTime[levelDn],
               dt, dt*invBlockTime[0]);

      empty1[level*2] = ' ';												// Restores prepared line for offset.
      empty2[profLongestName - name->nameLenght + 8 - 2*level] = ' ';							// Restores prepared line for column shift.
      event->accumulated = 0;												// Prepares for next session.
   }

#if defined (mc_debugChecks) && mc_debugChecks == 1
   ENSURE (profLevel == -1, "%d blocks are not closed properly", profLevel + 1);
#endif

   profEvent      =  profEvents;
   profLevel      = -1;
   profAccumEvent =  0;
}
