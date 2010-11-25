/** \file log.h
  * \Todo:
  * Ideas:
  *    + Merge consequtive debug messages to have the save origin.
  *    + Add 'SAY_DDEBUG', 'SAY_DDDEBUG' and optionally enable with '-v' '-vv'.
  *    + Strip them out of file/func/line info to gain the speed of 'printf'.
  */

#ifndef MC_LOG_HEADER
#define MC_LOG_HEADER

#include <stdlib.h>

void log_open    (int be_silent, int continue_log, const char *logname);
void log_close   (void);

void say         (const char *format, ...);
void say_doing   (const char *format, ...);
void say_warning (const char *func, const char *file, int line,
                  char *format, ...);
void say_debug   (const char *func, const char *file, int line,
                  char *format, ...);
void say_OMG     (const char *func, const char *file, int line,
                  char *format, ...);
void die         (const char *func, const char *file, int line,
                  const char *fail, char *format, ...);

/// "sprintf" wrapper with result persisting for at least 31 call.
const char* _ (const char *format, ...);

// Wrappers around 'say_...' to automatically add exact location of caller.
#define SAY_WARNING(...) say_warning(__func__, __FILE__, __LINE__, __VA_ARGS__)
#define SAY_DEBUG(...)   say_debug  (__func__, __FILE__, __LINE__, __VA_ARGS__)
#define SAY_OMG(...)     say_OMG    (__func__, __FILE__, __LINE__, __VA_ARGS__)

// ---------------------------------------------------------------------------
/// Macro-function to assure that critical condition is not violated.
///
/// Same like 'assert', but provides user-friendly report and thus retires
/// all feed-back message generation.
///
/// \TODO If MPI/OpenMP defines any macro to indicate it's precence then I can
/// use '#if' to replace 'exit()' by 'MPI_Abort()' to exit gracefully.
// ---------------------------------------------------------------------------
#define ENSURE(test, ...)						\
{									\
   if (! (test)) { 							\
      die  (__func__, __FILE__, __LINE__, #test, __VA_ARGS__);		\
      exit (EXIT_FAILURE);						\
   }									\
}

// ---------------------------------------------------------------------------
/// Explicit termination to pin-point unfinished backend, etc.
// ---------------------------------------------------------------------------
#define DIE(...)     							\
{									\
   die  (__func__, __FILE__, __LINE__, NULL, __VA_ARGS__);		\
   exit (EXIT_FAILURE);							\
}

#endif // MC_LOG_HEADER
