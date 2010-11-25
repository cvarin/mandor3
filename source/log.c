/** \file log.c
  * Message subsystem which lets you:
  *   + print messages
  *   + copy output into log file
  *   + do 'debug print' (data goes to the log file w/o garbaging the screen)
  *   + explicitly silenced logging (to use under MPI on slave nodes)
  *   + message level helps to separate debug, warning, and info messages
  *
  * \warning: lines are automatically terminated with '\n'.
  *
  * Note: 'say_doing' updates it's own last message with reduced frequency to
  *       protect terminal IO bandwidth (unprinted data are cached in 'cache').
  */

#include <sys/time.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>
#include <assert.h>

#include "log.h"

/// ANSI terminals support this escape-sequence commands.
#define ANSII_CURSOR_UP		"\033[A"
#define ANSII_ERASE_LINE	"\033[2K"

#define BUF_SIZE  	200
#define KEEP_MSG_POWER  8	// '_' keeps (2<<KEEP_MSG_POWER) messages.

enum { DEBUG = 0, INFO, STATUS, WARNING, ERROR };

static FILE *logfile         = NULL;
static int   be_quiet        = 0;	// <== To suppress screen output on
static int   warnings_posted = 0;   	//     slave nodes under MPI.

static int   lines_to_clean = 0;   	// Flag to rule the clean-up.

static char  msg  [BUF_SIZE + 2] = "";	///< Common 'vsnprintf' buffer.
static char  cache[BUF_SIZE + 2] = "";	///< Unflushed status string.

// Prints data to msg/cache/_ buffer, tests against truncation, adds tail.
#define VSNPRINTF_TO(BUF, TAIL_CHAR) {					\
   va_list  argptr;							\
   va_start (argptr, format);						\
   int printed = vsnprintf (BUF, BUF_SIZE, format, argptr);		\
   va_end   (argptr);							\
									\
   if (printed >= BUF_SIZE - 1) {					\
      SAY_WARNING ("message was truncated to %d chars", BUF_SIZE);	\
      BUF[BUF_SIZE]   = 0;				        	\
      BUF[BUF_SIZE-1] = 0;				        	\
   } else {								\
      BUF[printed]   = TAIL_CHAR;				       	\
      BUF[printed+1] = 0;         					\
   }									\
}

// ---------------------------------------------------------------------------
/// Saves message to the logfile with automatic prefixing of each line.
// ---------------------------------------------------------------------------
static void
add_to_log (int level, const char *buf)
{
   if (!logfile)
      return;

   assert (buf[BUF_SIZE] == 0);

   static char* prefix[5] = {"DEB: ", "INF: ", "   : ", "WRN: ", "ERR: "};

   const char *pos         = buf;
   int         need_prefix = 1;
   while (*pos) {
      if (need_prefix)
         fputs (prefix[level], logfile);
      fputc (*pos, logfile);
      need_prefix = (*pos == '\n');
      pos += 1;
   }
}

// ---------------------------------------------------------------------------
/// Flushes cache and prints message.
// ---------------------------------------------------------------------------
void
say (const char *format, ...)
{
   // Flushes cached message.
   if (!be_quiet && cache[0]) {
      fputs (cache, stdout);
      cache[0] = 0;
   }

   lines_to_clean = 0;

   // Prints the message.
   VSNPRINTF_TO (msg, '\n');
   add_to_log (INFO, msg);

   if (!be_quiet)
      fputs (msg,  stdout);
}

// ---------------------------------------------------------------------------
/// Updates progress bar (overwrites the old message using \ESC commands).
/// Updating frequency is limited to avoid printing bottlenecks.
// ---------------------------------------------------------------------------
void
say_doing (const char *format, ...)
{
   // Caches the message.
   VSNPRINTF_TO (cache, '\n');
   add_to_log (STATUS, cache);

   if (be_quiet)
      return;

   static struct timeval old = {0, 0}, now;
   gettimeofday (&now, NULL);
   if (now.tv_sec - old.tv_sec + (now.tv_usec - old.tv_usec)*1e-6 < 0.25)
      return;

   while (lines_to_clean--)
      fputs (ANSII_CURSOR_UP ANSII_ERASE_LINE, stdout);

   lines_to_clean = 0;
   char *pos = cache;
   while ((pos = strchr(pos, '\n'))) {
      lines_to_clean += 1;
      pos += 1;
   }

   fputs (cache, stdout);
   cache[0] = 0;
   old = now;
}

// ---------------------------------------------------------------------------
/// Issues warning message and counts it (to draw attention at the end of log).
// ---------------------------------------------------------------------------
void
say_warning (const char *func, const char *file, int line, char *format, ...)
{
   VSNPRINTF_TO (msg, '\0');
   add_to_log (WARNING, _("%s\n[ %s @ %s:%d ]\n", msg, func, file, line));
   warnings_posted++;

   if (!be_quiet)
      printf ("WRN: %s\n", msg);
}

// ---------------------------------------------------------------------------
/// Sends debug message directly into the file.
// ---------------------------------------------------------------------------
void
say_debug (const char *func, const char *file, int line, char *format, ...)
{
   // Generates message string.
   VSNPRINTF_TO (msg, '\n');
   add_to_log (DEBUG, msg);
}

// ---------------------------------------------------------------------------
/// Issues error message (the same warning counter is reused).
// ---------------------------------------------------------------------------
void
say_OMG (const char *func, const char *file, int line, char *format, ...)
{
   // Generates message string.
   VSNPRINTF_TO (msg, '\0');
   add_to_log (ERROR, _("%s\n[ %s @ %s:%d ]\n", msg, func, file, line));
   warnings_posted += 1000000;

   fflush (logfile);
   fsync (fileno (logfile));

   if (!be_quiet)
      printf ("ERR: %s\n", msg);
}

// ---------------------------------------------------------------------------
/// Issues final note and dies.
// ---------------------------------------------------------------------------
void
die (const char *func, const char *file, int line,
     const char *fail, char *format, ...)
{
   const char *info = fail ? _("condition '%s' is failed", fail)
                           : "explicit exit on critical error";

   VSNPRINTF_TO (msg, '\0');
   add_to_log (ERROR, _(">>>>> Epic fail: %s.\n", msg));
   add_to_log (ERROR, _("  Info: %s.\n"
                        "  File: %s\n"
                        "  Func: %s\n"
                        "  Line: %d\n", info, file, func, line));
   add_to_log (ERROR, _("<<<<< Epic fail: %s.\n", msg));
   warnings_posted += 1000000;

   fflush (logfile);
   fsync (fileno (logfile));

   printf ("\n"
           "ERROR: %s\n"
           "    Info: %s.\n"
           "    File: %s\n"
           "    Func: %s\n"
           "    Line: %d\n\n", msg, info, file, func, line);
}

// ---------------------------------------------------------------------------
/// Closes all logfiles (intended for 'atexit').
// ---------------------------------------------------------------------------
void
log_close (void)
{
   // 'say' automatically flushes cached buffer, if necessary.
   say ("-------------- END OF LOG --------------");
   say ("  %d error(s), %d warning(s)",
         warnings_posted/1000000, warnings_posted % 1000000);

   if (logfile)
      fclose (logfile);
   logfile = NULL;
}

// ---------------------------------------------------------------------------
/// Sets name of the log-file, submits 'msg_shutdown' and so on.
// ---------------------------------------------------------------------------
void
log_open (int be_silent, int continue_log, const char *logname)
{
   be_quiet = be_silent;

   if (logname) {
      logfile = fopen (logname, continue_log ? "at" : "wt");
      ENSURE (logfile, "cannot open '%s' logfile for writing", logname);

#if defined(MC_FORCE_FILE_FLUSHES) && MC_FORCE_FILE_FLUSHES
      setvbuf(logfile, NULL, _IONBF, 0);
      WARN ("unbuffered log IO may cause performance penalty");
#endif
   }

   static int submitted = 0;
   if (!submitted)
      ENSURE (!atexit (log_close), "cannot submit 'log_close'");
   submitted = 1;
}


#if KEEP_MSG_POWER < 5 || (2<<KEEP_MSG_POWER)*BUF_SIZE > 1<<20
#  error "Too small (or too big) buffer for '_'-wrapped messages."
#endif

// ---------------------------------------------------------------------------
/// Simular to 'sprintf' but prints into internal buffer which is safe to use
/// for at least 2^KEEP_MSG_POWER calls to this function: we can use the
/// returned string in nested calls few times :-).
// ---------------------------------------------------------------------------
const char*
_ (const char *format, ...)
{
   static char str[1<<KEEP_MSG_POWER][BUF_SIZE + 1];
   static int  slot = -1;
   slot = (slot + 1) & ((1<<KEEP_MSG_POWER) - 1);

   VSNPRINTF_TO (str[slot], '\0');
   return str[slot];
}
