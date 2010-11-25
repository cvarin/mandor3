/** \file timer.h
  * \brief Run-time timer with high-accuracy. If you need more multi-platform variants
  * take a look on the fftw3 library (check ./kernel/cycle.h file).
  *
  * \todo It is possible to measure real time and per-process time using, for example,
  * \b CLOCK_REALTIME and \b CLOCK_PROCESS_CPUTIME_ID. It may points to a regions with
  * bottlenecks in the system part of the execution process.
  */

#ifndef MC_TIMER_HEADER
#define MC_TIMER_HEADER							///< \internal Guard.

#include <math.h>
#include <sys/time.h>

/**
  * High accuracy timer structure.
  */
typedef struct timeval timeTick_t;

/**
  * Writes current time into the object provided.
  */
static inline void
time_get (timeTick_t *time)
{
  gettimeofday (time, NULL);
}

/**
  * Converts time to seconds.
  */
static inline double
time_seconds (timeTick_t *time)
{
  return time->tv_sec + time->tv_usec*1e-6;
}

/**
  * Returns time elapsed between 2 events (in seconds).
  */
static inline double
time_elapsed (timeTick_t *timeStart, timeTick_t *timeEnd)
{
  return timeEnd->tv_sec - timeStart->tv_sec + (timeEnd->tv_usec - timeStart->tv_usec)*1e-6;
}

/**
  * Converts time to seconds.
  */
static inline double
time_resolution (double *sigma)
{
  timeTick_t res[101];													// Reports resolution of the timer used.

  for (int i = 0 ; i < 101 ; ++i)
    time_get (res + i);

  double avr = 0, avr2 = 0;
  for (int i = 1 ; i < 101 ; ++i)
  {
    double t = time_elapsed (res + i - 1, res + i);
    avr += t;
    avr2 += t*t;
  }

  avr *= 0.01;
  avr2 *= 0.01;
  *sigma = sqrt (avr2 - avr*avr);

  return avr;
}

#endif
