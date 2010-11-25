/** \file plasma_parallel.h
  * Plasma particles parallel exchange framework.
  */

#ifndef MC_PLASMA_PARALLEL_HEADER
#define MC_PLASMA_PARALLEL_HEADER

void comm_plasma_configure  (void);
void comm_plasma_start      (void);
void comm_plasma_send       (void);
void comm_plasma_test_inbox (void);
void comm_plasma_complete   (int nanosecs);
void comm_plasma_dump       (void);

#endif
