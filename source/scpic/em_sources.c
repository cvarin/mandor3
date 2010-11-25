/** \file em_sources.c
  * Global interfaces to manage sources of EM waves.
  *
  * For now it manages only laser beams, focused using parabolic mirrors.
  * Later I can add plane waves or something else, so I keep this layer.
  *
  * CREDITS
  * =======
  * Mirror code was taken from Kostya Popov's fork of original Mandor.
  * Code was seriously improved (see 'em_mirror.h' for details), but thanks to
  * him I had motivation, working implementation of the idea in a production
  * code (which is always a very solid foundation for the work).
  *
  * Kostya's code is SCPIC: http://www.phys.ualberta.ca/~kpopov/scpic.html.
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "em_sources.h"

#include "log.h"
#include "misc_units.h"
#include "misc_cfgReader.h"
#include "misc_parameters.h"

/// Pack of laser beams focused by parabolic mirrors.
static int       mirrors_N = 0;
static mirror_t *mirrors   = NULL;

// ----------------------------------------------------------------------------
/// Initializes all local sources (loads cached Stratton-Chu integrals, etc).
// ----------------------------------------------------------------------------
void
init_sources (void)
{
   ENSURE (!mirrors_N && !mirrors, "double initialization");

   // Loads all laser beams focused by parabolic mirrors.
   while (1) {
      mirror_t src;
      if (! mirror_source_is_loaded (mirrors_N, &src))
         break;

      mirrors = (mirror_t*) realloc (mirrors, (++mirrors_N)*sizeof (mirror_t));
      mirrors[mirrors_N-1] = src;

      print_mirror_parameters (mirrors[mirrors_N-1]);
   }
}

// ----------------------------------------------------------------------------
/// Adds new laser beam focused by parabolic mirror.
/// Used by 'setup.out' to compute and save data on HDD (see 'load_mirrors').
// ----------------------------------------------------------------------------
void
add_new_mirror (mirror_t mirror)
{
   mirrors_N++;
   mirrors = realloc (mirrors, mirrors_N*sizeof (mirror_t));
   assert (mirrors && "cannot allocate tiny array");
   mirrors[mirrors_N-1] = mirror;

   save_mirror_source (mirrors_N - 1, mirror);
}

// ----------------------------------------------------------------------------
/// Adds contribution from all soft sources after completing Yee E timestep.
// ----------------------------------------------------------------------------
void
add_E_sources (meshVec_p E, double time, reg_t *to_update)
{
   for (int m = 0; m < mirrors_N ; m++) {
      add_mirror_E (mirrors + m, E, time, to_update);
   }
}

// ----------------------------------------------------------------------------
/// Adds contribution from all soft sources after completing Yee H timestep.
// ----------------------------------------------------------------------------
void
add_H_sources (meshVec_p H, double time, reg_t *to_update)
{
   for (int m = 0; m < mirrors_N ; m++) {
      add_mirror_H (mirrors + m, H, time, to_update);
   }
}
