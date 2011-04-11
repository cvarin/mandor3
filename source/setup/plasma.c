/** \file plasma.c
  * All particles are created by this interfaces. This way plasma objects can
  * be accessed independently (they are enumerated). I think that plasma objects
  * would be managed in Python in near future, so this module will be removed.
  */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <limits.h>

#include "type_marker.h"
#include "core/plasma_rho.h"		// For charge density evaluation.

#include "log.h"
#include "misc_cfgReader.h"

// Page size; storage size grows by big (~11Mb) steps to avoid many reallocs.
#define MC_PLASMA_PAGE 200000

typedef struct plasmaObject_s {
   marker_t *p;				///< Markers.
   long int  N;				///< Number of particles in the object.
   long int  i;				///< Current end of chunk.
   struct plasmaObject_s *next;		///< Pointer to the next object.
} plasmaObject_t;

// One-way linked list of all plasma chunks.
static plasmaObject_t *plasma = NULL;

// ---------------------------------------------------------------------------
/// Creates new plasma group to fill.
// ---------------------------------------------------------------------------
void
plasma_newObject (void)
{
   // Returns free memory back to the system.
   if (plasma && plasma->N - plasma->i > 0.1*MC_PLASMA_PAGE) {
      plasma->N = plasma->i;
      plasma->p = (marker_t *) realloc (plasma->p, plasma->i*sizeof (marker_t));
      assert (plasma->p && "cannot realloc with smaller RAM-block size!");
   }

   plasmaObject_t *obj = (plasmaObject_t*) calloc (1, sizeof (plasmaObject_t));
   assert (obj);
   obj->next = plasma;
   plasma = obj;
}

// ---------------------------------------------------------------------------
/// Gets pointer to the next unitialized particle.
// ---------------------------------------------------------------------------
marker_t*
plasma_marker (void)
{
   ENSURE (plasma, "plasma storage is not created");
   if (plasma->i >= plasma->N) {
      plasma->N += MC_PLASMA_PAGE;
      plasma->p  = (marker_t *) realloc (plasma->p,
					 plasma->N*sizeof (marker_t));
      ENSURE (plasma->p, "cannot grow storage up to %d particles", plasma->N);
   }
   return plasma->p + (plasma->i++);
}

// ---------------------------------------------------------------------------
/// Returns object number 'count' counting backward from last to the first.
// ---------------------------------------------------------------------------
marker_t*
plasma_getObject (int count, long int *N)
{
   plasmaObject_t *obj = plasma;
   for ( ; obj && count > 0 ; --count)
      obj = obj->next;

   if (!obj) {
      // \TODO: make N = -1, it will screw the result to signal the error.
      *N = 0;
      return NULL;
   }

   *N = obj->i;
   return obj->p;
}

// ---------------------------------------------------------------------------
/// Saves all markers into a file with given 'name'.
// ---------------------------------------------------------------------------
void
plasma_save (const char *name)
{
   long int N  = 0;
   FILE    *fp = cfg_open (name, "wb", __func__);

   // Dummy writing to reserve slot for total number of particles.
   fwrite (&N, sizeof (long int), 1, fp);

#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
   SAY_WARNING ("Tracer info is erased and regenerated.");
#endif

   // Saves all plasma chunks.
   plasmaObject_t *obj = plasma;
   while (obj) {
#if defined(ACTIVATE_TRACER) && ACTIVATE_TRACER
      for (long int i = 0 ; i < obj->i ; i++)
	 obj->p[i].id = N + i + 1;
#endif
      if (obj->i)
	 fwrite (obj->p, sizeof (marker_t), obj->i, fp);
      N  += obj->i;
      ENSURE (N < INT_MAX, "too many particles, integer overflow in tracer");
      obj = obj->next;
   }

   // Writes real number of markers in the beginning of file.
   fseek  (fp, 0, SEEK_SET);
   fwrite (&N, sizeof (long int), 1, fp);
   fclose (fp);
}

// ---------------------------------------------------------------------------
/// Saves charge density as tecplot wants.
// ---------------------------------------------------------------------------
void
plasma_rho (meshDouble_p rho)
{
   mf_mesh_clean (rho);		// Clears mesh accumulator.
   plasmaRho_setupConnections ();	// Resets plasma system.

   plasmaObject_t *obj = plasma;
   while (obj) {
      if (obj->i)
	 plasmaRho_add (rho, obj->p, obj->i);
      obj = obj->next;
   }

   // Parallel exchange and boundary conditions inside.
   plasmaRho_startExchange  (rho);
   plasmaRho_finishExchange (rho);
}
