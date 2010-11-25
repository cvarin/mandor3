/** \file em_sources.h
  * Global interfaces to cache and use soft sources of EM waves.
  *
  * Mirrors:
  *   + 'setup.out' calls 'add_new_mirror' to add new mirror to the system.
  *   + Stratton-Chu integrals are precomputed during setup and saved on HDD.
  *   + The solver ('core.out') loads this data and uses it.
  */

#ifndef EM_SOURCES_HEADER
#define EM_SOURCES_HEADER

#include "type_reg.h"
#include "type_mesh.h"

#include "em_mirror.h"

#ifdef hz       // For IBM AIX compiler.
#undef hz
#endif

// Setup interfaces.
void add_new_mirror (mirror_t mirror);

// Core interfaces.
void init_sources  (void);
void add_E_sources (meshVec_p E, double time, reg_t *to_update);
void add_H_sources (meshVec_p H, double time, reg_t *to_update);

#endif
