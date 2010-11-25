#ifndef MC_MISC_TAGS_HEADER
#define MC_MISC_TAGS_HEADER

/// MPI tags are gathered here to ensure uniqueness.
enum {TAG_PLASMA_RHO = 10, TAG_PLASMA_J, 		// Plasma timestep exchange.
      TAG_SYNC_REGLISTS, 				// Reg-list syncronization.
      TAG_STOP_REQUEST, 				// Core.out timestep thing.
      TAG_EM_SPLITTER_E, TAG_EM_SPLITTER_H, 		// EM-solver timestep exhanges.
      TAG_GHOST_SYNC_E, TAG_GHOST_SYNC_H,		// Ghost cell explicit syncronization.
      TAG_IO_BUFFER_N, TAG_IO_BUFFER_DATA,		// Plasma loading in core.out.
      TAG_MARKERS_N, TAG_MARKERS_DATA, 			// Plasma markers time-step exchange.
};

#endif
