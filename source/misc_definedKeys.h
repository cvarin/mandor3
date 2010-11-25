#ifndef MC_DEBUGSWITCHES_HEADER
#define MC_DEBUGSWITCHES_HEADER

#define mc_debug_currentsKernelRangeCheck	0			///< Activates (if != 0) checks of the range in the current density kernel.

/* ============================== EM part ============================== */
/* Switch to activate checking of the boundaries at "total/scatter" interface */
#define player_regionBoundariesTest
#undef  player_regionBoundariesTest

#endif
