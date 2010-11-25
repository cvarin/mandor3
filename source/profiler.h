/** \file profiler.h
  * Run-time profiler to measure execution time of big logical blocks.
  *
  * Using pair of function calls to mark begin and end of the interval. All
  * safety checks can be removed to reduce perturbations. Main goal is to
  * briefly review the big picture and to guess the performance issues.
  *
  * Key points:
  * - begin/end approach;
  * - use it in non-performance critical sections only (profiler may be
  *   deactivated at compile time to estimate overheads);
  * - all sample intervals are named with enumerated constants at compile-time;
  * - accumulation helps to collect total statistic;
  * - simple stack is used for events, additional stack exists to accumulate.
  *
  * \todo It is possible to measure real time and per-process time using,
  *       for example, \b CLOCK_REALTIME and \b CLOCK_PROCESS_CPUTIME_ID.
  *       It may points to a regions with bottlenecks in the system rather
  *       than in the code.
  *
  * REFACTORING
  * ===========
  * 1) Use macros to call 'prof_enter (const char *__func__)' to substitute
  *    __func__ as argument.
  * 2) Save raw timestamp '{ .name = __func__, .time = time () }'. Do not
  *    analyse or strcopy func-name, save only the pointer and use it in final
  *    call where we do dump data.
  * 3) Use profiler by adding 'PROFILE_ME' at the beginning of the function.
  * 4) Use gtk analyzer.
  */

#ifndef MC_PROFILER_HEADER
#define MC_PROFILER_HEADER

/// Use this branching to disable profiling (to estimate the intrusion).
#if defined(mc_blockDummyMacroDeclarations) || 1
  void profiler_init       (const char *name);
  void profiler_syncClocks (void);
  void profiler_begin      (int id);
  void profiler_end        (void);
  void profiler_endBegin   (int id);
  void profiler_accumulate (void);
  void profiler_startLoop  (void);
  void profiler_finishLoop (void);
#else
  /// Dummy prototypes to efficiently strip the profiler calls.
  #define profiler_init(name) 		{}
  #define profiler_syncClocks()		{}
  #define profiler_begin(id)		{}
  #define profiler_end()		{}
  #define profiler_endBegin(id)		{}
  #define profiler_accumulate()		{}
  #define profiler_startLoop()		{}
  #define profiler_finishLoop()		{}
#endif

/// Size limit for the name of the block.
#define mc_prof_nameSize 32

#define mc_prof_main			0	///< ID of the root interval.
#define mc_prof_emh1			1	///< em_HHalfStep().
#define mc_prof_eme			2	///< em_EStep().
#define mc_prof_emh2			3	///< em_HStep().
#define mc_prof_throwStopFlag		4	///< Sending stop flag request.
#define mc_prof_EMenergy		5	///< Evaluating total energy of EM field.
#define mc_prof_probes			6	///< Probes diagnostic.
#define mc_prof_divJ_prep		7	///< Preparation of the auxilary data for charge conservation test.
#define mc_prof_plasma_move		8	///< Pushing of particles.
#define mc_prof_plasma_Txyz		9	///< Getting precalculated energy distribution over degrees of freedom.
#define mc_prof_tecRho			10	///< Calculating rho for tecplot output.
#define mc_prof_tecEH			11	///< Saving tecplot output (E and H).
#define mc_prof_spectr			12	///< Saving spectr dump.
#define mc_prof_wDensity		13	///< Saving energy density.
#define mc_prof_jbcFinish		14	///< Receiving currents.
#define mc_prof_EFinish			15	///< Receiving fields.
#define mc_prof_divJ			16	///< Dumping charge conservation test results.
#define mc_prof_tecRhoJ			17	///< Writing charge and current density for tecplot/openDX.
#define mc_prof_sysSave			18	///< System check-point.
#define mc_prof_wDensityFlush		19	///< Flushing total energy distribution data.
#define mc_prof_defrag			20	///< Defragmenting particles storages.
#define mc_prof_catchStop		21	///< Catching stop flag.

#define mc_prof_HHalf_TFSF_preH		22	///< emh1: TFSF.
#define mc_prof_HHalf_cacheEH		23	///< emh1: cap_cache.
#define mc_prof_HHalf			24	///< emh1 itself.
#define mc_prof_HHalf_flushH		25	///< emh1: bc.

#define mc_prof_plasma_cleanJ		26	///< Plasma - initial cleaning of the J array.
#define mc_prof_plasma_mainLoop		27	///< Main loop - time step.
#define mc_prof_plasma_pbcPop		28	///< Main loop - sending particles.
#define mc_prof_plasma_pbcRecv		29	///< Main loop - receiving particles.
#define mc_prof_plasma_jbc		30	///< Plasma - current boundary conditions.

#define mc_prof_plasma_jbc_mesh		31	///< Exchange by j-mesh.
#define mc_prof_plasma_jbc_particles 	32	///< Receiving of the particles send.
#define mc_prof_plasma_jbc_merge 	33	///< Merging chapter with received particles with the body of the plasma.
#define mc_prof_plasma_jbc_sweep 	34	///< VSP sweep pass.
#define mc_prof_plasma_jbc_localBC 	35	///< Local BC applications before sweep pass.

#define mc_prof_plasma_ebc_cap		36	///< Exchange by EM-mesh data.
#define mc_prof_plasma_ebc_recv 	37	///< Receiving of the particles send.
#define mc_prof_plasma_ebc_jPass 	38	///< Merging chapter with received particles with the body of the plasma.
#define mc_prof_plasma_ebc_misc 	39	///< VSP sweep pass.

#define mc_prof_N			40	///< Total number of the hand-coded events.

#endif
