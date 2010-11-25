/** \file mainPage.c
  * \brief Just main page formatted for \b Doxygen.
  */

/** \mainpage Setup module of the "Mandor".
  * This is the setup utility of the package. Initial state is created using configuration file containing basic elements - \b tags.
  * Tags are used to add physically solid element of the simulation - wave, plasma component, perturbation and so on. In the code
  * each \b tag is processed independently by call to the interface function which also performs all necessary checks of the input
  * parameters and diagnostic output.
  *
  * \b setup.out may be executed in <em>memory estimate</em> mode. In this mode initial state is not created but memory usage
  * is estimated and reported end few checks of the input parameters are performed if possible to check correctness of the tags parameters.
  * Main goal of this mode is to tune parameters to satisfy hardware restrictions on the site. For example, sequence of runs may be:
  * - choose set of tags, run setup in <em>memory estimate</em> mode;
  * - fix errors in the parameters (if any reported);
  * - see what is the memory consumptions, time/spatial steps, parameters of the objects created and so on;
  * - refine parameters (using full picture of the simulation run sketched), check again;
  * - create initial conditions using updated config file.
  *
  * The most important design notes:
  * - Output of the tag must be informative - this output is the main source of the information about physical system created, better
  *   even than config file itself. Remember to use human attention wisely - don't output too much of an irrelevant technical information.
  * - All problem related stuff seats in this module only. No hacks in the \b core.out (modelling binary). Keeping kernel unmodified
  *   through all history of testing guarantees that code used is the code tested.
  * - Remember to account for the initial time shift of the magnetic field and velocity.
  *
  * <b>P.S.</b> The latest list of all supported tags is simply a list of all predefined constants in tag.h.
  */
