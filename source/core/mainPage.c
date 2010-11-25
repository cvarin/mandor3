/** \file mainPage.c
  * \brief \b Doxygen formatted main page of the project.
  */

/** \mainpage Mandor: numerical core.
  * Main module reads distribution function and fields from disk and performs
  * specified number of time steps. Results are written to hard drive for
  * subsequent analysis.
  *
  * Diagnostics are build-in:
  * - W(t) - full energies of the fields and particles;
  * - snapshot of meshes (currents, EM-fields) for subsequent visualisation;
  * - checkpoints for visualisation and analisis of the total system state.
  *
  * Structure of sources (major groups and libraries):
  *
  * - Geometry manager: collection of routines for fast partitioning,
  *   overlapping, analizing, exchanging, printing, adding to lists of
  *   geometric entities which may be point, line, rectangular of volume with
  *   faces and edges parallel to coordinate axises (type_reg.c).
  *
  * - Mesh manager: collection of routines for allocating, resizing, loading,
  *   saving, accessing, etc of the unrolled 1D/2D/3D uniform mesh. Index
  *   ranges are arbitrary, for 1D/2D simulations efficient use of macroses
  *   removes unnecessary stuff without sacrificing readability or speed
  *   (type_mesh.c).
  *
  * - socket layer - collection of routines to build and maintain non-blocking
  *   out-of-order (fist in - first served) exchange (misc_socket.c).
  *
  * - IO library: IO_sys.c, IO_tecplot.c.
  *
  * - Hierarchical domain decomposition: misc_partition.c.
  *
  * - Maxwell solver (em.c, em_caps.c, em_TFSF.c, etc)
  *
  * - Boris particle pusher (plasma.c, plasma_walls.c, plasma_currentWalls.c, etc).
  *
  * Parameters of the node are const global variables accessible from anywhere
  * (see misc_parameters.c).
  *
  * \warning type_mesh.c, misc_partition.c and few old libraries will be
  * refined using new geometry engine.
  *
  * Structure of the program is well hinted by the body of the main() function.
  */
