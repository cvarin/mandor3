.POSIX:

CCOMPILER=mpicc


# Use this keys to activate debug/fallback options inside of the code
# (recompile the code to apply all changes).

# MC_FILL_MESH_WITH_NAN: fills allocated mesh arrays with qNaNs to catch
# uninitialized mesh cell contributions

# MC_GAUSS_TEST: checks charge conservation by computing new and old
# charge density and testing d(rho)/dt + div(j)

# MC_GAUSS_DEBUG: in addition to saving 'max (d(rho)/dt + div(j))' it
# outputs full array of actual charge density balance. Use it only if
# the previous test is failed: it generates full dump on every timestep <=>
# VERY BIG FILES. This options activates MC_GAUSS_TEST automatically.

# MC_NO_MPI_ALLOC_MEM: forces to use malloc/free instead of MPI_Alloc_mem
# and MPI_Free_mem (MPI on one of the clusters doesn't have it yet!)

# MC_CURRENT_KERNEL_RANGE_TEST: LOW-LEVEL DEBUG TEST of the array boundaries
# violations in plasma timestep, performance penalties!

# MC_FORCE_FILE_FLUSHES: forces to flush log-file after each logged message
# (to catch messages issued just before the code crushes).

# Apply it like 'TEST = -DMC_FILL_MESH_WITH_NAN -DMC_NO_MPI_ALLOC_MEM'.

TESTS = -DMC_FILL_MESH_WITH_NAN
        #-DMC_GAUSS_TEST
        #-DMC_CURRENT_KERNEL_RANGE_TEST
        #-DMC_FORCE_FILE_FLUSHES
COMMON = -Isource -Wuninitialized -Wunused -pedantic -std=gnu99 -DLINUX -pipe $(TESTS)

# Debug options
  OPTIMIZE = -O1 -g -Wall -Wno-unknown-pragmas

# Profile options - gmon
  OPTIMIZE = -O2 -g -pg -Wall -Wno-unknown-pragmas -ffast-math -fno-math-errno -funsafe-math-optimizations -fno-trapping-math

# Profile options - callgrind
  OPTIMIZE = -O2 -g -Wall -Wno-unknown-pragmas -ffast-math -fno-math-errno -funsafe-math-optimizations -fno-trapping-math

# Performance options
  OPTIMIZE = -O4 -g0 -march=native -mtune=native  -ffast-math -fno-math-errno -funsafe-math-optimizations -fno-trapping-math		# Speed options' set

LFLAGS = $(OPTIMIZE) $(COMMON)
CFLAGS = -c $(LFLAGS)
VPATH = source:o

# my_all: all
my_all: setup.out core.out TFSFdefrag.out distiller math.out tecplot.out
	@echo ;
	@echo "All files are compiled using '$(CCOMPILER) $(CFLAGS)'."
	@echo ;

.PHONY: distiller
distiller: o/distiller.out

#
# Information for project dependency scanner.
#

# [megamake] [folders] (src: source/) (compiled: o/)

# [megamake] [scan] @$(CCOMPILER) $(CFLAGS) -MM $<
# [megamake] [compile] @echo "Compiling $@ ..."
# [megamake] [compile] @$(CCOMPILER) $(CFLAGS) $< -o $@
# [megamake] [link] @echo "  Linking $@..."
# [megamake] [link] @$(CCOMPILER) $^ -o $@ $(LFLAGS) -lm ___libs___
# [megamake] [link] @echo "  Executable $@ is ready."
# [megamake] [link] @echo ""
# [megamake] [all] @echo ;
# [megamake] [all] @echo "All files are compiled using '$(CCOMPILER) $(CFLAGS)'."
# [megamake] [all] @echo ;

# [megamake] [header only] (source/misc_PIC.h) (source/misc_MPItags.h) 		\
#   (source/misc_definedKeys.h) (source/type_particle.h) 	\
#   (source/timer.h) (source/mf_vector.h) (source/frame.h) 			\
#   (source/type_particle.h) (source/probe.h) (source/types.h)			\
#   (source/type_marker.h)

# [megamake] [project] (root: source/core/main.c) (exec: core.out)
# [megamake] [project] (root: source/setup/main.c) (exec: setup.out)
# [megamake] [project] (root: source/TFSFdefrag/main.c) (exec: TFSFdefrag.out)
# [megamake] [project] (root: source/distiller/main.c)  (exec: o/distiller.out)
# [megamake] [project] (root: source/mathematica/main.c) (exec: math.out)
# [megamake] [project] (root: source/log_test.c) (exec: test_logger.out)

# [megamake] [project] (root: source/tecplot/main.c) (exec: tecplot.out)
# [megamake-] [project] (root: source/spectr/main.c) (exec: spectr.out) (libs: -lfftw3)
# [megamake-] [project] (root: source/phasedot/main.c) (exec: phasedot.out)
# [megamake-] [project] (root: source/distrib/main.c) (exec: distrib.out)
# [megamake-] [project] (root: source/probe2tec/main.c) (exec: probe2tec.out) (libs: -lfftw3)
# [megamake-] [project] (root: source/test/main_test.c) (exec: test.out)
# [megamake-] [project] (root: source/fft.c) (exec: testFFT.out) (libs: -lfftw3)
# [megamake-] [project] (root: source/fft/main.c) (exec: fft.out) (libs: -lfftw3)
# [megamake-] [project] (root: source/markerExplorer/main.c) (exec: markerExplorer.out)

# Clean up. Wipes everything.
.PHONY: help DEPEND cleanTFSF cleanData cleanAll doc test

help:
	@echo 'That is "Mandor" makefile rules for you:'
	@echo '  o make depend - refreshes dependency information.'
	@echo '  o make all - compiles all executable.'
	@echo '  o make clean - removes compiled objects files and'
	@echo '                 executables (data files are safe!).'
	@echo '  o make cleanTFSF - removes TFSF written data.'
	@echo '  o make cleanData - removes written data files, error messages, etc.'
	@echo '  o make cleanAll - removes EVERYTHING but sources.'
	@echo ;

depend:
	@echo "Making dependencies ..."
	@megamake.pl Makefile --verbose --make:gmake
	@echo "Done."

strip: all
	@echo "Making stripped executables:"
	@for f in *.out ; do echo "    stripping $$f ..." ; strip -s $$f ; done
	@echo "Done."

test: test_logger.out
	@echo "Making all tests..."

cleanTFSF:
	@echo "Cleaning 'Total field / scattering field' interface data..."
	@find ./.EM_sources -type f -name "*" | grep -v .svn | grep -v CVS | while read f ; do rm -f $$f ; done

# File list might be too long for command line, so one erases files one by one.
cleanData:
	@echo "Cleaning ./tmp folder ..."
	@find ./tmp -type f -name "*" | grep -v .svn | grep -v CVS | while read f ; do rm -f $$f ; done
	@echo "Cleaning ./output folder ..."
	@find ./output -type f -name "*" | grep -v .svn | grep -v CVS | while read f ; do rm -f $$f ; done
	@echo "Cleaning ./binData folder ..."
	@find ./binData -type f -name "*" | grep -v .svn | grep -v CVS | while read f ; do rm -f $$f ; done
	@echo "Cleaning created files in the project folder ..."
	@rm -f core.[0-9]* .gaussSpot.cfg errors.*

# File list might be too long for command line, so erase files one by one.
cleanAll: clean cleanTFSF cleanData
	@echo "Cleaning ./doc folder ..."
	@rm -fr ./doc/*
	@rm -f ./o/dependSignalFile
	@echo "Everything but distribution is wiped out."

# Builds documentation in doc folder.
doc:
	doxygen .Doxyfile.setup
	doxygen .Doxyfile.core
	doxygen .Doxyfile.probes
	doxygen .Doxyfile.distrib
	doxygen .Doxyfile.spectr
	@echo ;
	@echo "Documentation is ready. Please examine your doc folder."
	@echo ;
	@cd doc/core/html ; firefox index.html

# [megamake] [end of header]



########################################################################
#
# This part is AUTOMATICALLY generated by '/home/henaro/bin/Mandor2/megaMake/megamake.pl' perl script.
# IT WILL BE REMOVED WITHOUT ANY WARNING AFTER THE NEXT CALL TO SCRIPT.
# ADD YOUR VARIABLES/COMMANDS/STUFF IN THE HEADER PART IF YOU NEED IT.
# Date: Fri Oct 15 23:10:46 KRAST 2010
# Project(s) folder: /home/henaro/smyth/Mandor2
#
########################################################################



o/setup/tag_jBalancer.o: source/setup/tag_jBalancer.c source/type_marker.h source/log.h source/misc_cfgReader.h source/misc_parameters.h source/frame.h source/setup/main.h source/setup/plasma.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_foil.o: source/setup/tag_foil.c source/type_marker.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/misc_parameters.h source/setup/main.h source/setup/plasma.h source/type_mesh.h source/misc_partition.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_meanVNoise.o: source/setup/tag_meanVNoise.c source/log.h source/misc_PIC.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/misc_cfgReader.h source/misc_parameters.h source/misc_markerPlacer.h source/setup/main.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_twoStream.o: source/setup/tag_twoStream.c source/type_marker.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_parameters.h source/frame.h source/setup/main.h source/setup/plasma.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag.o: source/setup/tag.c source/setup/tag.h source/setup/tag_units.h source/setup/tag_mesh.h source/setup/tag_boundary.h source/setup/tag_TFSF.h source/setup/tag_TFSFFaces.h source/setup/tag_laser.h source/setup/tag_plasma.h source/setup/tag_DF_uniform.h source/setup/tag_point.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/setup/tag_EMWave.h source/setup/tag_ringDF.h source/setup/tag_vShift.h source/setup/tag_photoelectrons.h source/setup/tag_jBalancer.h source/setup/tag_twoStream.h source/setup/tag_gaussSpot.h source/setup/tag_maxwell.h source/setup/tag_planePulse.h source/setup/tag_plasmaWave.h source/setup/tag_EMResonator.h source/setup/tag_scissors.h source/setup/tag_meanVNoise.h source/setup/tag_seedWeibel.h source/setup/tag_seed_PITS.h source/setup/tag_cluster.h source/setup/tag_foil.h source/setup/tag_scales.h source/setup/tag_trianglePrizm.h source/log.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_vShift.o: source/setup/tag_vShift.c source/log.h source/misc_cfgReader.h source/misc_parameters.h source/frame.h source/setup/plasma.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/type_marker.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_scales.o: source/setup/tag_scales.c source/misc_units.h source/log.h source/misc_cfgReader.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/IO_mesh.o: source/IO_mesh.c source/log.h source/IO_mesh.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/IO_fileMap.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_mesh.o: source/setup/tag_mesh.c source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/main.o: source/core/main.c source/timer.h source/type_reg.h source/type_marker.h source/log.h source/test_mesh.h source/misc_units.h source/misc_MPItags.h source/misc_cfgReader.h source/misc_definedKeys.h source/IO_sys.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/IO_tecplot.h source/IO_fileMap.h source/profiler.h source/spectr_dump.h source/commandLine.h source/core/em.h source/type_mesh.h source/core/plasma.h source/core/plasma_IO.h source/core/plasma_VSP.h source/core/plasma_parallel.h source/core/plasma_gaussTest.h source/core/diag_probe.h source/probe.h source/core/diag_WDensity.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_plasmaWave.o: source/setup/tag_plasmaWave.c source/type_marker.h source/setup/main.h source/setup/plasma.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_parameters.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/plasma.o: source/setup/plasma.c source/type_marker.h source/core/plasma_rho.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_marker.h source/log.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/plasma_gaussTest.o: source/core/plasma_gaussTest.c source/profiler.h source/log.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/misc_cfgReader.h source/core/plasma.h source/type_mesh.h source/misc_partition.h source/type_marker.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_scissors.o: source/setup/tag_scissors.c source/log.h source/misc_units.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/misc_markerPacker.o: source/misc_markerPacker.c source/log.h source/misc_markerPacker.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/plasma_VSP.o: source/core/plasma_VSP.c source/core/plasma_VSP.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/em_caps.o: source/core/em_caps.c source/core/em_TFSF.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/scpic/em_sources.h source/type_reg.h source/type_mesh.h source/scpic/em_mirror.h source/log.h source/parr_meshes.h source/misc_socket.h source/misc_MPItags.h source/parr_regLists.h source/parr_ghostCells.h source/type_mesh.h source/core/emCap_Mur.c source/core/emCap_mirror.c source/core/emCap_periodic.c source/core/emCap_decomposition.c  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/type_reg.o: source/type_reg.c source/type_reg.h source/log.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/plasma_IO.o: source/core/plasma_IO.c source/log.h source/core/plasma.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_marker.h source/IO_names.h source/misc_MPItags.h source/misc_parameters.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/test_mesh.o: source/test_mesh.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/type_vector.o: source/type_vector.c source/type_vector.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/type_mesh.o: source/type_mesh.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/misc_cfgReader.o: source/misc_cfgReader.c source/log.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/tecplot/main.o: source/tecplot/main.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_CFile.h source/log.h source/misc_units.h source/misc_partition.h source/misc_cfgReader.h source/misc_parameters.h source/IO_tecplot.h source/type_mesh.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/em.o: source/core/em.c source/log.h source/profiler.h source/core/em.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/core/em_caps.h source/core/em_TFSF.h source/core/em_gaussSpot.h source/scpic/em_sources.h source/type_reg.h source/type_mesh.h source/scpic/em_mirror.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_DF_uniform.o: source/setup/tag_DF_uniform.c source/setup/plasma.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_marker.h source/setup/tag_plasma.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_partition.h source/misc_parameters.h source/setup/main.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_gaussSpot.o: source/setup/tag_gaussSpot.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/misc_units.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_boundary.o: source/setup/tag_boundary.c source/misc_cfgReader.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/TFSFdefrag/main.o: source/TFSFdefrag/main.c source/type_reg.h source/log.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/parr_regLists.o: source/parr_regLists.c source/misc_MPItags.h source/log.h source/type_reg.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_TFSF.o: source/setup/tag_TFSF.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/core/em_TFSF.h source/type_mesh.h source/setup/main.h source/log.h source/misc_units.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/type_mesh2.o: source/type_mesh2.c source/type_mesh2.h source/misc_parameters.h source/frame.h source/log.h source/type_reg.h source/type_mesh.h source/misc_partition.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_planePulse.o: source/setup/tag_planePulse.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_parameters.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_TFSFFaces.o: source/setup/tag_TFSFFaces.c source/log.h source/misc_cfgReader.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/IO_names.o: source/IO_names.c source/log.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/em_gaussSpot.o: source/core/em_gaussSpot.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_parameters.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_photoelectrons.o: source/setup/tag_photoelectrons.c source/misc_units.h source/log.h source/misc_cfgReader.h source/setup/main.h source/setup/plasma.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_marker.h source/setup/tag_photoelectrons.h source/setup/setup_denavit.h source/setup/setup_distrMapper.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/plasma_rho.o: source/core/plasma_rho.c source/type_marker.h source/misc_PIC.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/profiler.h source/parr_meshes.h source/misc_socket.h source/misc_MPItags.h source/core/plasma_VSP.h source/type_mesh.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/diag_WDensity.o: source/core/diag_WDensity.c source/log.h source/misc_MPItags.h source/misc_cfgReader.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/diag_probe.o: source/core/diag_probe.c source/core/diag_probe.h source/probe.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/misc_PIC.h source/type_mesh.h source/misc_units.h source/log.h source/misc_MPItags.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/IO_fileMap.o: source/IO_fileMap.c source/log.h source/IO_names.h source/IO_fileMap.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/misc_partition.o: source/misc_partition.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/commandLine.h source/misc_cfgReader.h source/misc_partitionPack.c  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/profiler.o: source/profiler.c source/timer.h source/profiler.h source/log.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_EMWave.o: source/setup/tag_EMWave.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/mf_vector.h source/log.h source/setup/main.h source/misc_units.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_units.o: source/setup/tag_units.c source/misc_units.h source/log.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_laser.o: source/setup/tag_laser.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/setup/main.h source/log.h source/misc_units.h source/misc_cfgReader.h source/scpic/em_sources.h source/type_reg.h source/type_mesh.h source/scpic/em_mirror.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_EMResonator.o: source/setup/tag_EMResonator.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/mf_vector.h source/setup/main.h source/log.h source/misc_units.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/plasma.o: source/core/plasma.c source/type_marker.h source/misc_PIC.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/misc_cfgReader.h source/misc_markerPlacer.h source/timer.h source/core/plasma.h source/type_mesh.h source/profiler.h source/parr_meshes.h source/misc_socket.h source/misc_MPItags.h source/core/plasma_rho.h source/core/plasma_VSP.h source/core/plasma_parallel.h source/core/plasma_walls.c source/core/plasma_currentWalls.c source/misc_definedKeys.h source/core/plasma_currentKernel.c  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/mathematica/main.o: source/mathematica/main.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/misc_units.h source/misc_partition.h source/misc_cfgReader.h source/misc_parameters.h source/IO_tecplot.h source/type_mesh.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_cluster.o: source/setup/tag_cluster.c source/log.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_marker.h source/misc_units.h source/misc_cfgReader.h source/setup/main.h source/setup/tag_cluster.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/em_TFSF.o: source/core/em_TFSF.c source/misc_units.h source/log.h source/misc_cfgReader.h source/misc_definedKeys.h source/core/em_TFSF.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_seedWeibel.o: source/setup/tag_seedWeibel.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_marker.h source/log.h source/misc_cfgReader.h source/misc_parameters.h source/setup/tag_photoelectrons.h source/setup/main.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/setup_denavit.o: source/setup/setup_denavit.c source/log.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_ringDF.o: source/setup/tag_ringDF.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_marker.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_parameters.h source/setup/main.h source/setup/plasma.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/log.o: source/log.c source/log.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/IO_sys.o: source/IO_sys.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_marker.h source/IO_sys.h source/IO_mesh.h source/IO_fileMap.h source/IO_names.h source/log.h source/misc_MPItags.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/misc_socket.o: source/misc_socket.c source/log.h source/misc_socket.h source/misc_MPItags.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/type_CFile.o: source/type_CFile.c source/log.h source/type_CFile.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_point.o: source/setup/tag_point.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/misc_units.h source/log.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_maxwell.o: source/setup/tag_maxwell.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/misc_PIC.h source/type_mesh.h source/misc_units.h source/log.h source/misc_cfgReader.h source/setup/main.h source/setup/setup_denavit.h source/setup/plasma.h source/type_marker.h source/setup/tag_plasma.h source/setup/setup_distrMapper.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/parr_ghostCells.o: source/parr_ghostCells.c source/log.h source/parr_meshes.h source/type_reg.h source/misc_socket.h source/misc_MPItags.h source/core/em_caps.h source/type_mesh.h source/misc_partition.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/main.o: source/setup/main.c source/type_marker.h source/misc_units.h source/misc_cfgReader.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/misc_parameters.h source/misc_markerPlacer.h source/log.h source/commandLine.h source/IO_sys.h source/type_mesh.h source/misc_partition.h source/IO_tecplot.h source/IO_names.h source/IO_fileMap.h source/setup/tag.h source/setup/tag_units.h source/setup/tag_mesh.h source/setup/tag_boundary.h source/setup/tag_TFSF.h source/setup/tag_TFSFFaces.h source/setup/tag_laser.h source/setup/tag_plasma.h source/setup/tag_DF_uniform.h source/setup/tag_point.h source/type_mesh.h source/setup/tag_EMWave.h source/setup/tag_ringDF.h source/setup/tag_vShift.h source/setup/tag_photoelectrons.h source/setup/tag_jBalancer.h source/setup/tag_twoStream.h source/setup/tag_gaussSpot.h source/setup/tag_maxwell.h source/setup/tag_planePulse.h source/setup/tag_plasmaWave.h source/setup/tag_EMResonator.h source/setup/tag_scissors.h source/setup/tag_meanVNoise.h source/setup/tag_seedWeibel.h source/setup/tag_seed_PITS.h source/setup/tag_cluster.h source/setup/tag_foil.h source/setup/tag_scales.h source/setup/tag_trianglePrizm.h source/setup/plasma.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/misc_parameters.o: source/misc_parameters.c source/log.h source/misc_cfgReader.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/scpic/em_sources.o: source/scpic/em_sources.c source/scpic/em_sources.h source/type_reg.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/scpic/em_mirror.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_parameters.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/scpic/em_mirror.o: source/scpic/em_mirror.c source/scpic/em_mirror.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_parameters.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/IO_tecplot.o: source/IO_tecplot.c source/IO_mesh.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/IO_fileMap.h source/IO_names.h source/IO_unroll.c source/log.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_seed_PITS.o: source/setup/tag_seed_PITS.c source/misc_units.h source/log.h source/misc_cfgReader.h source/setup/main.h source/setup/tag_photoelectrons.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_plasma.o: source/setup/tag_plasma.c source/type_vector.h source/log.h source/misc_cfgReader.h source/setup/main.h source/setup/plasma.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/type_marker.h source/setup/tag_plasma_box.c source/setup/tag_plasma_convex.c source/setup/tag_plasma_spheres.c source/setup/tag_plasma_cylinders.c  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/misc_markerPlacer.o: source/misc_markerPlacer.c source/misc_PIC.h source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/misc_cfgReader.h source/misc_markerPacker.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/core/plasma_parallel.o: source/core/plasma_parallel.c source/log.h source/misc_MPItags.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/core/plasma.h source/type_mesh.h source/misc_partition.h source/type_marker.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/misc_units.o: source/misc_units.c source/log.h source/misc_units.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/log_test.o: source/log_test.c source/log.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/tag_trianglePrizm.o: source/setup/tag_trianglePrizm.c source/type_marker.h source/log.h source/misc_units.h source/misc_cfgReader.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/misc_parameters.h source/setup/main.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/distiller/main.o: source/distiller/main.c source/misc_units.h source/log.h source/distiller/dll.h source/type_marker.h source/IO_names.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/commandLine.o: source/commandLine.c source/log.h source/type_vector.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/parr_meshes.o: source/parr_meshes.c source/type_mesh.h source/misc_partition.h source/type_reg.h source/misc_parameters.h source/frame.h source/log.h source/parr_meshes.h source/misc_socket.h source/misc_MPItags.h source/parr_regLists.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/distiller/dll.o: source/distiller/dll.c source/type_marker.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/setup/setup_distrMapper.o: source/setup/setup_distrMapper.c source/log.h source/setup/setup_distrMapper.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@

o/spectr_dump.o: source/spectr_dump.c source/type_mesh2.h source/misc_parameters.h source/frame.h source/log.h source/type_reg.h source/type_mesh.h source/misc_partition.h source/type_CFile.h source/misc_cfgReader.h  Makefile
	@echo "Compiling $@ ..."
	@$(CCOMPILER) $(CFLAGS) $< -o $@


core.out: o/core/main.o o/log.o o/test_mesh.o o/misc_units.o o/misc_cfgReader.o o/IO_sys.o o/IO_mesh.o o/IO_names.o o/misc_partition.o o/type_reg.o o/misc_parameters.o o/IO_tecplot.o o/IO_fileMap.o o/profiler.o o/spectr_dump.o o/type_mesh2.o o/type_CFile.o o/commandLine.o o/type_vector.o o/core/em.o o/core/em_caps.o o/parr_meshes.o o/misc_socket.o o/parr_regLists.o o/parr_ghostCells.o o/core/em_TFSF.o o/core/em_gaussSpot.o o/scpic/em_sources.o o/scpic/em_mirror.o o/type_mesh.o o/core/plasma.o o/misc_markerPlacer.o o/misc_markerPacker.o o/core/plasma_rho.o o/core/plasma_IO.o o/core/plasma_VSP.o o/core/plasma_parallel.o o/core/plasma_gaussTest.o o/core/diag_probe.o o/core/diag_WDensity.o 
	@echo "  Linking $@..."
	@$(CCOMPILER) $^ -o $@ $(LFLAGS) -lm 
	@echo "  Executable $@ is ready."
	@echo ""

.PHONY: autoCleanByMegaMake_core.out
autoCleanByMegaMake_core.out:
	@echo "Removing compiled files for project core.out ..."
	@rm -f core.out o/core/main.o o/log.o o/test_mesh.o o/misc_units.o o/misc_cfgReader.o o/IO_sys.o o/IO_mesh.o o/IO_names.o o/misc_partition.o o/type_reg.o o/misc_parameters.o o/IO_tecplot.o o/IO_fileMap.o o/profiler.o o/spectr_dump.o o/type_mesh2.o o/type_CFile.o o/commandLine.o o/type_vector.o o/core/em.o o/core/em_caps.o o/parr_meshes.o o/misc_socket.o o/parr_regLists.o o/parr_ghostCells.o o/core/em_TFSF.o o/core/em_gaussSpot.o o/scpic/em_sources.o o/scpic/em_mirror.o o/type_mesh.o o/core/plasma.o o/misc_markerPlacer.o o/misc_markerPacker.o o/core/plasma_rho.o o/core/plasma_IO.o o/core/plasma_VSP.o o/core/plasma_parallel.o o/core/plasma_gaussTest.o o/core/diag_probe.o o/core/diag_WDensity.o

setup.out: o/setup/main.o o/misc_units.o o/misc_cfgReader.o o/type_reg.o o/misc_parameters.o o/misc_markerPlacer.o o/misc_markerPacker.o o/log.o o/commandLine.o o/type_vector.o o/IO_sys.o o/IO_mesh.o o/misc_partition.o o/IO_tecplot.o o/IO_names.o o/IO_fileMap.o o/setup/tag.o o/setup/tag_units.o o/setup/tag_mesh.o o/setup/tag_boundary.o o/setup/tag_TFSF.o o/core/em_TFSF.o o/setup/main.o o/setup/tag_TFSFFaces.o o/setup/tag_laser.o o/scpic/em_sources.o o/scpic/em_mirror.o o/setup/tag_plasma.o o/setup/tag_DF_uniform.o o/setup/tag_point.o o/type_mesh.o o/setup/tag_EMWave.o o/setup/tag_ringDF.o o/setup/tag_vShift.o o/setup/tag_photoelectrons.o o/setup/setup_denavit.o o/setup/setup_distrMapper.o o/setup/tag_jBalancer.o o/setup/tag_twoStream.o o/setup/tag_gaussSpot.o o/setup/tag_maxwell.o o/setup/tag_planePulse.o o/setup/tag_plasmaWave.o o/setup/tag_EMResonator.o o/setup/tag_scissors.o o/setup/tag_meanVNoise.o o/setup/tag_seedWeibel.o o/setup/tag_seed_PITS.o o/setup/tag_cluster.o o/setup/tag_foil.o o/setup/tag_scales.o o/setup/tag_trianglePrizm.o o/setup/plasma.o o/core/plasma_rho.o o/profiler.o o/parr_meshes.o o/parr_regLists.o o/misc_socket.o o/core/plasma_VSP.o 
	@echo "  Linking $@..."
	@$(CCOMPILER) $^ -o $@ $(LFLAGS) -lm 
	@echo "  Executable $@ is ready."
	@echo ""

.PHONY: autoCleanByMegaMake_setup.out
autoCleanByMegaMake_setup.out:
	@echo "Removing compiled files for project setup.out ..."
	@rm -f setup.out o/setup/main.o o/misc_units.o o/misc_cfgReader.o o/type_reg.o o/misc_parameters.o o/misc_markerPlacer.o o/misc_markerPacker.o o/log.o o/commandLine.o o/type_vector.o o/IO_sys.o o/IO_mesh.o o/misc_partition.o o/IO_tecplot.o o/IO_names.o o/IO_fileMap.o o/setup/tag.o o/setup/tag_units.o o/setup/tag_mesh.o o/setup/tag_boundary.o o/setup/tag_TFSF.o o/core/em_TFSF.o o/setup/main.o o/setup/tag_TFSFFaces.o o/setup/tag_laser.o o/scpic/em_sources.o o/scpic/em_mirror.o o/setup/tag_plasma.o o/setup/tag_DF_uniform.o o/setup/tag_point.o o/type_mesh.o o/setup/tag_EMWave.o o/setup/tag_ringDF.o o/setup/tag_vShift.o o/setup/tag_photoelectrons.o o/setup/setup_denavit.o o/setup/setup_distrMapper.o o/setup/tag_jBalancer.o o/setup/tag_twoStream.o o/setup/tag_gaussSpot.o o/setup/tag_maxwell.o o/setup/tag_planePulse.o o/setup/tag_plasmaWave.o o/setup/tag_EMResonator.o o/setup/tag_scissors.o o/setup/tag_meanVNoise.o o/setup/tag_seedWeibel.o o/setup/tag_seed_PITS.o o/setup/tag_cluster.o o/setup/tag_foil.o o/setup/tag_scales.o o/setup/tag_trianglePrizm.o o/setup/plasma.o o/core/plasma_rho.o o/profiler.o o/parr_meshes.o o/parr_regLists.o o/misc_socket.o o/core/plasma_VSP.o

TFSFdefrag.out: o/TFSFdefrag/main.o o/type_reg.o o/misc_parameters.o o/log.o o/misc_cfgReader.o 
	@echo "  Linking $@..."
	@$(CCOMPILER) $^ -o $@ $(LFLAGS) -lm 
	@echo "  Executable $@ is ready."
	@echo ""

.PHONY: autoCleanByMegaMake_TFSFdefrag.out
autoCleanByMegaMake_TFSFdefrag.out:
	@echo "Removing compiled files for project TFSFdefrag.out ..."
	@rm -f TFSFdefrag.out o/TFSFdefrag/main.o o/type_reg.o o/misc_parameters.o o/log.o o/misc_cfgReader.o

o/distiller.out: o/distiller/main.o o/misc_units.o o/misc_cfgReader.o o/log.o o/distiller/dll.o o/IO_names.o 
	@echo "  Linking $@..."
	@$(CCOMPILER) $^ -o $@ $(LFLAGS) -lm 
	@echo "  Executable $@ is ready."
	@echo ""

.PHONY: autoCleanByMegaMake_o/distiller.out
autoCleanByMegaMake_o/distiller.out:
	@echo "Removing compiled files for project o/distiller.out ..."
	@rm -f o/distiller.out o/distiller/main.o o/misc_units.o o/misc_cfgReader.o o/log.o o/distiller/dll.o o/IO_names.o

math.out: o/mathematica/main.o o/type_reg.o o/log.o o/misc_units.o o/misc_partition.o o/commandLine.o o/type_vector.o o/misc_cfgReader.o o/misc_parameters.o o/IO_tecplot.o o/IO_mesh.o o/IO_fileMap.o o/IO_names.o o/type_mesh.o 
	@echo "  Linking $@..."
	@$(CCOMPILER) $^ -o $@ $(LFLAGS) -lm 
	@echo "  Executable $@ is ready."
	@echo ""

.PHONY: autoCleanByMegaMake_math.out
autoCleanByMegaMake_math.out:
	@echo "Removing compiled files for project math.out ..."
	@rm -f math.out o/mathematica/main.o o/type_reg.o o/log.o o/misc_units.o o/misc_partition.o o/commandLine.o o/type_vector.o o/misc_cfgReader.o o/misc_parameters.o o/IO_tecplot.o o/IO_mesh.o o/IO_fileMap.o o/IO_names.o o/type_mesh.o

test_logger.out: o/log_test.o o/log.o 
	@echo "  Linking $@..."
	@$(CCOMPILER) $^ -o $@ $(LFLAGS) -lm 
	@echo "  Executable $@ is ready."
	@echo ""

.PHONY: autoCleanByMegaMake_test_logger.out
autoCleanByMegaMake_test_logger.out:
	@echo "Removing compiled files for project test_logger.out ..."
	@rm -f test_logger.out o/log_test.o o/log.o

tecplot.out: o/tecplot/main.o o/type_reg.o o/type_CFile.o o/log.o o/misc_units.o o/misc_partition.o o/commandLine.o o/type_vector.o o/misc_cfgReader.o o/misc_parameters.o o/IO_tecplot.o o/IO_mesh.o o/IO_fileMap.o o/IO_names.o o/type_mesh.o 
	@echo "  Linking $@..."
	@$(CCOMPILER) $^ -o $@ $(LFLAGS) -lm 
	@echo "  Executable $@ is ready."
	@echo ""

.PHONY: autoCleanByMegaMake_tecplot.out
autoCleanByMegaMake_tecplot.out:
	@echo "Removing compiled files for project tecplot.out ..."
	@rm -f tecplot.out o/tecplot/main.o o/type_reg.o o/type_CFile.o o/log.o o/misc_units.o o/misc_partition.o o/commandLine.o o/type_vector.o o/misc_cfgReader.o o/misc_parameters.o o/IO_tecplot.o o/IO_mesh.o o/IO_fileMap.o o/IO_names.o o/type_mesh.o

.PHONY: clean
clean: autoCleanByMegaMake_core.out autoCleanByMegaMake_setup.out autoCleanByMegaMake_TFSFdefrag.out autoCleanByMegaMake_o/distiller.out autoCleanByMegaMake_math.out autoCleanByMegaMake_test_logger.out autoCleanByMegaMake_tecplot.out 


all: core.out setup.out TFSFdefrag.out o/distiller.out math.out test_logger.out tecplot.out 
	@echo ;
	@echo "All files are compiled using '$(CCOMPILER) $(CFLAGS)'."
	@echo ;

