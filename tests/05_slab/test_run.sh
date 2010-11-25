#!/bin/sh

echo "----------------------------------------------------"
echo "Please activate core_continuity test by using"
echo " -DMC_GAUSS_TEST flag in Makefile!"
echo "----------------------------------------------------"

ncpu=8

.mn_set.sh slab.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out

cd visualization
alchemist gauss_simple.alch
