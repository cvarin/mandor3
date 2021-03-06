#!/bin/sh

echo "----------------------------------------------------"
echo "Please activate core_continuity test by using"
echo " -DMC_GAUSS_TEST flag in Makefile!"
echo "----------------------------------------------------"

ncpu=8

make -j 6 && mpirun -np $ncpu ./setup.out slab.input create
mpirun -np $ncpu ./core.out

cd visualization
alchemist gauss_simple.alch
