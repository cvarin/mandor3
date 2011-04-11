#!/bin/sh

echo "----------------------------------------------------"
echo "Please activate core_continuity test by using"
echo " -DMC_GAUSS_TEST flag in Makefile!"
echo "----------------------------------------------------"

ncpu=1

# Quick erase to speed-up 'make cleanAll'.
rm -f output/markers/*.*
rm -f binData/*.*
make cleanAll

make -j 6 && mpirun -np $ncpu ./setup.out slab.input create
mpirun -np $ncpu ./core.out

./tracer markers.csv
./core.out
./distiller -a
./tracer -v
