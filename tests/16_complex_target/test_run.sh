#!/bin/sh

ncpu=1

make cleanData
make -j4

mpirun -np $ncpu ./setup.out cylinders.input create
mpirun -np $ncpu ./core.out

#vProbes.sh
# vW.sh

cd visualization
alchemist EM_3D.alch

#cd visualization
#alchemist dbg_parallel.alch

