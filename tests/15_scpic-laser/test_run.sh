#!/bin/sh

ncpu=1

make cleanData
make -j4

mpirun -np $ncpu ./setup.out min_workable.input create
mpirun -np $ncpu ./core.out

#vProbes.sh
lm 1 4
# vW.sh

#cd visualization
#alchemist dbg_parallel.alch

