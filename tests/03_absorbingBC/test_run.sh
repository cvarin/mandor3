#!/bin/sh

ncpu=1

.mn_set.sh absorbingBC.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out

vProbes.sh
vW.sh
