#!/bin/sh

ncpu=1

.mn_set.sh twoStream.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out
vPhase.sh
vW.sh
