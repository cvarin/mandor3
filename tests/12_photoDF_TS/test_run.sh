#!/bin/sh

ncpu=1

.mn_set.sh photoDF_TS.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out
vW.sh
vSpectr.sh
