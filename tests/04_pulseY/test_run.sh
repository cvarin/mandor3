#!/bin/sh

ncpu=10
.mn_set.sh pulseY.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out
vTec1D.v2.sh
vW.sh
