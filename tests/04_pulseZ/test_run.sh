#!/bin/sh

ncpu=10
.mn_set.sh pulseZ.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out
vTec1D.v2.sh
vW.sh
