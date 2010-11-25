#!/bin/sh

ncpu=1						# Serial
#ncpu=4	 					# Parallel

.mn_set.sh cube.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out

vProbes.sh
vW.sh
