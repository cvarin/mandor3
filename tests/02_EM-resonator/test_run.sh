#!/bin/sh

ncpu=2

cfg=cube
cfg=ABC
cfg=resonator

.mn_set.sh ${cfg}.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out

#./zzz.sh
vW.sh
#vTec3D.v2.sh
#vPhase.sh
#vChargeConserv.sh
vProbes.sh
