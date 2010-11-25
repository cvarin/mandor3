#!/bin/sh

ncpu=1								# Set number of cpus here.

echo "Start to play recorded signal."
echo;

cp -f play.run_mandor.cfg run_mandor.cfg
cp -f play.run_probes.cfg run_probes.cfg

.mn_set.sh laser.play.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out

vTec2D.v2.sh
