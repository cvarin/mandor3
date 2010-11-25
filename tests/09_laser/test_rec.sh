#!/bin/sh

ncpu=4								# Set number of cpus here.

echo "Start to record signal."
echo;

cp -f rec.run_mandor.cfg run_mandor.cfg
cp -f rec.run_probes.cfg run_probes.cfg

.mn_set.sh laser.rec.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out

./TFSFdefrag.out

vTec2D.v2.sh
