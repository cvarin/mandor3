#!/bin/sh

ncpu=2

echo "Start to record signal."
echo;

cp -f rec.run_mandor.cfg run_mandor.cfg

.mn_set.sh 2mirror.rec.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out

./TFSFdefrag.out

vTec2D.v2.sh
