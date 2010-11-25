#!/bin/sh

ncpu=4

echo "Start to play recorded signal."
echo;

cp -f play.run_mandor.cfg run_mandor.cfg

.mn_set.sh 2mirror.play.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out

vTec2D.v2.sh
