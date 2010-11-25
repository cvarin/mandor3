#!/bin/sh

ncpu=4								# Set number of cpus here.

echo "Start to record signal."
echo;

.mn_set.sh interface.rec.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out
vProbes.sh

./TFSFdefrag.out

echo "Please save output of the probe diagnostic using File/Export"
echo "menu of the TecPlot (you will compare it later with reproduction."

