#!/bin/sh

ncpu=4								# Set number of cpus here.

.mn_set.sh interface.play.input "mpirun -np $ncpu"
mpirun -np $ncpu ./core.out
vProbes.sh

echo "Please save output of the probe diagnostic using File/Export"
echo "menu of the TecPlot to compare with the result of record session."

