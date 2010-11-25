#!/bin/sh

NCPU=2

./recorder && ./laser.rec_all.sh 

mpirun -np $NCPU ./setup.out laser.play.input create && mpirun -np $NCPU ./core.out

cd visualization && alchemist EM_2D.alch && cd ..
