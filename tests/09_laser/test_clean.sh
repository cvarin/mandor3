#!/bin/sh

mn_workArea.sh

if [ -e Makefile ]
then
  make cleanAll
  rm -f laser.*.input test_rec.sh test_play.sh test_rec2.sh test_play2.sh test_run.sh test_clean.sh ./lay/tecplot2D.rec.lay ./lay/tecplot2D.play.lay
  rm -f res play.run_*.cfg rec.run_*.cfg
else
  echo "In this folder Makefile does not exist."
  echo "Probably you run command in the test folder (which is not wise)."
fi


