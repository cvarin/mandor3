#!/bin/sh

mn_workArea.sh

if [ -e Makefile ]
then
  make cleanAll
  rm -f interface.*.input test_rec.sh test_play.sh test_run.sh test_clean.sh
  rm -f -r res
else
  echo "In this folder Makefile does not exist."
  echo "Probably you run command in the test folder (which is not wise)."
fi


