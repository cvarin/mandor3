#!/bin/sh

mn_workArea.sh 

if [ -e Makefile ]
then
  make cleanAll
  rm -f pulseZ.input test_run.sh test_clean.sh
else
  echo "In this folder Makefile does not exist."
  echo "Probably you run command in the test folder (which is not wise)."
fi


