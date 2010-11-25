#!/bin/sh

until [ -z "$1" ]  # Until all parameters used up...do
do
  cp -f $1 ./tmpFile1
  echo "  " > ./tmpFile2
  
  while [ `cmp ./tmpFile1 ./tmpFile2 | wc -l ` -ne 0 ]
  do
    cat ./tmpFile1 | sed 's/mc_bound_periodic/BC_PERIODIC/g' > ./tmpFile2
    cat ./tmpFile2 | sed 's/mc_bound_periodic/BC_PERIODIC/g' > ./tmpFile1
    cat ./tmpFile1 | sed 's/mc_bound_mirror/BC_MIRROR/g' > ./tmpFile2
    cat ./tmpFile2 | sed 's/mc_bound_mirror/BC_MIRROR/g' > ./tmpFile1
    cat ./tmpFile1 | sed 's/mc_bound_openMur/BC_OPEN/g' > ./tmpFile2
    cat ./tmpFile2 | sed 's/mc_bound_openMur/BC_OPEN/g' > ./tmpFile1
    cat ./tmpFile1 | sed 's/mc_bound_splitter/BC_SPLITTER/g' > ./tmpFile2
    cat ./tmpFile2 | sed 's/mc_bound_splitter/BC_SPLITTER/g' > ./tmpFile1
  done
  
  mv -f ./tmpFile1 ./$1
  rm -f ./tmpFile2
  
  shift
done
