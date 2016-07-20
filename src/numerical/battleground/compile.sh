#!/usr/bin/env bash

./generate_cpp.py

i=""
l=""

if [ "$HOSTNAME" = carson ]; then 
    i=-I/cm/shared/apps/gsl/gcc/1.16/include
    j=-L/cm/shared/apps/gsl/gcc/1.16/lib
fi

g++ -Wall $i $j -o xconflict_numerical conflict_numerical2.cpp -lm -lrt -lgsl -lgslcblas

