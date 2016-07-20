#!/usr/bin/env bash

./generate_cpp.py

g++ -Wall -o xconflict_numerical conflict_numerical2.cpp -lm -lrt -lgsl -lgslcblas

