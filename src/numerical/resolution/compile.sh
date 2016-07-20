#!/usr/bin/env bash

./generate_cpp.py

g++ -Wall -o xresolution_numerical resolution_numerical2.cpp -lm -lrt -lgsl -lgslcblas

