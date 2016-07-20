#!/usr/bin/env python3

import numpy as np


s12 = [ 0, 0.05, 0.25, 0.5, 0.75, 1.0 ]
s21 = [ 0, 0.05, 0.25, 0.5, 0.75, 1.0 ]
d = [ 0, 0.05, 0.25, 0.5, 0.75, 1.0 ]
n = [ 2] 

ctr = 0

for s12_i in s12:
    for s21_i in s21:
        for d_i in d:
            for n_i in n:
                print("echo " + str(ctr))
                ctr += 1
                print("./xconflict_numerical  0.5 " + str(n_i) + " " + str(d_i) + " " + str(s12_i) + " "  +str(s21_i))
