#!/usr/bin/env python

n = [ 2, 4 ]
d = [ 0.05, 0.1, 0.25, 0.5, 0.75]
sigma12 = [ 0.05, 0.25, 0.5, 0.75 ]
sigma21 = [ 0.05, 0.25, 0.5, 0.75 ]

ctr = 0

exe = "xresolution_numerical"

for n_i in n:
    for d_i in d:
        for sigma12_i in sigma12:
            for sigma21_i in sigma21:

                p1 = sigma21_i / (sigma12_i + sigma21_i)
                print("echo " + str(ctr))
                ctr+=1
                print(exe + " " + str(p1) + " " + str(n_i) + " " + str(d_i) + " " + str(sigma_12_i) + " " + str(sigma_21_i))
