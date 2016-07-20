#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.cm as cm
import  matplotlib.pyplot as plt
import math

import sys, csv

def readfile(filename):
    fo = csv.DictReader(open(filename,"r"), delimiter=";")
    data = [[], [], [], []]

    histdict = { "x": [],"s0": [],"s1": [],"dS": [], "dNS": [] }

    generation = 0 

    bins = np.arange(0,1.01,0.01)
    binmids = np.arange(0,1.01,0.01)

    print "generation;val;s0;s1;dS;dNS;"

    for row in fo:
        data[0].append(float(row["s0"]))
        data[1].append(float(row["s1"]))
        data[2].append(float(row["dS"]))
        data[3].append(float(row["dNS"]))

        if generation != row["generation"]:
            histdict["s0"].append([ math.log10(x+1) for x in list(np.histogram(data[0],bins=bins)[0])])
            histdict["s1"].append([ math.log10(x+1) for x in list(np.histogram(data[1],bins=bins)[0])])
            histdict["dS"].append([ math.log10(x+1) for x in list(np.histogram(data[2],bins=bins)[0])])
            histdict["dNS"].append([ math.log10(x+1) for x in list(np.histogram(data[3],bins=bins)[0])])

            data = [[],[],[],[],[]]
            generation = row["generation"]
            #for i in range(0,(len(histdict["s0"][-1]))):
            #    print str(generation) + ";" + str(bins[i]) + ";" + str(histdict["s0"][-1][i]) + ";" + str(histdict["s1"][-1][i]) + ";" + str(histdict["dS"][-1][i]) + ";" + str(histdict["dNS"][-1][i]) + ";"

            #    if bins[i] == 0.99:
            #        print str(generation) + ";1.0;0;0;0;0;" 


#            print generation

#            if int(generation) == 1000:
#                break;


    return(histdict)

histfile = sys.argv[1]

data = readfile(histfile)

nrow = 4

colormapje = cm.jet

fig = plt.figure(figsize=(8,8))
fig.subplots_adjust(bottom=0.1,top=0.9)
ax = fig.add_subplot(1,nrow,1)
cax = ax.imshow(data["s0"], cmap=colormapje)
ax.set_title("$s_{0}$")

ax = fig.add_subplot(1,nrow,2)
cax = ax.imshow(data["s1"],cmap=colormapje)
ax.set_title("$s_{1}$")

ax = fig.add_subplot(1,nrow,3)
cax = ax.imshow(data["dS"],cmap=colormapje)
ax.set_title("$d_{\mathrm{S}}$")

ax = fig.add_subplot(1,nrow,4)
cax = ax.imshow(data["dNS"],cmap=colormapje)
ax.set_title("$d_{\mathrm{NS}}$")

plt.savefig("branchplot_" + histfile + ".pdf",format="pdf")
