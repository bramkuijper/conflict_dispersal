#!/usr/bin/env python

import sys, re, os.path, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# read in the csv file
dictdat = csv.DictReader(open(filename,"r"), delimiter=";")


# data can now only accessed through looping row by row
# whereas we want lists of each column
# this function does that
def get_csvdict_by_column(the_raw_dict):

    # initialize a empty dict to contain the data
    by_column_dat = {}

    rowctr = 0

    # loop through the rows of the csv file and
    # put data in the dictionary
    for row in dictdat:

        if rowctr == 0:
            for key in row.keys():
                by_column_dat[key] = []
        
        rowctr += 1

        if None in row.values():
            break

        for key, val in row.iteritems():
            if key != "":
                by_column_dat[key].append(float(val))


    return by_column_dat

histdat = get_csvdict_by_column(dictdat)

# generate the figure

# initialize and specify size 
fig = plt.figure(figsize=(10,18))

num_rows = 5

# add first subplot
plt.subplot(num_rows,1,1)
plt.plot(histdat["generation"],histdat["s0"],'r',histdat["generation"],histdat["s1"],'b',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'maternal signal, $\bar{s}_{i}$')
plt.legend((r'$\bar{s}_{0}$',r'$\bar{s}_{1}$'))
plt.ylim(-0.05,1.05)

plt.subplot(num_rows,1,2)
plt.plot(histdat["generation"],histdat["dsignal"],'g',histdat["generation"],histdat["dnosignal"],'#0EB2FF',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'offspring dispersal response, $\bar{d}_{j}$')
plt.legend((r'$\bar{d}_{\mathrm{S}}$',r'$\bar{d}_{\mathrm{NS}}$'))
plt.ylim(-0.05,1.05)

plt.subplot(num_rows,1,3)
plt.plot(histdat["generation"],histdat["d0"],'#FF6700',histdat["generation"],histdat["d1"],'#C13CFF',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'resulting dispersal, $\bar{d}_{j}$')
plt.legend((r'$\bar{d}_{\mathrm{0}}$',r'$\bar{d}_{\mathrm{1}}$'))
plt.ylim(-0.05,1.05)

plt.subplot(num_rows,1,4)
plt.plot(histdat["generation"],histdat["vars0"],'r',histdat["generation"],histdat["vars1"],'b',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'var maternal signal, $\sigma_{{s}_{i}}^{2}$')
plt.legend((r'$\sigma_{{s}_{0}}^{2}$',r'$\sigma_{{s}_{1}}^{2}$'))

plt.subplot(num_rows,1,5)
plt.plot(histdat["generation"],histdat["vardsignal"],'g',histdat["generation"],histdat["vardnosignal"],'#0EB2FF',linewidth=1)
plt.tick_params(axis='x',which='both',bottom='on',top='on')
plt.ylabel(r'var maternal signal, $\sigma_{{d}_{i}}^{2}$')
plt.legend((r'$\sigma_{{d}_{\mathrm{S}}}^{2}$',r'$\sigma_{{d}_{\mathrm{NS}}}^{2}$'))

#plt.subplot(num_rows,1,2)
#plt.plot(histdat["time"],histdat["p00"],'g',histdat["time"],histdat["p01"],'y',histdat["time"],histdat["p10"],'r',histdat["time"],histdat["p11"],('#8900ff'))
#plt.tick_params(axis='x',which='both',bottom='on',top='on',labelbottom='off')
#plt.ylabel(r'$\bar{p}_{ij}$')
#plt.legend((r'$p_{q_{1},e_{1}}$',r'$p_{q_{1},e_{2}}$',r'$p_{q_{2},e_{1}}$',r'$p_{q_{2},e_{2}}$'))
#plt.ylim(-0.05,1.05)
#
#plt.subplot(num_rows,1,3)
#plt.plot(histdat["time"],histdat["v00"],'g',histdat["time"],histdat["v01"],'y',histdat["time"],histdat["v10"],'r',histdat["time"],histdat["v11"],('#8900ff'))
#plt.ylabel(r'$\bar{v}_{ij}$')
#plt.legend((r'$v_{q_{1},e_{1}}$',r'$v_{q_{1},e_{2}}$',r'$v_{q_{2},e_{1}}$',r'$v_{q_{2},e_{2}}$'))
#plt.ylim(-0.05,1.05)


graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf",)
