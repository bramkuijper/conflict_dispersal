#!/usr/bin/env python3

# plot the battleground between offspring and parent strategies

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams

rcParams['text.usetex'] = True
rcParams['font.family'] = 'sans-serif'

# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]  

data = pd.read_csv("total_output.csv", sep=";")
subset = data[(data["d"] == 0.05) & (data["sigma21"] == 0.5) & (data["sigma12"] == 0.05)]

subset_off = subset[subset["type"] == "offspring"][["c1","c2","s1"]]
subset_mother = subset[subset["type"] == "mother"][["c1","c2","s1"]]

# make a pivot table
pivot_mom = subset_mother.pivot_table(values="s1", index="c1", columns="c2")
pivot_off = subset_off.pivot_table(values="s1", index="c1", columns="c2")

xo = pivot_mom.columns.values
yo = pivot_mom.index.values

x, y = np.meshgrid(xo, yo)
z = pivot_mom.values
z2 = pivot_off.values

fig = plt.figure()

ax = fig.add_subplot(111, projection="3d")
ax.plot_wireframe(x, y, z, color="red")
ax.plot_wireframe(x, y, z2, color="blue")
ax.set_xlabel(r'Survival cost envt 2, $c_{2}$')
ax.set_ylabel(r'Survival cost envt 1, $c_{1}$')
ax.set_zlabel(r'Proportion $z_{1}$ offspring envt 1')

format = "pdf"
plt.savefig("battleground3d_s1" + format, format=format)

subset_off = subset[subset["type"] == "offspring"][["c1","c2","s2"]]
subset_mother = subset[subset["type"] == "mother"][["c1","c2","s2"]]

# make a pivot table
pivot_mom = subset_mother.pivot_table(values="s2", index="c1", columns="c2")
pivot_off = subset_off.pivot_table(values="s2", index="c1", columns="c2")

xo = pivot_mom.columns.values
yo = pivot_mom.index.values

x, y = np.meshgrid(xo, yo)
z = pivot_mom.values
z2 = pivot_off.values





fig = plt.figure()

ax = fig.add_subplot(111, projection="3d")
ax.plot_wireframe(x, y, z2, color="blue")
ax.plot_wireframe(x, y, z, color="red")
ax.set_xlabel(r'Survival cost envt 2, $c_{2}$')
ax.set_ylabel(r'Survival cost envt 1, $c_{1}$')
ax.set_zlabel(r'Proportion $z_{1}$ offspring envt 2')

format = "pdf"
plt.savefig("battleground3d_s2" + format, format=format)
