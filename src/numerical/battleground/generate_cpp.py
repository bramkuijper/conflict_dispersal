#!/usr/bin/env python
# -*- coding: utf-8 -*-
# script to get mathematica output into the *.cpp file

import os, re


# regular expressions to replace variables 
regexps = [ (r"s\((\d)\)", "s\\1"), (r"c\((\d)\)", "c\\1"), (r"p\((\d)\)", "p\\1"), (r"([u|v])\((\d)\)", "\\1\\2"), (r"σ\((\d),(\d)\)", "sigma\\1\\2"), ("Power","pow"),("Sqrt","sqrt"), (r"λ","ev"), (r"List","") ]

# open base *.cpp file
f = open("conflict_numerical.cpp")
fl = f.read()
f.close()

# open eigenvectors file
f = open("eigenvectors.txt")
fl_eigenvec = f.read()
f.close()

# selection gradients mom
f = open("selgrad_mother.txt")
fl_selgrad_mom = f.read()
f.close()

# selection gradients offspring
f = open("selgrad_offspring.txt")
fl_selgrad_off = f.read()
f.close()

fl = re.sub("EIGENVECS",fl_eigenvec, fl) 
fl = re.sub("SELGRADS_MOM",fl_selgrad_mom, fl) 
fl = re.sub("SELGRADS_OFF",fl_selgrad_off, fl) 

# change all indices etc to a c-able format
for regexp in regexps:
    fl = re.sub(regexp[0], regexp[1], fl)

f = open("conflict_numerical2.cpp","w")
f.write(fl)
f.close()
