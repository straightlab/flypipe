#!/usr/bin/env python

import sys
import os
import argparse
parser = argparse.ArgumentParser(description="""
Takes in two sorted files and two column values (one for each file) and does a line by line comparison based on the value in the selected columns. Both input files must be sorted by the values in the respective selected column in ascending order. The output is full lines of the first file that don't have the same value in the selected column of the second file. 

sortedfilter.py
ChAR-seq codebase, Straight/Greenleaf/Skotheim Labs
Viviana Risca, vrisca@stanford.edu
Version 1.0, September 2016
""")

# Read in and parse arguments
parser.add_argument('-i', help="input file to be filtered - MUST BE SORTED BY KEY COLUMN")
parser.add_argument('-v', help="file that will be used as the filter - MUST BE SORTED BY KEY COLUMN")
parser.add_argument('-c', type=int, help="key column in file to be filtered (1-based)")
parser.add_argument('-d', type=int, default=1, help="key column in filtering file (1-based), default=1")
args = parser.parse_args(sys.argv[1:])

# Prepare files and variables
ffilt = open(args.v, 'r')
# fin opened below

# Iterate through files, writing to stdout the lines of the input file 
# that don't have a selected column value matching the selected column 
# value in a line of the filter file.

filtline = ffilt.readline()
filtcols = filtline.split()
currentpat = filtcols[args.d - 1]
morefiltpatterns = True
with open(args.i, 'r') as fin:
    for line in fin:
        cols = line.split()
        tocheck = cols[args.c - 1]
        # Debug line:
        # print("tocheck:", tocheck, "currentpat:", currentpat, "morefiltpat:", morefiltpatterns)
        if ((tocheck == currentpat) and morefiltpatterns):
            # if there's a match, advance to the next line in the filter file
            filtline = ffilt.readline()
            if filtline:
                filtcols = filtline.split()
                currentpat = filtcols[args.d - 1]
            else:
                morefiltpatterns = False
        elif ((tocheck > currentpat) and morefiltpatterns):
            # special case where an entry is in the filtering file but not the input file to be filtered
            sys.stdout.write(line)
            filtline = ffilt.readline()
            if filtline:
                filtcols = filtline.split()
                currentpat = filtcols[args.d - 1]
            else:
                morefiltpatterns = False
        else:            
            # if there's no match, write the input file line to stdout
            sys.stdout.write(line)
        prevcheck = tocheck

# Clean up
ffilt.close()
