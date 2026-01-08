#!/usr/bin/env pypy3
import sys

outfile = open('positions.tsv','w')

with open(sys.argv[1], 'r') as infile:
    for line in infile:
        fields = line.split('\t')
        bases = fields[4].replace(',','.') #make all matching fields the same symbol
        # dont make a single base a variable position
        if len(bases) < 2:
            continue
        # if there is more than 49% support for the reference, not a variable position
        # 1 out of 2 -> not variable
        # 2 out of 3 -> variable
        # double check that!
        if sum(1 for x in bases if x == '.') / len(bases) > 0.49:
            continue

        print(fields[1], file=outfile)

outfile.close()