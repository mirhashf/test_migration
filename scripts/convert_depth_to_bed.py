#!/usr/bin/python

# Reads a depth file from samtools depth 
# depth file position is 1-based and sorted
# <chr> <pos> <depth>
# and turns it into a BED file of the form
# <chr> <start_pos> <end_pos> <depth>
# remember that BED is 0-based and end position is not inclusive

import fileinput
import sys


def print_bed(chr, start, end, depth):
    print '\t'.join([chr, str(start), str(end), str(depth)])


curr_chr = None
curr_start = None
curr_end = None
curr_depth = None
for line in fileinput.input(sys.argv[1:]):
    line = line.strip()
    if len(line) == 0:
        continue

    ll = line.split()
    chrname = ll[0]
    pos = int(ll[1])
    depth = int(ll[2])

    if not curr_chr:
        curr_chr = chrname

    if curr_chr != chrname:
        print_bed(curr_chr, curr_start, curr_end, curr_depth)
        curr_chr = chrname
        curr_start = None
        curr_end = None
        curr_depth = None

    if not curr_depth:
        curr_depth = depth

    if curr_depth != depth:
        print_bed(curr_chr, curr_start, curr_end, curr_depth)
        curr_start = None
        curr_end = None
        curr_depth = depth

    if not curr_start:
        curr_start = pos - 1
        curr_end = pos
    else:
        if pos > curr_end + 1:
            print_bed(curr_chr, curr_start, curr_end, curr_depth)
        else:
            curr_end += 1

if curr_chr:
    # if we have at least one line print it or print the last one
    print_bed(curr_chr, curr_start, curr_end, curr_depth)
