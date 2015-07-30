#!/usr/bin/python

# Reads a depth file from samtools depth 
# depth file position is 1-based and sorted
# <chr> <pos> <depth>
# and turns it into a BED file of the form
# <chr> <start_pos> <end_pos> <depth>
# remember that BED is 0-based and end position is not inclusive
#
# Note: There are some modifications here to handle pegah's specific depth file

import fileinput
import sys
import argparse

parser = argparse.ArgumentParser("Convert depth file to BED",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--minlen", help="Minimum BED region length", default=200, type=int, required=False)
parser.add_argument("--merge_dist", help="Maximum distance to merge intervals", default=200, type=int, required=False)
parser.add_argument("--mincov", help="Minimum coverage", default=2, type=int, required=False)
parser.add_argument("files", help="Depth file(s), if none given will read from stdin", nargs='*')

args = parser.parse_args()

def print_bed(chrname_str, start, end, depth, min_len, min_cov):
    if start and end - start >= min_len and depth >= min_cov:
        print '\t'.join([chrname_str, str(start), str(end), str(depth)])


curr_chr = None
curr_start = None
curr_end = None
curr_depth = None
for line in fileinput.input(args.files):
    line = line.strip()
    if len(line) == 0:
        continue

    ll = line.split()
    chrname = ll[0]
    pos = int(ll[1])
    depth_ll = map(int, ll[2].split(','))
    if max(depth_ll[1:]) > 0:
        # this is set if any of pegah's other features are non-zero
        depth = 0
    else:
        depth = depth_ll[0]

    # Make this optional
    if depth < args.mincov:
        depth = 0
    else:
        depth = args.mincov

    if not curr_chr:
        curr_chr = chrname

    if curr_chr != chrname:
        print_bed(curr_chr, curr_start, curr_end, curr_depth, args.minlen, args.mincov)
        curr_chr = chrname
        curr_start = None
        curr_end = None
        curr_depth = None

    if not curr_depth:
        curr_depth = depth

    if curr_depth != depth:
        print_bed(curr_chr, curr_start, curr_end, curr_depth, args.minlen, args.mincov)
        curr_start = None
        curr_end = None
        curr_depth = depth

    if not curr_start:
        if depth >= args.mincov:
            curr_start = pos - 1
            curr_end = pos
    else:
        if pos > curr_end + 1 + args.merge_dist:
            print_bed(curr_chr, curr_start, curr_end, curr_depth, args.minlen, args.mincov)
            curr_start = pos - 1
            curr_end = pos
        else:
            curr_end = pos

if curr_chr and curr_start:
    # if we have at least one line print it or print the last one
    print_bed(curr_chr, curr_start, curr_end, curr_depth,args.minlen, args.mincov)
