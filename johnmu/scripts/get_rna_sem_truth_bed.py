#!/usr/bin/python

# Reads in .ti and .sim.isoforms.results from RNA-SEM simulation 
# and outputs a bed file of all the junctions in transcripts with 
# at least one read count. This is not exactly the truth set of junctions
# but it is close enough. 

import sys
import argparse

parser = argparse.ArgumentParser(sys.argv[0], formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--ti", help="Minimum BED region length", type=argparse.FileType('r'), required=True)
parser.add_argument("--isoform", help="The .sim.isoforms.results file", type=argparse.FileType('r'), required=True)

args = parser.parse_args()

# read in isoform file and store in dictionary
good_isoforms = set()
first_line = True # skip the header line
for line in args.isoform:
	if first_line:
		first_line = False
		continue
	
	line = line.strip()
	if len(line) == 0:
		continue
	
	ll = line.split()
	
	if float(ll[3]) > 0.0:
		good_isoforms.add(ll[0]); # Add the transcript name
	

# go through ti file and output junctions to standard out
first_line = True # skip the header line
line_idx = 0
chr = None
output_bed = False
for line in args.ti:
	if first_line:
		first_line = False
		continue
		
	line = line.strip()
	if len(line) == 0:
		continue
		
	if line_idx % 6 == 0 and line in good_isoforms:
		output_bed = True
	
	if line_idx % 6 == 2:
		chr = line
		
	if line_idx % 6 == 4 and output_bed:
		ll = line.split()
		if len(ll) > 3:
			start = [int(v) for v in ll[2::2]]
			end = [int(v) for v in ll[3::2]]
			for i in xrange(0,len(end)):
				print '\t'.join([chr, str(start[i]), str(end[i]-1)])
		
	
	if line_idx % 6 == 5:
		chr = None
		output_bed = False
	
	line_idx += 1

