#!/usr/bin/python

# Formats TopHat bed file to the actual junctions

import sys
import fileinput

if len(sys.argv) != 2:
	print "Usage: " + sys.argv[0] + " tophat_junctions.bed"
	sys.exit(1)

first_line = True # skip the header line
for line in fileinput.input(sys.argv[1]):
	if first_line:
		first_line = False
		continue

	line = line.strip()
	if len(line) == 0:
		continue

	ll = line.split()

	chr = ll[0]
	start = int(ll[1])
	end = int(ll[2])
	lens = [int(i) for i in ll[10].split(',')]

	print '\t'.join([chr,str(start + lens[0]),str(end - lens[1])])

