#!/usr/bin/python3

# Filter out read pairs where there are inconsistent quality score lengths in at least one pair

import sys
import fileinput

if len(sys.argv) < 2:
	print("Usage: " + sys.argv[0] + " read1.fq read2.fq")
	sys.exit(1)

with open(sys.argv[1]) as fq1, open(sys.argv[2]) as fq2:
	line_counter = 0
	for line1, line2 in zip(fq1,fq2):
		line1 = line1.strip()
		line2 = line2.strip()


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

