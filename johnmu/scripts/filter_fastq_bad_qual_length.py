#!/usr/bin/python3

# Filter out read pairs where there are inconsistent quality score lengths in at least one pair

import sys

if len(sys.argv) < 2:
	print("Usage: " + sys.argv[0] + " read1.fq read2.fq")
	sys.exit(1)

with open(sys.argv[1]) as fq1, open(sys.argv[2]) as fq2, open("fixed_"+sys.argv[1],'w') as fixed_fq1, open("fixed_"+sys.argv[2],'w') as fixed_fq2:
	line_counter = 0
	line_buffer = list()
	for line1, line2 in zip(fq1,fq2):
		line1 = line1.strip()
		line2 = line2.strip()
		line_buffer.append((line1,line2))
		line_counter += 1
		if line_counter == 4:
			# Check that the line lengths match
			good_line = True
			if len(line_buffer[1][0]) != len(line_buffer[3][0]):
				good_line = False
			if len(line_buffer[1][1]) != len(line_buffer[3][1]):
				good_line = False
			if good_line:
				# output the files
				for rec in line_buffer:
					fixed_fq1.write(rec[0] + '\n')
					fixed_fq2.write(rec[1] + '\n')
			line_counter = 0
			del line_buffer[:]
