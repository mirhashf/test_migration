#!/usr/bin/env python
'''Converts a PCR-validated variants text file to a BED file, appropriately
extending each region by the requested number of bases.

Usage: python pcr_to_bed.py [pcr_text_file] --output regions.bed --size 100
'''

from optparse import OptionParser
import sys

DEFAULT_REGION_SIZE = 100

def process_pcr_file(input_file, output_file = None, size = DEFAULT_REGION_SIZE):

    print "Processing file '" + input_file + "'."
    
    f = open(input_file, 'r')
    if output_file: o = open(output_file, 'w')
    
    # Ignore first line
    f.readline()
    
    for line in f:
        entries = line.split('\t')
        
        name = entries[0]
        contig = "chr" + entries[1]
        location = int(entries[2])
        region_start = location - size
        region_end = location + size
        
        bed_entry = '\t'.join([contig, str(region_start), str(region_end), name])
        
        if output_file:
            o.write(bed_entry + '\n')
        else:
            print bed_entry
        
    f.close()
    if output_file: o.close()

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-o", "--output", dest="output",
        help="destination BED file", metavar="FILE",
        default=None)
    parser.add_option("-s", "--size", dest="size",
        help="number of bases to extend each region by", metavar="INT",
        default=DEFAULT_REGION_SIZE)
    (options, args) = parser.parse_args()
    
    if len(args) < 1:
        print __doc__
        sys.exit()
        
    process_pcr_file(args[0], options.output, int(options.size))