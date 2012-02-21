#!/usr/bin/python

'''
    Project Kakashi
    run_test.py
    
    Runs a complete SeqAlto test by performing an alignment and performing
    several analyses.
    
    Written by AJ Minich, January 2012    
'''

import argparse
from alignment_runner import alignment_runner
from analysis_runner import analysis_runner

# Default Parameters
DATA_DIRECTORY_DEFAULT = "/mnt/scratch0/ajminich"
INDEX_DEFAULT = "chr21.fa_22.midx"
READS1_DEFAULT = "reads_1.fq"
READS2_DEFAULT = "reads_2.fq"
MAP_FILES_DEFAULT = ["randCEUma.fa.map", "randCEUpa.fa.map"]
MAP_NAMES_DEFAULT = ["ma", "pa"]

if __name__ == '__main__':
    
    index_file = '/'.join([DATA_DIRECTORY_DEFAULT, INDEX_DEFAULT])
    reads1_file = '/'.join([DATA_DIRECTORY_DEFAULT, READS1_DEFAULT])
    reads2_file = '/'.join([DATA_DIRECTORY_DEFAULT, READS2_DEFAULT])
    map_files = ['/'.join([DATA_DIRECTORY_DEFAULT, map_file]) for map_file in MAP_FILES_DEFAULT]
    
    # Run the alignment
    aligner = alignment_runner()
    
    #align_output = aligner.run_alignment(index_file, reads1_file, reads2_file)
    align_output = "results.sam"
    
    # Get the number of reads
    num_reads = aligner.get_num_reads(reads1_file)
    
    # Run the analysis
    analyzer = analysis_runner()
    
    analyzer.run_analyzer(map_files, MAP_NAMES_DEFAULT, align_output, num_reads)
    
    print "Test complete."
    
    
    