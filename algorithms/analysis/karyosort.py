#!/usr/bin/python
'''Reorders entries in a file in karyotypic order.

Usage: karyosort.py <file> <fasta_index> <outfile>

Note that this script will temporarily generate several
files in the current folder. These files will be deleted
after the sorting is complete.
'''

import os, sys

def get_contig(text):
    return text.replace('\n', '').split('\t')[0]
    
def get_contig_file(text):
    return text + ".edits"

def karyosort(input_file, output_file, fasta_index):
    
    contigs = [ ]
    
    # Get the FASTA index order
    for line in open(fasta_index, 'r'):
        contigs.append(get_contig(line))
        
    files_open = { }
    
    print "Splitting edits into separate contig files."    
    for line in open(input_file, 'r'):
        contig = get_contig(line)
        contig_file = get_contig_file(contig)
        
        if contig not in contigs:
            print "Warning: contig '" + contig + "' not recognized."
            continue
        
        if contig not in files_open:
            current_file = open(contig_file, 'w')
            files_open[contig] = current_file
            
        files_open[contig].write(line)
        
    output_writer = open(output_file, 'w')
        
    for contig in contigs:
        
        if contig not in files_open:
            print "Warning: reference contig '" + contig + "' not found in edits file."
            continue
        
        current_file = files_open[contig]
        current_file.flush()
        current_file.close()
        current_file = open(get_contig_file(contig), 'r')
        
        for line in current_file.readlines():
            output_writer.write(line)
            
        # Close and delete the file
        current_file.close()
        os.remove(get_contig_file(contig))
        
    # Close the final output file
    output_writer.flush()
    output_writer.close()
    
    print "Sorted into karyotypic order."


if __name__ == '__main__':
    
    if len(sys.argv) < 4:
        print __doc__
        exit(1)
        
    input_file = sys.argv[1]
    fasta_index = sys.argv[2]
    output_file = sys.argv[3]
    
    karyosort(input_file, output_file, fasta_index)
    
    
    