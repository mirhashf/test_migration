#!/usr/bin/env python
"""Runs BWA paired-end alignment on a given reads set from a single command.

Notes:
- Requires paired reads.
"""
import os, sys, subprocess
import logging
from optparse import OptionParser

BWA="/home/ajminich/programs/bwa-0.6.1/bwa"
#BWA="/home/ajminich/programs/bwa-0.5.9/bwa"

# Setup logger
FORMAT = '%(asctime)-15s %(levelname)s [%(funcName)s:%(lineno)d] %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
    
def bwa(reference, file1, file2, outFilePrefix, numThreads):
    
    logger.info("Performing first alignment on file '" + file1 + "'.")
    
    saiFile1 = outFilePrefix + "_1.sai"
    saiWriter1 = open(saiFile1, 'w')
    subprocess.call([BWA, "aln", "-t", str(numThreads), reference, file1], stdout=saiWriter1)
    saiWriter1.close()
    
    logger.info("Performing second alignment on file '" + file2 + "'.")
    
    saiFile2 = outFilePrefix + "_2.sai"
    saiWriter2 = open(saiFile2, 'w')
    subprocess.call([BWA, "aln", "-t", str(numThreads), reference, file2], stdout=saiWriter2)
    saiWriter2.close()
        
    logger.info("Running SAM pairing on SAI files '" + saiFile1 + "' and '" + saiFile2 + "'.")        
    
    finalSam = outFilePrefix + ".sam"
    samWriter = open(finalSam, 'w')
    subprocess.call([BWA, "sampe", reference, saiFile1, saiFile2, file1, file2], stdout=samWriter)
    samWriter.close()
        
    logger.info("BWA alignment complete: alignment file available as '" + finalSam + "'.")
 
if __name__ == "__main__":

    usage = "Usage: %prog [options] <reference> <file1> <file2> <output_prefix>"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-t", "--threads", dest="threads",
                      help="The number of threads to run.",
                      type="int",
                      default=1,
                      metavar="<int>")
    
    (options, args) = parser.parse_args()

    if (len(args) == 4):
        bwa(args[0], args[1], args[2], args[3], options.threads)
        sys.exit()
    else:
        parser.print_usage(file=sys.stderr)
