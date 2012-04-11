#!/usr/bin/env python
"""%prog [options] <index.fa> <fastq1.fq> <fastq2.fq>

Parameter Search v1.0.0

Runs SeqAlto over a range of parameter values, measuring key results
(perfect mappings, far mappings, unmapped reads) for each parameter
combination."""
import os
import sys
import subprocess
from optparse import OptionParser
import itertools
import ConfigParser
import logging

modes = ["combinatorial", "sequential"]

# Setup logger
FORMAT = '%(asctime)-15s %(levelname)s [%(funcName)s:%(lineno)d] %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def runLoop(seqalto, alignstats, index, fastq1, fastq2, numThreads, kmerSize, paramsListMap):

    for paramCombination in paramsListMap:
        print paramCombination
    
   # aligner -mode align -idx ../chr15_21.fa_"$kmersize".sidx -s -logtostderr -p $numthreads -1 ./ds4_trimmed_1.fq -2 ./ds4_trimmed_2.fq -max_occ $max_num_hash -u > $samfile 2>$samfile.log
           
def parseItems(items, mode):
    
    # Convert values to list of lists
    keys = []
    values = []
    for entry in items:
        keys.append(entry[0])
        
        paramsList = [int(value) for value in entry[1].split(',')]
        values.append(paramsList)
    
    # Perform appropriate parameter combination based on mode
    paramCombos = []
    
    if mode == modes[0]:
        logger.info("Using combinatorial method of searching parameters space.")
        
        # Generate Cartesian product of all selected parameter values
        paramCombos = [element for element in itertools.product(*values)]
            
    elif mode == modes[1]:
        logger.info("Using sequential method of searching parameters space.")
        
        minLen = len(min(values, key=len))
        maxLen = len(max(values, key=len))
        
        if (minLen != maxLen):
            logger.warn("Warning: some parameter lists are as long as " + str(maxLen) + " elements, "
                        "but only creating a total of " + str(minLen) + " params lists.")
            
        for index in range(minLen):
            paramCombos.append([entry[index] for entry in values])
        
    else:
        logger.error("Mode not understood: '" + mode + "'. Exiting.")
        sys.exit(2)
        
    # Place parameters into a list of SeqAlto options
    paramsListMap = []
    
    for paramCombo in paramCombos:
        paramsMap = {}
        for index in range(len(paramCombo)):
            paramsMap[keys[index]] = paramCombo[index]
    
        paramsListMap.append(paramsMap)
        
    return paramsListMap

if __name__ == "__main__":

    parser = OptionParser(usage=__doc__)
    
    # Parse inputs
    parser.add_option("-c", "--config", dest="config",
                      help="Parameter search configuration to use",
                      type="string",
                      default="paramSearch.cfg",
                      metavar="<cfg_file>")
    parser.add_option("-d", "--outputDir", dest="outputDir",
                      help="Output directory for alignment files",
                      type="string",
                      default="./",
                      metavar="<outputDir>")
    parser.add_option("-m", "--mode", dest="mode",
                      help="Combination mode to use for parameter space",
                      choices=modes,
                      default=modes[0],
                      metavar=modes)
    
    (options, args) = parser.parse_args()

    if (len(args) == 3):
    
        (index, fastq1, fastq2) = args
        
        # Parse configuration file
        logger.info("Using configuration file '" + options.config + "'.")
        config = ConfigParser.SafeConfigParser()
        config.read(options.config)
    
        seqalto         = config.get('Executables', 'seqalto')
        alignstats      = config.get('Executables', 'alignstats')
        
        numThreads      = config.getint('Configuration', 'numThreads')
        kmerSize        = config.getint('Configuration', 'kmerSize')
        
        paramsListMap       = parseItems(config.items('Parameters'), options.mode)
        
        runLoop(seqalto, alignstats, index, fastq1, fastq2, numThreads, kmerSize, paramsListMap)
        sys.exit()
    else:
        parser.print_help(file=sys.stderr)