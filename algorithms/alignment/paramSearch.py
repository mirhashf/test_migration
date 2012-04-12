#!/usr/bin/env python
"""%prog [options] <index.fa> <fastq1.fq> <fastq2.fq> <golden_bamfile> <results_file>

Parameter Search v1.0.0

Runs SeqAlto over a range of parameter values, measuring key results
(perfect mappings, far mappings, unmapped reads) for each parameter
combination."""
import os
import sys
import subprocess
from seqalto import SeqAlto
from optparse import OptionParser
import itertools
import ConfigParser
import logging

DEFAULT_CONFIG_FILE = os.path.split(sys.argv[0])[0] + "/paramSearch.cfg"
modes = ["combinatorial", "sequential"]

# Setup logger
FORMAT = '%(asctime)-15s %(levelname)s [%(funcName)s:%(lineno)d] %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def runValidation(executable, alignedFile, goldenFile):
    ''' Runs the AlignStats validation and returns a tuple of (header, report).
    '''
    
    p = subprocess.Popen(["java", "-jar", executable,
        goldenFile, alignedFile,
        "--names", "golden,seqalto",
        "--basicReport",
        "--f1gold"
        ], stdout=subprocess.PIPE)
    
    header = ""
    result = ""
    
    # Poll process for new output until finished
    while True:
        nextline = p.stdout.readline()
        if nextline == '' and p.poll() != None:
            break
        
        # If line starts with "Report:", capture it as the header
        if nextline.startswith("Report:"):
            header = nextline.replace("Report:","")
        elif header != "" and nextline != "":
            result = nextline.replace('\n','')
            
        sys.stdout.write(nextline)
        sys.stdout.flush()
     
    return (header, result)

def runLoop(seqaltoExec, alignstats, index, fastq1, fastq2, goldenFile, outputDir,
            resultsFile, seqOutput, numThreads, paramsListMap):

    # Set up output directory
    if (outputDir is None or outputDir == ""):
        outputDir = "./"
    elif (not outputDir.endswith("/")):
        outputDir = outputDir + "/"
    
    # Create necessary directories        
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)    

    seqalto = SeqAlto(seqaltoExec)

    # Initialize output file
    outWriter = open(outputDir + resultsFile, 'w')
    outWriter.write(getFlags(paramsListMap[0]))
    outWriter.write("Identical,Close,Far,Mapped-Multi,Multi-Mapped,Mapped-Unmapped,Unmapped-Mapped,Time\n")
    
    for paramCombination in paramsListMap:
    
        paramCombination["-p"] = numThreads
        logger.info("Using parameter combination: " + str(paramCombination))
    
        alignedFile = outputDir + seqOutput
    
        # Run SeqAlto
        segfaulted = seqalto.align(index, fastq1, fastq2, alignedFile, paramCombination)
        
        runTime = seqalto.getElapsedTime()
        
        outString = getValues(paramCombination)
        
        if segfaulted:
            result = "segfault"
        else:
            
            # Run validation and get the result
            (header, result) = runValidation(alignstats, alignedFile, goldenFile)
        
        outWriter.write(outString + result + "," + str(runTime) + "\n")
        outWriter.flush()
        
    outWriter.close()
           
def getFlags(paramCombination):
    flags = ""
    for key in paramCombination:
        flags += key + ","
    return flags

def getValues(paramCombination):
    values = ""
    for key in paramCombination:
        values += str(paramCombination[key]) + ","
    return values
           
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
                      default=DEFAULT_CONFIG_FILE,
                      metavar="<cfg_file>")
    parser.add_option("-d", "--outputDir", dest="outputDir",
                      help="Output directory for alignment files",
                      type="string",
                      default=".",
                      metavar="<outputDir>")
    parser.add_option("-m", "--mode", dest="mode",
                      help="Combination mode to use for parameter space",
                      choices=modes,
                      default=modes[0],
                      metavar=modes)
    
    (options, args) = parser.parse_args()

    if (len(args) == 5):
    
        (index, fastq1, fastq2, goldenFile, resultsFile) = args
        
        # Parse configuration file
        logger.info("Using configuration file '" + options.config + "'.")
        config = ConfigParser.SafeConfigParser()
        config.read(options.config)
    
        seqalto         = config.get('Executables', 'seqalto')
        alignstats      = config.get('Executables', 'alignstats')
        
        outputSamFile   = config.get('Configuration', 'outputFile')
        numThreads      = config.getint('Configuration', 'numThreads')
        
        paramsListMap       = parseItems(config.items('Parameters'), options.mode)
        
        runLoop(seqalto, alignstats, index, fastq1, fastq2, goldenFile,
                options.outputDir, resultsFile, outputSamFile, numThreads,
                paramsListMap)
        
        sys.exit()
    else:
        parser.print_help(file=sys.stderr)