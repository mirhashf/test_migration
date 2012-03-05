#!/usr/bin/env python
"""Runs the GATK pipeline on an alignment.
"""
import os
import sys
import logging
import gatk
from optparse import OptionParser

# Setup logger
FORMAT = '%(asctime)-15s %(levelname)s [%(funcName)s:%(lineno)d] %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
    
def run(reference, alignmentFile, outputDir="", numThreads=1, fastMode=False, chrom=None, configFile=None):
    '''Runs the GATK pipeline, placing the resulting files into the specified output directory.
    '''

    # Set up GATK runner    
    runner = GATK()
        
    runner.setNumThreads(numThreads)
    runner.setFastMode(fastMode)
    runner.setOutputDir(outputDir)
        
    if (chrom > 0):
        runner.setChrom(chrom)
        
    if (configFile != None):
        runner.setConfigFile(configFile)

    logger.info("Running GATK pipeline on alignment '" + alignmentFile + "'.")
    logger.info("Threads:          " + str(numThreads))
    logger.info("Fast Mode:        " + ("on" if fastMode else "off"))
    logger.info("Output Directory: " + outputDir)
    logger.info("Chromosomes:      " + (("chr" + str(chrom)) if (chrom != None) else "all"))

    # Get the full file path without the extension
    fileNameBase = outputDir + os.path.splitext(os.path.basename(alignmentFile))[0]

    # Create necessary directories        
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)        
    
    # Define file names
    groupsFile      = fileNameBase + ".groups.bam"
    markedFile      = fileNameBase + ".marked.bam"
    metricsFile     = fileNameBase + ".metrics"
    intervalsFile   = fileNameBase + ".realign.intervals"
    realignFile     = fileNameBase + ".realign.bam"
    recalFile       = fileNameBase + ".realign.csv"
    tableRecalFile  = fileNameBase + ".recal.bam"
    variantsFile    = fileNameBase + ".vcf"
    
    # Add/Remove Groups
    runner.addOrReplaceGroups(alignmentFile, groupsFile, 'coordinate')
    
    # Remove Duplicates
    runner.removeDuplicates(groupsFile, markedFile, metricsFile, True)
    runner.rebuildIndex(markedFile)

    # Realign Reads        
    if (fastMode):
        runner.fastRealign(reference, markedFile, realignFile, recalFile)
    else:
        runner.realign(reference, markedFile, realignFile, intervalsFile)
        runner.rebuildIndex(realignFile)
        runner.countCovariates(reference, realignFile, recalFile)
        
    runner.rebuildIndex(realignFile)
    
    # Table Recalibration
    runner.tableRecalibration(reference, realignFile, tableRecalFile, recalFile)
    runner.rebuildIndex(tableRecalFile)
    
    # Unified Genotyper Variant Calling
    runner.callVariants(reference, tableRecalFile, variantsFile)

    # Output runtimes
    for entry in runner.getRuntimes():
        logger.info(entry[0] + ": " + str(entry[1]) + " seconds")
        
    logger.info("GATK pipeline complete: results files available in '" + outputDir + "'.")
         
if __name__ == "__main__":

    usage = "Usage: %prog [options] <reference> <alignment_file> <file2>\n" + \
        "        <reference>           - the genome reference in FASTA format (pre-indexed)\n" + \
        "        <alignment_file>      - the alignment file in SAM or BAM format\n" + \
        "        -h                    - show extended help"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-t", "--threads", dest="threads",
                      help="number of threads to run",
                      type="int",
                      default=gatk.GATK.DEFAULT_NUM_THREADS,
                      metavar="<int>")
    parser.add_option("-f", "--fast", dest="fastMode",
                      help="run FastRealign",
                      action="store_true")
    parser.add_option("-c", "--chrom", dest="chrom",
                      help="specify chromosome to align against",
                      type="int",
                      default=None,
                      metavar="<chr#>")
    parser.add_option("-d", "--dir", dest="outputDir",
                      help="directory to place generated files",
                      type="string",
                      default=gatk.GATK.DEFAULT_OUTPUT_DIR,
                      metavar="<outputDir>")
    parser.add_option("--config", dest="configFile",
                      help="configuration file",
                      type="string",
                      default=None,
                      metavar="<configFile>")
    
    (options, args) = parser.parse_args()

    if (len(args) == 2):

        reference = args[0]
        alignFile = args[1]
        
        run(reference,
            alignFile,
            outputDir=options.outputDir,
            numThreads=options.numThreads,
            fastMode=options.fastMode,
            chrom=options.chrom,
            configFile=options.configFile)
        
        sys.exit()
    else:
        parser.print_usage(file=sys.stderr)
