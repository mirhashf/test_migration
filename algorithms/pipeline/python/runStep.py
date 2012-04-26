#!/usr/bin/env python
"""Runs a single step in the variant calling pipeline.
"""
import os
import sys
import logging
import json
from gatk import GATK
from picard import Picard
from optparse import OptionParser
from time import time

# Setup logger
FORMAT = '%(asctime)-15s %(levelname)s [%(funcName)s:%(lineno)d] %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
    
if __name__ == "__main__":

    usage = "Usage: %prog [options] <reference> <alignment_file> <config_json>\n" + \
        "        <reference>           - the genome reference in FASTA format (pre-indexed)\n" + \
        "        <alignment_file>      - the alignment file in SAM or BAM format\n" + \
        "        <config_json>         - JSON file containing pipeline configuration\n" + \
        "        -h                    - show extended help"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-d", "--dir", dest="outputDir",
                      help="directory to place generated files",
                      type="string",
                      default="./",
                      metavar="<outputDir>")
    parser.add_option("-c", "--chrom", dest="chrom",
                      help="specify chromosomal region to align against",
                      type="string",
                      default=None,
                      metavar="<region>")
    parser.add_option("--recalFile", dest="recalFile",
                      help="recalibration file (if Table Recalibration or "
                      "Combined Genotyper is selected mode)",
                      type="string",
                      default=None,
                      metavar="<recalFile>")
    
    # Running options
    parser.add_option("--groups", dest="groups",
                      help="Sort and add/replace read groups",
                      action="store_true")
    parser.add_option("--removeDuplicates", dest="removeDups",
                      help="Remove duplicates",
                      action="store_true")
    parser.add_option("--realign", dest="realign",
                      help="Standard realign",
                      action="store_true")
    parser.add_option("--fastRealign", dest="fastRealign",
                      help="Fast realign",
                      action="store_true")
    parser.add_option("--countCovars", dest="countCovars",
                      help="Count covariates",
                      action="store_true")
    parser.add_option("--recal", dest="tableRecalibration",
                      help="Table recalibration",
                      action="store_true")
    parser.add_option("--callVariants", dest="callVariants",
                      help="Call variants using Unified Genotyper",
                      action="store_true")
    parser.add_option("--combinedVariants", dest="combinedVariants",
                      help="Recalibrate and call variants using Combined Genotyper",
                      action="store_true")
    
    (options, args) = parser.parse_args()

    if (len(args) < 3):
        parser.print_usage(file=sys.stderr)
        sys.exit()
                
    reference   = args[0]
    inFile   = args[1]
    configFile  = args[2]
        
    # Get configuration
    configData = open(configFile).read()
    config = json.loads(configData)
    
    name = config['name']
    
    # Set up runners
    picardRunner = Picard(configFile)
    gatkRunner = GATK(configFile)
    
    # Set up output directory
    outputDir = options.outputDir
    if (outputDir == ""):
        outputDir = "./"
    elif (not outputDir.endswith("/")):
        outputDir = outputDir + "/"
    
    # Create necessary directories        
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)    
    
    # Get the full file path without the extension
    fileNameBase = outputDir + os.path.splitext(os.path.basename(inFile))[0]
    
    if (options.chrom > 0):
        gatkRunner.setChrom(options.chrom)
        
    logger.info("Configuration:      " + name)
    logger.info("Output Directory:   " + outputDir)
    logger.info("Chromosome Filter:  " +
        (("chr" + str(options.chrom)) if (options.chrom != None) else "none"))
    
    # Define file names
    groupsFile      = fileNameBase + ".groups.bam"
    markedFile      = fileNameBase + ".marked.bam"
    metricsFile     = fileNameBase + ".metrics"
    intervalsFile   = fileNameBase + ".realign.intervals"
    realignFile     = fileNameBase + ".realign.bam"
    recalFile       = fileNameBase + ".realign.csv"
    tableRecalFile  = fileNameBase + ".recal.bam"
    variantsFile    = fileNameBase + ".vcf"
    
    startTime = time()
    
    if options.groups:
        picardRunner.addOrReplaceGroups(inFile, groupsFile)
        logger.info("Add or Replace Groups completed in " + str(time() - startTime) + "seconds.")
        
    elif options.removeDups:
        picardRunner.removeDuplicates(inFile, markedFile, metricsFile, True)
        picardRunner.rebuildIndex(markedFile)
        logger.info("Duplicate Removal completed in " + str(time() - startTime) + "seconds.")
        
    elif options.realign:
        gatkRunner.realign(reference, inFile, realignFile, intervalsFile)
        picardRunner.rebuildIndex(realignFile)
        logger.info("Realignment completed in " + str(time() - startTime) + "seconds.")
        
    elif options.countCovars:
        gatkRunner.countCovariates(reference, inFile, recalFile)
        logger.info("Realignment completed in " + str(time() - startTime) + "seconds.")
        
    elif options.fastRealign:
        gatkRunner.fastRealign(reference, inFile, realignFile, recalFile)
        picardRunner.rebuildIndex(realignFile)
        logger.info("Fast Realignment completed in " + str(time() - startTime) + "seconds.")
        
    elif options.tableRecalibration:
        
        if (options.recalFile is None):
            logger.error("Please specify a recalibration file using the --recalFile option.")
        else:
            gatkRunner.tableRecalibration(reference, inFile, tableRecalFile, options.recalFile)
            picardRunner.rebuildIndex(tableRecalFile)
            logger.info("Table Recalibration completed in " + str(time() - startTime) + "seconds.")
    
    elif options.callVariants:
        gatkRunner.callVariants(reference, inFile, variantsFile)
        logger.info("Variant Calling completed in " + str(time() - startTime) + "seconds.")
        
    elif options.combinedVariants:
        
        if (options.recalFile is None):
            logger.error("Please specify a recalibration file using the --recalFile option.")
        else:
            gatkRunner.combinedGenotyper(reference, inFile, variantsFile, recalFile)
            logger.info("Combined Variant Calling completed in " + str(time() - startTime) + "seconds.")
    
    else:
        logger.warn("No operation selected! Please select from one of our fine pipeline steps.")
        