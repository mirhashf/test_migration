#!/usr/bin/env python
"""Runs the variant calling pipeline on an alignment.
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

trueVals = ['on', 'enable', 'true', 1]
    
def run(reference, alignmentFile, jsonConfigFile, outputDir="", chrom=None):
    '''Runs the variant calling pipeline, placing the resulting files into the specified output directory.
    '''

    # Create timer structure
    runTimes = []        
        
    # Get configuration
    configData = open(jsonConfigFile).read()
    config = json.loads(configData)
    
    name                = config['name']
    fastRealign         = config['gatk']['fastRealignMode'] in trueVals
    combinedGenotyper   = config['gatk']['combinedMode'] in trueVals
    runAddReplaceGroups = config['picard']['runAddReplaceGroups'] in trueVals
    runRemoveDuplicates = config['picard']['runRemoveDuplicates'] in trueVals
    runRealign          = config['gatk']['runRealign'] in trueVals
    runRecal            = config['gatk']['runRecal'] in trueVals
    
    # Set up runners
    picardRunner = Picard(jsonConfigFile)
    gatkRunner = GATK(jsonConfigFile)
    
    # Set up output directory
    if (outputDir is None or outputDir == ""):
        outputDir = "./"
    elif (not outputDir.endswith("/")):
        outputDir = outputDir + "/"
    
    # Create necessary directories        
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)    
    
    if (chrom > 0):
        gatkRunner.setChrom(chrom)
        
    logger.info("Running variant calling pipeline on alignment '" + alignmentFile + "'.")
    logger.info("Configuration:      " + name)
    logger.info("Output Directory:   " + outputDir)
    logger.info("Chromosome Filter:  " + (chrom if (chrom != None) else "none"))
    logger.info("Add/Replace Groups: " + ("enabled" if runAddReplaceGroups else "disabled"))
    logger.info("Remove Duplicates:  " + ("enabled" if runRemoveDuplicates else "disabled"))
    logger.info("Realignment:        " + ("enabled" if runRealign else "disabled"))
    logger.info("Fast Realigner:     " + ("on" if fastRealign else "off"))
    logger.info("Recalibration:      " + ("enabled" if runRecal else "disabled"))
    logger.info("Combined Genotyper: " + ("on" if combinedGenotyper else "off"))
    
    # Get the full file path without the extension
    fileNameBase = outputDir + os.path.splitext(os.path.basename(alignmentFile))[0]

    # Define file names
    groupsFile      = fileNameBase + ".groups.bam"
    markedFile      = fileNameBase + ".marked.bam"
    metricsFile     = fileNameBase + ".metrics"
    intervalsFile   = fileNameBase + ".realign.intervals"
    realignFile     = fileNameBase + ".realign.bam"
    recalFile       = fileNameBase + ".realign.csv"
    tableRecalFile  = fileNameBase + ".recal.bam"
    variantsFile    = fileNameBase + ".vcf"
    
    pipelineStartTime = time()
    
    try:

        # Add or Replace Read Groups
        if runAddReplaceGroups:
            startTime = time()
            picardRunner.addOrReplaceGroups(alignmentFile, groupsFile)
            runTimes.append(("Add or Replace Groups", time() - startTime))
        else:
            groupsFile = alignmentFile
            
        # Remove Duplicates
        if runRemoveDuplicates:
            startTime = time()
            picardRunner.removeDuplicates(groupsFile, markedFile, metricsFile, True)
            picardRunner.rebuildIndex(markedFile)
            runTimes.append(("Remove Duplicates", time() - startTime))
        else:
            markedFile = groupsFile
            
        if runRealign:
            
            # Realign Reads        
            if (fastRealign):
                startTime = time()
                gatkRunner.fastRealign(reference, markedFile, realignFile, recalFile)
                picardRunner.rebuildIndex(realignFile)
                runTimes.append(("Fast Realign", time() - startTime))
            else:
                startTime = time()
                gatkRunner.realign(reference, markedFile, realignFile, intervalsFile)
                picardRunner.rebuildIndex(realignFile)
                runTimes.append(("Realignment", time() - startTime))
                
                startTime = time()
                gatkRunner.countCovariates(reference, realignFile, recalFile)
                runTimes.append(("Counting Covariates", time() - startTime))
                
        else:
            realignFile = markedFile
            
            if combinedGenotyper or runRecal:
                # Must count covariates and create recal file
                startTime = time()
                gatkRunner.countCovariates(reference, realignFile, recalFile)
                runTimes.append(("Counting Covariates", time() - startTime))        
        
        # Genotyping
        if (combinedGenotyper):
            startTime = time()
            gatkRunner.combinedGenotyper(reference, realignFile, variantsFile, recalFile)
            runTimes.append(("Combined Genotyping", time() - startTime))
        else:
            
            if runRecal:
                # Table Recalibration
                startTime = time()
                gatkRunner.tableRecalibration(reference, realignFile, tableRecalFile, recalFile)
                picardRunner.rebuildIndex(tableRecalFile)
                runTimes.append(("Table Recalibration", time() - startTime))
            else:
                tableRecalFile = realignFile
            
            # Unified Genotyper Variant Calling
            startTime = time()
            gatkRunner.callVariants(reference, tableRecalFile, variantsFile)
            runTimes.append(("Variant Calling", time() - startTime))
        
    except RuntimeError as error:
        logger.error(error.args)

    # Output runtimes
    for entry in runTimes:
        logger.info(entry[0] + ": " + str(entry[1]) + " seconds")
        
    logger.info("Variant calling pipeline completed in " + str(time() - pipelineStartTime) + " seconds.");
    logger.info("Results files available in '" + outputDir + "'.")
         
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
                      default="",
                      metavar="<outputDir>")
    parser.add_option("-c", "--chrom", dest="chrom",
                      help="specify chromosomal region to align against",
                      type="string",
                      default=None,
                      metavar="<region>")
    
    (options, args) = parser.parse_args()

    if (len(args) == 3):

        reference   = args[0]
        alignFile   = args[1]
        configFile  = args[2]
        
        run(reference,
            alignFile,
            outputDir=options.outputDir,
            chrom=options.chrom,
            jsonConfigFile=configFile
            )
        
    else:
        parser.print_usage(file=sys.stderr)
