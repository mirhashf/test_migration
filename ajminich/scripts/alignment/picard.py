#!/usr/bin/env python
"""Runs the GATK pipeline on an alignment.
- Requires paired reads.

Needed Fixes:
- do not perform realignment if the .intervals file is empty
- do not perform add/remove read groups if the RG tag is already in the header
- do not perform sorting if the input BAM file is already sorted by coordinate
"""
import os, sys, subprocess
import logging
from optparse import OptionParser
from ConfigParser import ConfigParser
from time import time

# Setup logger
FORMAT = '%(asctime)-15s %(levelname)s [%(funcName)s:%(lineno)d] %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger()
logger.setLevel(logging.INFO)
    
class PipelineRunner:

    # Defaults
    DEFAULT_NUM_THREADS         = 1
    DEFAULT_FAST_MODE           = False
    DEFAULT_SINGLE_CHROM        = False
    DEFAULT_OUTPUT_DIR          = ""
    DEFAULT_CONFIG_FILE         = os.path.dirname(__file__) + "/" + "pipeline.cfg"
    
    def __init__(self):

        self.__numThreads__     = self.DEFAULT_NUM_THREADS
        self.__fastMode__       = self.DEFAULT_FAST_MODE
        self.__singleChrom__    = self.DEFAULT_SINGLE_CHROM
        self.__outputDir__      = self.DEFAULT_OUTPUT_DIR
        
        # Import configuration
        self.__loadConfig__(self.DEFAULT_CONFIG_FILE)
        
        # Create timer structure
        self.__runTimes__ = []        
        
        logger.info("Pipeline Runner created.")
        
    def setNumThreads(self, numThreads):
        self.__numThreads__ = numThreads
        
    def setFastMode(self, useFast):
        self.__fastMode__ = useFast
        
    def setChrom(self, chrom):
        self.__chrom__ = chrom
        self.__singleChrom__ = True
        
    def setOutputDir(self, outputDir):
        if (outputDir.endswith("/") or outputDir == ""):
            self.__outputDir__ = outputDir
        else:
            self.__outputDir__ = outputDir + "/"
        
    def setConfigFile(self, configFile):
        self.__loadConfig__(configFile)
    
    def run(self, reference, alignmentFile):
        '''Runs the GATK pipeline, placing the resulting files into the specified output directory.
        '''
        
        logger.info("Running GATK pipeline on alignment '" + alignmentFile + "'.")
        logger.info("Threads:          " + str(self.__numThreads__))
        logger.info("Fast Mode:        " + ("on" if self.__fastMode__ else "off"))
        logger.info("Output Directory: " + self.__outputDir__)
        logger.info("Chromosomes:      " + (("chr" + str(self.__chrom__)) if self.__singleChrom__ else "all"))

        # Get the full file path without the extension
        fileNameBase = self.__outputDir__ + os.path.splitext(os.path.basename(alignmentFile))[0]

        # Create necessary directories        
        if not os.path.exists(self.__outputDir__):
            os.makedirs(self.__outputDir__)        
        
        # Reset timer structure
        self.__runTimes__ = []
        
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
        self.addOrReplaceGroups(alignmentFile, groupsFile, 'coordinate')
        
        # Remove Duplicates
        self.removeDuplicates(groupsFile, markedFile, metricsFile, True)
        self.rebuildIndex(markedFile)

        # Realign Reads        
        if (self.__fastMode__):
            self.fastRealign(reference, markedFile, realignFile, recalFile)
        else:
            self.realign(reference, markedFile, realignFile, intervalsFile)
            self.rebuildIndex(realignFile)
            self.countCovariates(reference, realignFile, recalFile)
            
        self.rebuildIndex(realignFile)
        
        # Table Recalibration
        self.tableRecalibration(reference, realignFile, tableRecalFile, recalFile)
        self.rebuildIndex(tableRecalFile)
        
        # Unified Genotyper Variant Calling
        self.callVariants(reference, tableRecalFile, variantsFile)

        # Output runtimes
        for entry in self.__runTimes__:
            logger.info(entry[0] + ": " + str(entry[1]) + " seconds")
            
        logger.info("GATK pipeline complete: results files available in '" + self.__outputDir__ + "'.")
         
    def convertToBam(self, inFile, outFile):
        logger.info("Converting file '" + inFile + "' to BAM format.")

        startTime = time()        
        
        bamWriter = open(outFile, 'w')
        subprocess.call([self.__samtools__, "view", "-bS", inFile], stdout=bamWriter)
        bamWriter.close()
        self.__checkFile__(outFile)
        
        self.__runTimes__.append(("Convert to BAM", time() - startTime))
    
    def addOrReplaceGroups(self, inFile, outFile, sortOrder):
        logger.info("Adding/Replacing Read Groups in '" + inFile + "'.")
        
        startTime = time()
        
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__picard__ + "/AddOrReplaceReadGroups.jar",
             "I=" + inFile,
             "O=" + outFile,
             "SORT_ORDER=" + sortOrder,
             "RGID=" + self.__RGID__,
             "RGLB=" + self.__RGLB__,
             "RGPL=" + self.__RGPL__,
             "RGSM=" + self.__RGSM__,
             "RGPU=" + self.__RGPU__,
             "VALIDATION_STRINGENCY=" + self.__validation__,
             "TMP_DIR=" + self.__tmp__
            ]
        )
        self.__checkFile__(outFile)
        
        self.__runTimes__.append(("Add or Replace Groups", time() - startTime))
        
    def sort(self, inFile, outFile, sortOrder):
        logger.info("Sorting file '" + inFile + "' by coordinate.")
        
        startTime = time()
        
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__picard__ + "/SortSam.jar",
             "I=" + inFile,
             "O=" + outFile,
             "SORT_ORDER=" + sortOrder,
             "VALIDATION_STRINGENCY=" + self.__validation__,
             "TMP_DIR=" + self.__tmp__
            ]
        )
        self.__checkFile__(outFile)

        self.__runTimes__.append(("Sorting", time() - startTime))
        
    def removeDuplicates(self, inFile, outFile, metricsFile, assumeSorted):
        logger.info("Marking duplicates in '" + inFile + "'.")
        
        startTime = time()
        
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__picard__ + "/MarkDuplicates.jar",
             "I=" + inFile,
             "O=" + outFile,
             "M=" + metricsFile,
             "REMOVE_DUPLICATES=true",
             "ASSUME_SORTED=true",
             "VALIDATION_STRINGENCY=" + self.__validation__,
             "TMP_DIR=" + self.__tmp__
            ]
        )
        self.__checkFile__(outFile)
        
        self.__runTimes__.append(("Remove Duplicates", time() - startTime))
        
    def fastRealign(self, reference, inFile, outFile, recalFile):
        logger.info("Fast Realigning '" + inFile + "'.")
        
        startTime = time()
        
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__fastGatk__ + "/GenomeAnalysisTK.jar",
             "-T", "FastRealign",
             "-R", reference,
             "-I", inFile,
             "-o", outFile,             
             "-standard",
             "-knownSites", self.__dbsnp__,
             "-recalFile", recalFile
            ]
        )
        self.__checkFile__(outFile)
        
        self.__runTimes__.append(("Fast Realign", time() - startTime))
        
    def realign(self, reference, inFile, outFile, intervalsFile):
        logger.info("Realigning '" + inFile + "'.")
        
        startTime = time()
        
        # Realigner Target Creator
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__gatk__ + "/GenomeAnalysisTK.jar",
             "-T", "RealignerTargetCreator",
             "-R", reference,
             "-I", inFile,
             "-o", intervalsFile,
             "-et", self.__et__,
             "-nt", str(self.__numThreads__)
            ]
        )
        self.__checkFile__(intervalsFile)
        
        # Indel Realigner
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__gatk__ + "/GenomeAnalysisTK.jar",
             "-T", "IndelRealigner",
             "-R", reference,
             "-I", inFile,
             "-o", outFile,
             "-targetIntervals", intervalsFile,
             "-et", self.__et__
            ]
        )
        self.__checkFile__(outFile)
        
        self.__runTimes__.append(("Realignment", time() - startTime))
        
    def countCovariates(self, reference, inFile, recalFile):
        logger.info("Counting Covariates in '" + inFile + "'.")
        
        startTime = time()
        
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__gatk__ + "/GenomeAnalysisTK.jar",
             "-T", "CountCovariates",
             "-cov", "ReadGroupCovariate",
             "-cov", "QualityScoreCovariate",
             "-cov", "CycleCovariate",
             "-cov", "DinucCovariate",
             "-R", reference,
             "-I", inFile,
             "-knownSites", self.__dbsnp__,
             "-recalFile", recalFile,
             "-et", self.__et__,
             "-nt", str(self.__numThreads__)
            ]
        )
        self.__checkFile__(recalFile)
        
        self.__runTimes__.append(("Counting Covariates", time() - startTime))
        
    def tableRecalibration(self, reference, inFile, outFile, recalFile):
        logger.info("Recalibrating table in '" + inFile + "'.")
        
        startTime = time()
          
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__gatk__ + "/GenomeAnalysisTK.jar",
             "-T", "TableRecalibration",
             "-R", reference,
             "-I", inFile,
             "-o", outFile,
             "-baq", "RECALCULATE",
             "--doNotWriteOriginalQuals",
             "-recalFile", recalFile,
             "-et", self.__et__
            ]
        )
        self.__checkFile__(outFile)
        
        self.__runTimes__.append(("Table Recalibration", time() - startTime))
        
    def callVariants(self, reference, inFile, outFile):
        logger.info("Performing variant calling on '" + inFile + "'.")
          
        startTime = time()
          
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__gatk__ + "/GenomeAnalysisTK.jar",
             "-T", "UnifiedGenotyper",
             "-R", reference,
             "-I", inFile,
             "-o", outFile,
             "-D", self.__dbsnp__,
             "-A", "AlleleBalance",
             "-A", "DepthOfCoverage",
             "-A", "MappingQualityZero",
             "-baq", "CALCULATE_AS_NECESSARY",
             "-rf", "BadCigar",
             "-dcov", self.__coverageDepth__,
             "-stand_call_conf", self.__standCallConf__,
             "-stand_emit_conf", self.__standEmitConf__,
             "-glm", self.__glm__,
             "-et", self.__et__,
             "-nt", str(self.__numThreads__),
             "-S", self.__validation__
            ]
        )
        self.__checkFile__(outFile)
        
        self.__runTimes__.append(("Variant Calling", time() - startTime))
    
    def rebuildIndex(self, fileName):
        subprocess.call(self.__getJavaArgs__() + [
             "-jar", self.__picard__ + "/BuildBamIndex.jar",
             "I=" + fileName,
             "VALIDATION_STRINGENCY=" + self.__validation__,
             "TMP_DIR=" + self.__tmp__
            ]
        )
        
    def __checkFile__(self, filename):
        try:
            open(filename)
        except IOError as e:
            logger.error("File '" + filename + "' does not exist; aborting.")
            sys.exit()
        
    def __loadConfig__(self, configFile):
        
        logger.info("Loading configuration from '" + configFile + "'.")        
        
        config = ConfigParser()
        config.read(configFile)
        
        # Get executables
        self.__samtools__       = config.get('Executables', 'samtools')
        self.__picard__         = config.get('Executables', 'picard')
        self.__gatk__           = config.get('Executables', 'gatk')
        self.__fastGatk__       = config.get('Executables', 'gatk_fast')
        
        # Get files
        self.__tmp__            = config.get('Files', 'tmp')
        self.__dbsnp__          = config.get('Files', 'dbsnp')
        
        # Get Java settings
        self.__initialHeap__    = config.get('Java', 'initialHeap')
        self.__maximumHeap__    = config.get('Java', 'maximumHeap')
        self.__validation__     = config.get('Java', 'validation')
        
        # Add or Remove Groups
        self.__RGID__           = config.get('AddOrRemoveGroups', 'RGID')
        self.__RGLB__           = config.get('AddOrRemoveGroups', 'RGLB')
        self.__RGPL__           = config.get('AddOrRemoveGroups', 'RGPL')
        self.__RGSM__           = config.get('AddOrRemoveGroups', 'RGSM')
        self.__RGPU__           = config.get('AddOrRemoveGroups', 'RGPU')
        
        # Genome Analysis Toolkit
        self.__et__             = config.get('GenomeAnalysisTK', 'et')
        self.__coverageDepth__  = config.get('GenomeAnalysisTK', 'coverageDepth')
        self.__standCallConf__  = config.get('GenomeAnalysisTK', 'standCallConf')
        self.__standEmitConf__  = config.get('GenomeAnalysisTK', 'standEmitConf')
        self.__glm__            = config.get('GenomeAnalysisTK', 'glm')
        
    def __getJavaArgs__(self):
        return ["java", "-Xms" + self.__initialHeap__, "-Xmx" + self.__maximumHeap__]

 
if __name__ == "__main__":

    usage = "Usage: %prog [options] <reference> <alignment_file> <file2>\n" + \
        "        <reference>           - the genome reference in FASTA format (pre-indexed)\n" + \
        "        <alignment_file>      - the alignment file in SAM or BAM format\n" + \
        "        -h                    - show extended help"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-t", "--threads", dest="threads",
                      help="number of threads to run",
                      type="int",
                      default=PipelineRunner.DEFAULT_NUM_THREADS,
                      metavar="<int>")
    parser.add_option("-f", "--fast", dest="fastMode",
                      help="run FastRealign",
                      action="store_true")
    parser.add_option("-c", "--chrom", dest="chrom",
                      help="specify chromosome to align against",
                      type="int",
                      default=0,
                      metavar="<chr#>")
    parser.add_option("-d", "--dir", dest="outputDir",
                      help="directory to place generated files",
                      type="string",
                      default=PipelineRunner.DEFAULT_OUTPUT_DIR,
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
        
        runner = PipelineRunner()
        
        runner.setNumThreads(options.threads)
        runner.setFastMode(options.fastMode)
        runner.setOutputDir(options.outputDir)
        
        if (options.chrom > 0):
            runner.setChrom(options.chrom)
            
        if (options.configFile != None):
            runner.setConfigFile(options.configFile)

        runner.run(reference, alignFile)        
        
        sys.exit()
    else:
        parser.print_usage(file=sys.stderr)
