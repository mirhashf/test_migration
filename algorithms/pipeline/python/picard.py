#!/usr/bin/env python
"""Wrapper for the Picard postprocessing tools.

Issues:
- performs add/replace read groups even if the RG tag is already in the header.
"""
import subprocess
from runner import Runner

class Picard(Runner):

    def __init__(self, configFile):
        Runner.__init__(self)
        
        self.loadConfig(configFile)
        
        self.logger.info("Picard Runner created.")
        
    def addOrReplaceGroups(self, inFile, outFile):
        ''' Runs Add or Replace Groups from the Picard library.
        
        Determines if the file has already been sorted by coordinate. If it hasn't, the sorting will
        be performed as part of the group adding and replacement.
        '''
        
        self.logger.info("Adding/Replacing Read Groups in '" + inFile + "'.")
        
        sortOrder = "coordinate"
        if self.getSorted(inFile, sortOrder):
            self.logger.info("File '" + inFile + "' is already sorted by coordinate; no need to re-sort.")
            sortOrder = "unsorted"
        else:
            self.logger.info("File '" + inFile + "' is not already sorted; will sort during read group replacement.")
        
        requiredArgs = [
             "I=" + inFile,
             "O=" + outFile,
             "SORT_ORDER=" + sortOrder,
             "VALIDATION_STRINGENCY=" + self._validation,
             "TMP_DIR=" + self._tmp
            ]
        
        self.run(self._picard + "/AddOrReplaceReadGroups.jar", requiredArgs, self.getArgs('addReplaceGroups'))
        self.checkFile(outFile)
        
    def sort(self, inFile, outFile, sortOrder):
        ''' Runs SAM Sorter from the Picard library.
        
        Does not determine if the file has already been sorted in the requested order.
        '''
        
        self.logger.info("Sorting file '" + inFile + "' by coordinate.")
        
        requiredArgs = [
             "I=" + inFile,
             "O=" + outFile,
             "SORT_ORDER=" + sortOrder,
             "VALIDATION_STRINGENCY=" + self._validation,
             "TMP_DIR=" + self._tmp
            ]

        self.run(self._picard + "/SortSam.jar", requiredArgs, self.getArgs('sort'))
        self.checkFile(outFile)
        
    def removeDuplicates(self, inFile, outFile, metricsFile, assumeSorted=True):
        ''' Runs Duplicate Marking and Removal from the Picard library.
        
        Outputs the duplicate removal metrics to a file with name specified by 'metricsFile'.
        
        Assumes that the file has already been sorted by coordinate. If the assumeSorted flag is
        set to False, will first sort by coordinate.
        '''
        
        self.logger.info("Marking duplicates in '" + inFile + "'.")
        
        requiredArgs = [
             "I=" + inFile,
             "O=" + outFile,
             "M=" + metricsFile,
             "REMOVE_DUPLICATES=true",
             "ASSUME_SORTED=" + ("true" if assumeSorted else "false"),
             "VALIDATION_STRINGENCY=" + self._validation,
             "TMP_DIR=" + self._tmp
            ]
        
        self.run(self._picard + "/MarkDuplicates.jar", requiredArgs, self.getArgs('removeDuplicates'))
        self.checkFile(outFile)
        
    def getSorted(self, inFile, order):
        ''' Determines whether the specified file is sorted in the specified order.
        '''
        return ("SO:" + order) in self.getHeader(inFile)
        
    def getHeader(self, inFile):
        ''' Gets the file's SAM-format header.
        
        Assumes that the input file is in BAM format.
        '''
        
        p = subprocess.Popen([self._samtools, "view", "-H", inFile], stdout=subprocess.PIPE)
        (header, err) = p.communicate()
        return header
        
    def rebuildIndex(self, fileName):
        ''' Runs Build Bam Index from the Picard library.
        
        Assuming the file's name is <name>.bam, builds a BAM index and saves it to <name>.bai.
        '''
        
        requiredArgs = [
             "I=" + fileName,
             "VALIDATION_STRINGENCY=" + self._validation,
             "TMP_DIR=" + self._tmp
            ]
        
        self.run(self._picard + "/BuildBamIndex.jar", requiredArgs, self.getArgs('rebuildIndex'))
        
    def loadConfig(self, configFile):
        
        config = Runner.loadConfig(self, configFile)
        
        # Get executables
        self._picard = self._config['picard']['executable']
        self._samtools = self._config['picard']['samtools']
    
    def getArgs(self, tool):
        ''' Extracts the specified tool's arguments and options from the config and returns it as
        a subprocess-executable list.
        
        Returns an empty list if the configuration file has no arguments or options for the specified tool.
        '''
        
        if 'picard' in self._config:
            if tool in self._config['picard']:
                entry = self._config['picard'][tool]
                args = []
            
                if 'arguments' in entry:
                    args.extend(entry['arguments'])
                
                if 'options' in entry:
                    args.extend([opt + "=" + str(entry['options'][opt]) for opt in entry['options']])
                
                return args
        
        # Entry not available; return empty array
        return []
        
