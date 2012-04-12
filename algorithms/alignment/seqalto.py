"""A Python wrapper for SeqAlto."""
import os
import sys
import subprocess
import logging
from time import time

# Setup logger
FORMAT = '%(asctime)-15s %(levelname)s [%(funcName)s:%(lineno)d] %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

class SeqAlto:
    
    def __init__(self, executable):
        
        self._exec = executable
        
    def align(self, index, fastq1, fastq2, seqOutput, paramsMap):
        ''' Runs SeqAlto with the specified parameters.
        Returns true if SeqAlto encounters a segmentation fault.
        '''

        # SeqAlto setup
        seqaltoArgs = [self._exec, "-mode", "align",
            "-idx", index,
            "-1", fastq1,
            "-2", fastq2,
            "-logtostderr"]
    
        optionsList = []
        for option in paramsMap:
            optionsList.extend([option, str(paramsMap.get(option))])
            
        segfaulted=False
            
        startTime = time()
            
        # Run SeqAlto
        samWriter = open(seqOutput, 'w')
        p = subprocess.Popen(seqaltoArgs + optionsList, stdout=samWriter, stderr=subprocess.PIPE)
        
        # Poll process for new output until finished
        while True:
            nextline = p.stderr.readline()
            if nextline == '' and p.poll() != None:
                break
            
            if "segmentation fault" in nextline.lower():
                logger.warn("SeqAlto experienced a segmentation fault on the specified paramters.")
                segfaulted = True
                
            sys.stderr.write(nextline)
            sys.stderr.flush()
         
        exitCode = p.returncode
    
        if (exitCode != 0):
            logger.error("SeqAlto returned an error code of 0.")
            segfaulted = True
        
        # Wait for thread to finish
        p.wait()
        
        self._elapsedTime = time() - startTime
        
        samWriter.close()
        
        return segfaulted
    
    def getElapsedTime(self):
        return self._elapsedTime