#!/usr/bin/env python
"""Implements a base Runner class containing Runner utilities.

This class should never be run itself - it is an abstract class containing functionality
for other wrappers.
"""
import os
import sys
import logging
import json
import subprocess

class Runner:
    
    def __init__(self):
        
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)

    def loadConfig(self, configFile):
        ''' Loads a configuration from the specified file.
            Note that all of the following configuration values must be present
            in the configuration file; see the documentation for details.
        '''
        
        configData = open(configFile).read()
        self._config = json.loads(configData)
    
        # Get files
        self._tmp = self._config['files']['tmp']
        
        # Get logger format
        self._logformat = self._config['loggerFormat']
        logging.basicConfig(format=self._logformat)
        
        # Get Java settings
        self._initialHeap = self._config['java']['initialHeap']
        self._maximumHeap = self._config['java']['maximumHeap']
        self._validation = self._config['java']['validation']
    
    def checkFile(self, filename):
        ''' Checks that the specified file exists. If the file does not exist, the function
            will throw an IO error.
        '''
        
        try:
            open(filename)
        except IOError as e:
            raise RuntimeError('File ' + filename + ' does not exist; aborting.')
    
    def getVal(self, config, valueName, defaultValue):
        ''' Gets a specified value from the provided configuration dictionary.
            If the value is not in the dictionary, returns the default value.
        '''
        
        return config[valueName] if valueName in config else defaultValue

    def run(self, jarFile, requiredArgs, configurationArgs):
        ''' Builds the Java command and runs it in a subprocess.
        
            jarFile           - path of executable JAR file to run
            requiredArgs      - required arguments (mode to run, input files, output files, etc.)
            configurationArgs - optional arguments coming from the configuration parameters
            
            Note: there is no real difference between the way required and optional arguments are treated,
            but it's useful to separate the lists according to the source of the arguments.
        '''
        
        args = []
        
        # Run Java with appropriate initial heap and maximum heap sizes
        args.extend(["java", "-Xms" + self._initialHeap, "-Xmx" + self._maximumHeap])
        
        # Specify the JAR file
        args.extend(["-jar", jarFile])
        
        # Include the user-specified arguments
        args.extend(requiredArgs)
        args.extend(configurationArgs)
        
        # Run the command
        subprocess.call(args)
        
        
        
        
