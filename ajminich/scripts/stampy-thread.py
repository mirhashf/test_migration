#!/usr/bin/env python
"""Runs Stampy with multiple threads to speed up the alignment.

Usage:
    -1 <file1>, --file1 <file1>             The first FASTQ input.
    -2 <file2>, --file2 <file2>             The second FASTQ input.
    -t <threads>, --threads <threads>       The number of threads to run.
    
Notes:
- Using the -M tag is unnecessary, since the multi-threader takes care of
  providing the correct reads files to Stampy.
"""
import os
import sys
import subprocess
from threading import Thread
from optparse import OptionParser

STAMPY="/home/ajminich/programs/stampy"
INPUT1_FLAGS=["-1", "--file1"]
INPUT2_FLAGS=["-2", "--file2"]
THREAD_FLAGS=["-t", "--threads"]

def stampy_threaded(file1, file2, numThreads, args):
    
    print "Running Stampy on '%s' and '%s' with %i threads." % (file1, file2, numThreads)
    
    # Divide files up into equal parts
    lines1 = file_len(file1)
    lines2 = file_len(file2)
    
    if (lines1 != lines2):
        print "Error: file 1 (%i lines) is not the same length as file 2 (%i lines)." % \
            (lines1, lines2)
    else:
        print "Dividing %i-line files into %i equal parts." % (lines1, numThreads)
    
    for threadIndex in range(numThreads):
        
        print "Starting Stampy thread #%i." % (threadIndex + 1)

        #threadArgs = new list
        #subprocess.POpen(STAMPY, len(args), *args)

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def parse(args, flags):
    
    for flag in flags:
        if flag in args:
            idx = args.index(flag)
            value = args[idx+1]
            
            # Clean up args
            args.remove(flag)         
            args.remove(value)
            
            return (True, value, args)
    
    return (False, "-1", args)

if __name__ == "__main__":

    args = sys.argv
    del args[0]
    
    # Get input files
    (file1_def, file1, args) = parse(args, INPUT1_FLAGS)    
    (file2_def, file2, args) = parse(args, INPUT2_FLAGS)    
    (threads_def, threads, args) = parse(args, THREAD_FLAGS)

    if (file1_def and file2_def and threads_def):
        stampy_threaded(file1, file2, int(threads), args)
        sys.exit()
    else:
        print __doc__
        subprocess.call([STAMPY])
        
