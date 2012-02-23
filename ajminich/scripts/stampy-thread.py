#!/usr/bin/env python
"""Runs Stampy with multiple threads to speed up the alignment.

Usage:
    -1 <file1>, --file1 <file1>             The first FASTQ input.
    -2 <file2>, --file2 <file2>             The second FASTQ input.
    -o <outfile>, --output <outfile>        The prefix to use for final BAM output.
    -t <threads>, --threads <threads>       The number of threads to run.
    
Notes:
- Using the -M tag is unnecessary, since the multi-threader takes care of
  providing the correct reads files to Stampy.
"""
import os
import sys
import subprocess
from threading import Thread
from math import floor

STAMPY="/home/ajminich/programs/stampy"
MERGER="/home/ajminich/programs/picard/dist/MergeSamFiles.jar"
INPUT1_FLAGS=["-1", "--file1"]
INPUT2_FLAGS=["-2", "--file2"]
OUTPUT_FLAGS=["-o", "--output"]
THREAD_FLAGS=["-t", "--threads"]

class StampyThread(Thread):

    def __init__(self, readsFile1, readsFile2, outFilePrefix, args):
        Thread.__init__(self)
        self.readsFile1 = readsFile1
        self.readsFile2 = readsFile2
        self.outFilePrefix = outFilePrefix
        self.args = args

    def run(self):
        
        # Create read args
        reads = ["-M", self.readsFile1 + "," + self.readsFile2]
        
        # Perform the alignment
        samFile = self.outFilePrefix + ".sam"
        
        print >>sys.stderr, "Aligning partition files '%s' and '%s' into '%s'." % \
            (self.readsFile1, self.readsFile2, samFile)

        samWriter = open(samFile, 'w')
        subprocess.call([STAMPY] + self.args + reads, stdout=samWriter)
        samWriter.close()

        # Convert from SAM to BAM
        bamFile = self.outFilePrefix + ".bam"
        bamWriter = open(bamFile, 'w')
        subprocess.call(["samtools", "view", "-bS", samFile], stdout=bamWriter)
        bamWriter.close()
        

def stampy_threaded(file1, file2, outFilePrefix, numThreads, args):
    
    # Divide files up into equal parts
    
    print >>sys.stderr, "Getting size of each input file."
    lines1 = file_len(file1)
    lines2 = file_len(file2)
    
    if (lines1 != lines2):
        print >>sys.stderr, "Error: file 1 (%i lines) is not the same length as file 2 (%i lines)." % \
            (lines1, lines2)
        return
    
    print >>sys.stderr, "Dividing %i-line files into %i equal parts." % (lines1, numThreads)
    (divFiles1, divFiles2) = partitionFiles(file1, file2, lines1, numThreads)

    # Start multi-threaded Stampy
    threads = []
    outFiles = []
    for threadIndex in range(numThreads):

        reads1 = divFiles1[threadIndex]
        reads2 = divFiles2[threadIndex]

        threadOutFile = outFilePrefix + "_thread_" + str(threadIndex)
        outFiles.append(threadOutFile)
        
        # Activate thread here
        thread = StampyThread(reads1, reads2, threadOutFile, args)
        thread.start()
        threads.append(thread)
    
    print >>sys.stderr, "All threads have been spooled, now waiting for Stampy to complete."

    # Wait for threads to finish
    for threadIndex in range(numThreads):
        threads[threadIndex].join()
        print >>sys.stderr, "Stampy Thread #%i: finished" % (threadIndex + 1)
        
    print >>sys.stderr, "All threads finished."
        
    # Concatenate BAM files
    finalFile = outFilePrefix + ".bam"

    bamFiles = [file + ".bam" for file in outFiles]
    
    subprocess.call(["samtools", "merge", outFilePrefix + ".bam"] + bamFiles)
        
    # Cleanup
    print >>sys.stderr, "Cleaning up partitioned files."
    for fileIndex in range(numThreads):
        os.remove(divFiles1[fileIndex])
        os.remove(divFiles2[fileIndex])
        os.remove(outFiles[fileIndex] + ".sam")
        os.remove(outFiles[fileIndex] + ".bam")
    
    print >>sys.stderr, "Stampy alignment complete: alignment file available as "
    

def partitionFiles(file1, file2, numLines, numParts):
    
    # To determine the number of lines in the divided files, we will divide
    # the total number of lines by the number of parts we want. But we want
    # to maintain the FASTQ read data, so we will make sure the number is a
    # multiple of 4.
    fileLengths = 4 * floor(numLines / (4 * numParts))
    
    print >>sys.stderr, "Partitioned files will be %i lines each." % fileLengths
    
    f1 = open(file1, 'r')
    f2 = open(file2, 'r')

    dividedFiles1 = []
    dividedFiles2 = []
    
    for fileIndex in range(numParts):
        fout1 = file1 + "_" + str(fileIndex)
        fout2 = file2 + "_" + str(fileIndex)
        
        dividedFiles1.append(fout1)
        dividedFiles2.append(fout2)

        print >>sys.stderr, "Creating partitioned files '%s' and '%s'." % (fout1, fout2)

        fw1 = open(fout1, 'w')
        fw2 = open(fout2, 'w')
        
        for lineNum in range(fileLengths):
            fw1.write(f1.readline())
            fw2.write(f2.readline())
            
        fw1.close()
        fw2.close()
    
    # Finish up the files if necessary, just putting all the extra reads in the last file
    line1 = f1.readline()
    line2 = f2.readline()
        
    fw1 = open(dividedFiles1[-1], 'a')
    fw2 = open(dividedFiles2[-1], 'a')
        
    while (line1 != "" or line2 != ""):
        fw1.write(f1.readline())
        fw2.write(f2.readline())
        
    fw1.close()
    fw2.close()
        
    f1.close()
    f2.close()
    
    return (dividedFiles1, dividedFiles2)
    
def file_len(fname):
    i = -1  # in case file is empty
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
        
    print >>sys.stderr, "File '%s': %i lines." % (fname, i + 1)
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
    (outfile_def, outFilePrefix, args) = parse(args, OUTPUT_FLAGS)    
    (threads_def, threads, args) = parse(args, THREAD_FLAGS)

    if (file1_def and file2_def and outfile_def and threads_def):
        stampy_threaded(file1, file2, outFilePrefix, int(threads), args)
        sys.exit()
    else:
        print >>sys.stderr, __doc__
        subprocess.call([STAMPY])
        
