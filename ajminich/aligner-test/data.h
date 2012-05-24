/**
 * data.h
 * Declares data generators for alignment unit testing.
 *
 * Written by AJ Minich, Jan 2012
 */

#ifndef DATA_H

#define DATA_H

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdint.h>
#include "../src/fast.h"

/**
 * The total number of reads maintained in the data structure, including
 * mismatch reads.
 */
const int NUM_READS = 4;

/**
 * The number of mismatch reads available.
 */
const int NUM_MISMATCH_READS = 1;

using namespace std;

struct Read_t {
    string name;
    
    string forwardSequence;
    string forwardQuality;
    int forwardLength;
    
    string reverseSequence;
    string reverseQuality;
    int reverseLength;
    
    // Mismatch parameters
    bool isMismatched;
    vector<uint32_t> forwardStartLocs, reverseStartLocs;
    vector<int> forwardDiffs, reverseDiffs;
    
};

class testDataGenerator {

  public:

    /**
     * Class constructor
     */
    testDataGenerator();
    
    /**
     * Class destructor
     */
    ~testDataGenerator();
    
    int getNumReads();
    
    /**
     * Writes a paired read data block to the requested file.
     */
    void writePairedBlocksToFiles(string readsFile1, string readsFile2);
    
    // Read Accessors
    string getForwardReadName(int readNumber);
    string getForwardReadSequence(int readNumber);
    string getForwardReadQuality(int readNumber);
    int getForwardReadLength(int readNumber);
    
    string getReverseReadName(int readNumber);
    string getReverseReadSequence(int readNumber);
    string getReverseReadQuality(int readNumber);
    int getReverseReadLength(int readNumber);
    
    bool getIsMismatched(int readNumber);
    vector<uint32_t> getForwardStartLocs(int readNumber);
    vector<int> getForwardDiffs(int readNumber);
    
    vector<uint32_t> getReverseStartLocs(int readNumber);
    vector<int> getReverseDiffs(int readNumber);
    
  private:
    
    /**
     * Returns a paired read data block in FASTQ format.
     */
    string getPairedReadLine(string name, string sequence, string quality, string suffix);

    void populateReads();
    
    Read_t * reads;
    
};

#endif