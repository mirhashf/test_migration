/**
 * data.cc
 * Contains data generators for alignment unit testing.
 *
 * Written by AJ Minich, Jan 2012
 */

#include "data.h"

using namespace std;

testDataGenerator::testDataGenerator() {
    
    // Populate data
    populateReads();
}

testDataGenerator::~testDataGenerator() {
    
    // Delete data structures
    delete [] reads;
}

int testDataGenerator::getNumReads() {
    return NUM_READS;
}

/**
 * Writes the paired read data to the specified FASTQ files.
 */
void testDataGenerator::writePairedBlocksToFiles(string readsFile1, string readsFile2) {

    ofstream readsFileStream1, readsFileStream2;
    
    // Open dual filestreams
    readsFileStream1.open (readsFile1.c_str());
    readsFileStream2.open (readsFile2.c_str());
    
    for (int readIndex = 0; readIndex < getNumReads(); readIndex++) {
        
        Read_t currentRead = reads[readIndex];
        
        readsFileStream1 << getPairedReadLine(currentRead.name, currentRead.forwardSequence, currentRead.forwardQuality, "1");
        readsFileStream2 << getPairedReadLine(currentRead.name, currentRead.reverseSequence, currentRead.reverseQuality, "2");
    }
    
    // Close filestreams
    readsFileStream1.close();
    readsFileStream2.close();
}

// Accessor Functions
string testDataGenerator::getForwardReadName(int readNumber) {
    return reads[readNumber].name;
}

string testDataGenerator::getForwardReadSequence(int readNumber) {
    return reads[readNumber].forwardSequence;
}

string testDataGenerator::getForwardReadQuality(int readNumber) {
    return reads[readNumber].forwardQuality;
}

int testDataGenerator::getForwardReadLength(int readNumber) {
    return reads[readNumber].forwardLength;
}

string testDataGenerator::getReverseReadName(int readNumber) {
    return reads[readNumber].name;
}

string testDataGenerator::getReverseReadSequence(int readNumber) {
    return reads[readNumber].reverseSequence;
}

string testDataGenerator::getReverseReadQuality(int readNumber) {
    return reads[readNumber].reverseQuality;
}

int testDataGenerator::getReverseReadLength(int readNumber) {
    return reads[readNumber].reverseLength;
}

bool testDataGenerator::getIsMismatched(int readNumber) {
    return reads[readNumber].isMismatched;
}

vector<uint32_t> testDataGenerator::getForwardStartLocs(int readNumber) {
    return reads[readNumber].forwardStartLocs;
}

vector<int> testDataGenerator::getForwardDiffs(int readNumber) {
    return reads[readNumber].forwardDiffs;
}

vector<uint32_t> testDataGenerator::getReverseStartLocs(int readNumber) {
    return reads[readNumber].reverseStartLocs;
}

vector<int> testDataGenerator::getReverseDiffs(int readNumber) {
    return reads[readNumber].reverseDiffs;
}

/** PRIVATE METHODS */

/**
 * Returns a short paired read data block for forward reads in FASTQ format.
 *
 * The suffix is typically the end that the read is coming from.
 * - For reads from the primary end, the suffix is "1".
 * - For reads from the paired end, the suffix is "2".
 */
string testDataGenerator::getPairedReadLine(string name, string sequence, string quality, string suffix) {

    stringstream ss;
    ss << "@" << name << "/" << suffix << endl;
    ss << sequence << endl;
    ss << "+" << endl;
    ss << quality << endl;
    return ss.str();
}

void testDataGenerator::populateReads() {
    
    // Forward Reads
    reads = new Read_t[NUM_READS];
    
    reads[0].name = "1pa;chr1:501295762;501295996;334;0";
    reads[0].forwardSequence = "GACAATTTGCTCCAAGCCCTGGGTGCTGTTGGTGCAGCCAACTCTTTTCTCTGCTGCTATTAGTTTGACAATTTCTGATATGTCATATACATGGAATCAA";
    reads[0].forwardQuality = "GGFGGGGGGGGGGAGFGGGGGGE?GGGGGGG=FGGGGGGGGGFFDGCBDGGGGGG?ACGE>GGDGGGG.GGG@GCGGGGGFGGG<AGGD4=@-CGG<7GG";
    reads[0].forwardLength = 100;
    reads[0].reverseSequence = "TCATTTAACTACAAGTGTCAAACCAATGAACAATGGATAGAAACTTTCTGTTAGATGAATGGATAAAGAAAACACAGTAAATACATACAGTGGAATATTA";
    reads[0].reverseQuality = "GGGGGGGCDGGGDGGGGFGGGGGGG@GG?GFGFGGGGGGGG@GGGGGGG=GEG=GGGGGGGGAFGG?GFGGG2C%GGDGEDGB==GG:CGG>?C/GG!G?";
    reads[0].reverseLength = 100;
    reads[0].isMismatched = false;
    
    reads[1].name = "2pa;chr19_gl000209_random:411629069;411629307;338;1";
    reads[1].forwardSequence = "CTCCAATAATAAGATAAAGCCAATAAAGTATTAGTACACTAAAATGCATGCTCTTCAGTCTTCCACATTTATTAAAAGTAGTGGCGAGCTAAACACTAAA";
    reads[1].forwardQuality = "GGGGGAGGDGGDGGGGGGDBGGGG@GGAGGGGGGFGGGGGGCGFGGEE<G?GGG?GG><G6/BGG6GGGGFGFGGGGBGEGGGGGGAGE!5C>3;B>G4B";
    reads[1].forwardLength = 100;
    reads[1].reverseSequence = "CTCAGCCTCCTGCGTAGCTGGCATTACAGGCGCCCACCACCATGCCTGGCTAATTTTTGTATTTTTAGTAGAGACGGGGTTTCACCATGTTGGCCGGGCT";
    reads[1].reverseQuality = "GGCGFGGGGGGGGGGGGGGGG>FGGGGGDEGGG<GGGGGBGGGGBCGC:EFGGGGG7G5DGFGGFG@GCDG;8BAG!GGGGGGGG*GG3GGDB7GGA82E";
    reads[1].reverseLength = 100;
    reads[1].isMismatched = false;
    
    reads[2].name = "3pa;chr19_gl000209_random:39204063;39204283;320;1";
    reads[2].forwardSequence = "GGTGATTTATAAAGAAAAAGAAGTTTAACAGATTCACATTTCGACATGGCTGGGGAGGCCTCATAATCATGGCAGAAGGCAAAGGAGGAGCAAAGCCATG";
    reads[2].forwardQuality = "GGGG@GGGGGGGGBG>GGGGGGBGGGGGGGBGGGGEGGCGGDGG?=FGEGGDGGGGBEG<DFFGG0GGG>GGGGGGG,GF;GG:GG:.E);B/>3GGG5#";
    reads[2].forwardLength = 100;
    reads[2].reverseSequence = "AGTGATATGGTTTAGCTGTGTCCCCACCCAAATTTCATCTTGAACTGTAGCCCCCATAATCCCCACATGTCATGGGAGGGACCCAGTGGGAGGTAACTGA";
    reads[2].reverseQuality = "GGGGFGGGFGGGGBGGGGGGEGGG@FGGGGGGGGGGGGGGGGGGGAG>G=GDGGGGGG1DGGF?CGGA-BG;G77DAG<:=?9GG<<AGGGG!G!G)G@7";
    reads[2].reverseLength = 100;
    reads[2].isMismatched = false;
    
    reads[3].name = "505pa;chr19_gl000209_random:752881921;752882144;323;1";
    reads[3].forwardSequence = "AGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCACGGTGTATATGTGCCACCTTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTG";
    reads[3].forwardQuality = "GGCGFF?GGGGGGGGGGGGEGBGGGG@GGGGGGGGGGFGGGEGG?GGCGGGGA8GGFGGGCC9GBGF>GGGG:GGCGDGG?FG.C:G7:GFGG?GG22CG";
    reads[3].forwardLength = 100;
    reads[3].reverseSequence = "AGAAATAGGAACACTTTTACACTGTTGGTGGGACTGTAAACTAGTTCCACCATTGTGGAAGTCAGTGTGGCCATTCCTCAGGGATCTAGAACTCGAAATA";
    reads[3].reverseQuality = "FEEEGGGGGGFGGGGGGDAGGGGDAGGGGGEGB@GBGGGGGGCGGGGGAGGGGGGGGAGGGGGGGGFG>GGAG6G66GG7GGG<GGG>GG@71)3!G89!";
    reads[3].reverseLength = 100;
    
    // Mismatch Read
    // NOTE: this read is a mismatch only if the hg19 chr21 genome is used.
    reads[3].isMismatched = true;
    reads[3].forwardStartLocs.push_back(14521475);
    reads[3].forwardDiffs.push_back(-1);
    reads[3].forwardStartLocs.push_back(26407107);
    reads[3].forwardDiffs.push_back(-1);
    reads[3].forwardStartLocs.push_back(26406877);
    reads[3].forwardDiffs.push_back(-1);
    
    reads[3].reverseStartLocs.push_back(14521475);
    reads[3].reverseDiffs.push_back(-1);
    reads[3].reverseStartLocs.push_back(26407107);
    reads[3].reverseDiffs.push_back(-1);
    
    //reads[3].reverseStartLocs.push_back(14789589);
    //reads[3].reverseDiffs.push_back(5);
    
}