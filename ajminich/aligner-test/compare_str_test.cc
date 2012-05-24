/**
 * compare_str_test.cc
 * Tests the compare_str.cc module of SeqAlto aligner-fast.
 *
 * Written by AJ Minich, Jan 2012
 */

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <stdint.h>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace std;

// SIMD
#include <x86intrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <emmintrin.h>

#include "../src/constants.h"
#include "../src/load_index.h"
#include "../src/compare_str.h"
#include "data.h"

// Constants
const string TEST_TARGET = "str_compare.cc";

const char* const INDEX_DEFAULT = "/mnt/scratch0/public/data/annotation/hg19/chr21.fa_22.midx";
const int32_t DEFAULT_MAX_DIFF = 10;
const int32_t DEFAULT_MISMATCH_PENALTY = 1;
const int32_t DEFAULT_MAX_MISMATCH = 10;

const char* const SAMPLE_READ_1 = "ABCDEDFGHIJKLMNOPQRSTUVWXYZ";
const int SAMPLE_READ_1_LENGTH = 26;
const char* const SAMPLE_READ_2 = "ABCDEDFGHIJKLMNOPQRSTUVWYXZ";
const int SAMPLE_READ_2_LENGTH = 26;
const char* const SAMPLE_READ_3 = "MNOPQRSTUVWXYZABCDEDFGHIJKL";
const int SAMPLE_READ_3_LENGTH = 26;

const int EXPECTED_DIFFS_1_2 = 2;
const int SHIFTS_1_3 = 12;
const int EXPECTED_DIFFS_1_3 = 11;

// Directories
const char* const TMP_DIR_DEFAULT = "";

// Google Flags
DEFINE_string(test_tmp, TMP_DIR_DEFAULT, "tmp directory with sufficient space to store a small alignment.");
DEFINE_int32(test_max_template_size, DEFAULT_MAX_TEMPLATE_SIZE, "Max template size to use.");
DEFINE_int32(test_mismatch_penalty, DEFAULT_MISMATCH_PENALTY, "Penalty for each read mismatch.");
DEFINE_int32(test_max_mismatch, DEFAULT_MAX_MISMATCH, "Maximum allowable mismatch.");
DEFINE_string(test_index_file, INDEX_DEFAULT, "Genome index file to load.");


namespace {

    class CompareStrTest : public ::testing::Test {
      protected:
        
        // Called before the first test in this test case.
        static void SetUpTestCase() {
          dataGen = new testDataGenerator();
            
          loadGenome();
        }
    
        // Called after the last test in this test case.
        static void TearDownTestCase() {
            delete dataGen;
        }
    
        virtual void SetUp() {
            
        }
    
        virtual void TearDown() {
        
        }
        
        void checkDiffs(const char* sequence, int length, vector<uint32_t> startLocs,
                        vector<int> expectedDiffs) {
        
          // Get SIMD for the sequence
          __m128i simdRead[length];
          char * testRead = (char *) simdRead;
          strcpy(testRead, sequence);
          
          for (unsigned int index = 0; index < startLocs.size(); index++) {
          
            uint32_t startLoc = startLocs[index];
            uint32_t endLoc = startLoc + FLAGS_test_max_template_size - length;
            
            int actualDiff = mismatch_align(testRead, length, genomeData, genomeLength, startLoc, endLoc, FLAGS_test_max_mismatch);
            int expectedDiff = expectedDiffs[index];
  
            LOG (INFO) << "Difference was " << actualDiff << "; expected to see " << expectedDiff << ".";
            EXPECT_EQ(expectedDiff, actualDiff);
          
          }
          
        }
        
        static testDataGenerator* dataGen;
        
        static char * genomeData;
        static uint64_t genomeLength;
        
      private:
        
        static void loadGenome() {
          LOG (INFO) << "Loading genome hash file.";
          
          // Open hash file and check that it is not null
          gzFile hashFile = gzopen(FLAGS_test_index_file.c_str(), "rb");
          bool hashFileOpened = (hashFile != NULL);
          ASSERT_TRUE(hashFileOpened);        
  
          genomeLength = 0;
          genomeLength = get_genome_length (hashFile);
          LOG (INFO) << "Genome length measured at " << genomeLength << ".";
  
          LOG (INFO) << "Loading genome index from '" << FLAGS_test_index_file << "'.";
  
          // Load the genome and check that the load was successful
          genomeData = NULL;
          int fileLoadResult = read_index_file_genome(hashFile, genomeLength, &genomeData);
          bool genomeLoaded = (fileLoadResult == 0);
          ASSERT_TRUE(genomeLoaded);
          
          LOG (INFO) << "Genome index successfully loaded.";
        }
        
    };
    
    // Shared Resources
    testDataGenerator* CompareStrTest::dataGen = NULL;
    char * CompareStrTest::genomeData;
    uint64_t CompareStrTest::genomeLength;
    
    TEST_F(CompareStrTest, CompareStrRef) {
     
      string read1 (SAMPLE_READ_1);
      string read2 (SAMPLE_READ_2);
      string read3 (SAMPLE_READ_3);
      int actualDiff, expectedDiff;
     
      // Equal read test
      actualDiff = compare_str_ref(read1, read1, 0, FLAGS_test_max_mismatch, FLAGS_test_mismatch_penalty);
      expectedDiff = 0;
      
      LOG (INFO) << "Equal read test: expected " << expectedDiff << "; actual value returned was " << actualDiff << "." << endl;
      EXPECT_EQ(expectedDiff, actualDiff);
     
      // Unequal read test
      actualDiff = compare_str_ref(read1, read2, 0, FLAGS_test_max_mismatch, FLAGS_test_mismatch_penalty);
      expectedDiff = EXPECTED_DIFFS_1_2 * FLAGS_test_mismatch_penalty;
      
      LOG (INFO) << "Unequal read test: expected " << expectedDiff << "; actual value returned was " << actualDiff << "." << endl;
      EXPECT_EQ(expectedDiff, actualDiff);
      
      // Shifted read test
      actualDiff = compare_str_ref(read3, read1, SHIFTS_1_3, FLAGS_test_max_mismatch, FLAGS_test_mismatch_penalty);
      expectedDiff = EXPECTED_DIFFS_1_3 * FLAGS_test_mismatch_penalty;
      
      LOG (INFO) << "Shifted read test: expected " << expectedDiff << "; actual value returned was " << actualDiff << "." << endl;
      EXPECT_EQ(expectedDiff, actualDiff);
      
    }
    
    TEST_F(CompareStrTest, CompareStrUnequal) {
      
      int actualDiff, expectedDiff;
      
      // Equal read test
      actualDiff = compare_str_unequal(SAMPLE_READ_1, SAMPLE_READ_1, SAMPLE_READ_1_LENGTH, FLAGS_test_max_mismatch, FLAGS_test_mismatch_penalty);
      expectedDiff = 0;
      
      LOG (INFO) << "Equal read test: expected " << expectedDiff << "; actual value returned was " << actualDiff << "." << endl;
      EXPECT_EQ(expectedDiff, actualDiff);
      
      // Unequal read test
      actualDiff = compare_str_unequal(SAMPLE_READ_1, SAMPLE_READ_2, SAMPLE_READ_2_LENGTH, FLAGS_test_max_mismatch, FLAGS_test_mismatch_penalty);
      expectedDiff = EXPECTED_DIFFS_1_2 * FLAGS_test_mismatch_penalty;
      
      LOG (INFO) << "Unequal read test: expected " << expectedDiff << "; actual value returned was " << actualDiff << "." << endl;
      EXPECT_EQ(expectedDiff, actualDiff);
    }
    
    TEST_F(CompareStrTest, CompareStrSimd) {
      
      // Get SIMD pointers for test reads
      __m128i simdRead1[SAMPLE_READ_1_LENGTH];
      char * sampleRead1 = (char *) simdRead1;
      strcpy(sampleRead1, SAMPLE_READ_1);
      
      __m128i simdRead2[SAMPLE_READ_2_LENGTH];
      char * sampleRead2 = (char *) simdRead2;
      strcpy(sampleRead2, SAMPLE_READ_2);
      
      int actualDiff, expectedDiff;
      
      // Equal read test
      actualDiff = compare_str_simd(sampleRead1, SAMPLE_READ_1_LENGTH, sampleRead1, FLAGS_test_max_mismatch, FLAGS_test_mismatch_penalty);
      expectedDiff = 0;
      
      LOG (INFO) << "Equal read test: expected " << expectedDiff << "; actual value returned was " << actualDiff << "." << endl;
      EXPECT_EQ(expectedDiff, actualDiff);
      
      // Unequal read test
      actualDiff = compare_str_simd(sampleRead1, SAMPLE_READ_1_LENGTH, sampleRead2, FLAGS_test_max_mismatch, FLAGS_test_mismatch_penalty);
      expectedDiff = EXPECTED_DIFFS_1_2 * FLAGS_test_mismatch_penalty;
      
      LOG (INFO) << "Unequal read test: expected " << expectedDiff << "; actual value returned was " << actualDiff << "." << endl;
      EXPECT_EQ(expectedDiff, actualDiff);
      
    }
    
    TEST_F(CompareStrTest, MismatchAlign) {
      
      // Repeat test for all available mismatched reads
      for (int index = 0; index < dataGen->getNumReads(); index++) {
        
        // If the read is supposed to be a mismatch, run the test
        if (dataGen->getIsMismatched(index)) {
         
          int length = dataGen->getForwardReadLength(index);
          const char * sequence = dataGen->getForwardReadSequence(index).c_str();
          
          vector<uint32_t> forwardStartLocs = dataGen->getForwardStartLocs(index);
          vector<int> forwardDiffs = dataGen->getForwardDiffs(index);
          
          checkDiffs(sequence, length, forwardStartLocs, forwardDiffs);
          
          vector<uint32_t> reverseStartLocs = dataGen->getReverseStartLocs(index);
          vector<int> reverseDiffs = dataGen->getReverseDiffs(index);
          
          checkDiffs(sequence, length, reverseStartLocs, reverseDiffs);
        }
              
      }
        
    }
    
} // namespace

// Runs the tests contained in the LoadIndexTest class.
int main(int argc, char **argv) {
    
    // Initialize logging
    google::InitGoogleLogging (argv [0]);
    google::ParseCommandLineFlags (&argc, &argv, true);
    
    LOG (INFO) << "Starting unit tests on '" << TEST_TARGET << "'.";
    
    // Check that last character of location is a backslash; if not, add it
    if (FLAGS_test_tmp.find_last_of('/') != FLAGS_test_tmp.length() - 1) {
        FLAGS_test_tmp.append("/");
    }
    
    LOG (INFO) << "Using directory '" << FLAGS_test_tmp << "' as temp directory.";
    LOG (INFO) << "Using genome index '" << FLAGS_test_index_file << "' as reference.";
    
    ::testing::InitGoogleTest(&argc, argv);
    
    int result = RUN_ALL_TESTS();
    
    LOG_IF (ERROR, result > 0) << "Unit testing on '" << TEST_TARGET << "' failed.";
    
    return result;
}
