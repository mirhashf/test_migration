/**
 * index2_test.cc
 * Tests the index2.cc module of SeqAlto aligner-fast.
 *
 * Written by AJ Minich, Jan 2012
 */

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace std;

#include "index2.h"

// Constants
const string TEST_TARGET = "index2.cc";

const char* const GENOME_DEFAULT = "/mnt/scratch0/public/data/annotation/hg19/chr21.fa";
const char* const KMER_SIZE_DEFAULT = "22";
const char* const VCF_FILE_DEFAULT = "/mnt/scratch0/public/data/variants/dbSNP/CEU/CEU-1409-21.vcf";

// Temp Directory
const char* const TMP_DIR_DEFAULT = "";

// Google Flags
DEFINE_string(test_tmp, TMP_DIR_DEFAULT, "tmp directory with sufficient space to store a small alignment.");
DEFINE_string(test_reference, GENOME_DEFAULT, "Genome file to index.");
DEFINE_string(test_kmer_size, KMER_SIZE_DEFAULT, "k-mer size to use in the index.");
DEFINE_string(test_vcf_file, VCF_FILE_DEFAULT, "Variation call file (VCF) to use for index generation.");

namespace {

    class IndexTest : public ::testing::Test {
      protected:
        // Called before the first test in this test case.
        static void SetUpTestCase() {
            
            // Create index file path
            stringstream ss;
            ss << FLAGS_test_tmp << FLAGS_test_reference << "_" << FLAGS_test_kmer_size << ".midx";
            indexFile = ss.str();
            
            // Check that tmp directory exists
            bool tmp_directory_exists = ( access( FLAGS_test_tmp.c_str(), 0 ) == 0 );
            ASSERT_TRUE (tmp_directory_exists);
        }
    
        virtual ~IndexTest() {
            
        }
    
        virtual void SetUp() {

        }
    
        virtual void TearDown() {
            
            // Remove index file if it exists
            ifstream index_file_exists (indexFile.c_str());
            if (index_file_exists) {
                remove (indexFile.c_str());
            }
        }
        
        static string indexFile;
        
    };
    
    // Shared Resources
    string IndexTest::indexFile;
    
    /**
     * Tests that the get_reads_from_file function works correctly.
     */
    TEST_F(IndexTest, CreateIndex) {

        // Generate the params
        // index: seqalto index [options] <genome_file> <k-mer_size> <output_file> <vcf_file>
        vector<string> params;
        params.push_back("-I");
        params.push_back("0");
        params.push_back(FLAGS_test_reference);
        params.push_back(FLAGS_test_kmer_size);
        params.push_back(indexFile);
        params.push_back(FLAGS_test_vcf_file);
        
        // Call using the same format as main
        int error_num = index2(params);
        
        // Aligner should return an error number of 0 (no error)
        EXPECT_EQ(0, error_num);
        
        // Check that the index file exists
        ifstream index_file_exists (indexFile.c_str());
        ASSERT_TRUE(index_file_exists);
            
        // Check that the index file is of the appropriate size
        struct stat filestatus;
        stat( indexFile.c_str(), &filestatus );
        long fileSize = (long long) filestatus.st_size;
        
        LOG (INFO) << "Final index size is " << fileSize << ".";
    }
    
} // namespace

// Runs the tests contained in the AlignTest class.
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
    LOG (INFO) << "Using genome '" << FLAGS_test_reference << "' as genome reference to index.";
    
    ::testing::InitGoogleTest(&argc, argv);
    
    int result = RUN_ALL_TESTS();
    
    LOG_IF (ERROR, result > 0) << "Unit testing on '" << TEST_TARGET << "' failed.";
    
    return result;
}
