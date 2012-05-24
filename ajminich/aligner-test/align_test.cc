/**
 * align_test.cc
 * Tests the align.cc module of SeqAlto aligner-fast.
 *
 * Written by AJ Minich, Jan 2012
 */

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <sstream>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

using namespace std;

#include "../src/constants.h"
#include "../src/align.h"
#include "data.h"
#include "testUtils.h"

// Switch to enable running or not running the aligner, which currently generates
// a ton of messages to stderr and stdout.
//#define ENABLE_ALIGN

// Constants
const string TEST_TARGET = "align.cc";

const char* const INDEX_DEFAULT = "/mnt/scratch0/public/data/annotation/hg19/chr21.fa_22.midx";
const int INDEX_ID = 1337;
const char* const READS_FILE_NAME_1 = "test_reads_1.fq";
const char* const READS_FILE_NAME_2 = "test_reads_2.fq";

// Temp Directory
const char* const TMP_DIR_DEFAULT = "";

// Google Flags
DEFINE_string(test_tmp, TMP_DIR_DEFAULT, "tmp directory with sufficient space to store a small alignment.");
DEFINE_string(test_index_file, INDEX_DEFAULT, "Genome index file to use for the aligner unit test.");
DEFINE_int32(test_mean_template_size, DEFAULT_MEAN_TEMPLATE_SIZE, "The mean template size to use for pairing.");

// const char* ALIGN_TMP = "testAlignment.sam";
const char* NUM_THREADS = "16";
const int BUFFER_SIZE = DEFAULT_BUFFER_SIZE;

namespace {

    class AlignTest : public ::testing::Test {
      protected:
        
        // Called before the first test in this test case.
        static void SetUpTestCase() {
          
            dataGen = new testDataGenerator();
            
            // Check that tmp directory exists
            bool tmp_directory_exists = ( access( FLAGS_test_tmp.c_str(), 0 ) == 0 );
            ASSERT_TRUE (tmp_directory_exists);
            
            // Generate testÊreads files
            reads1.assign(FLAGS_test_tmp);
            reads1.append(READS_FILE_NAME_1);
            
            reads2.assign(FLAGS_test_tmp);
            reads2.append(READS_FILE_NAME_2);
            
            dataGen->writePairedBlocksToFiles(reads1, reads2);
            
            // Check that the reads files exist
            ifstream reads1_file_exists (reads1.c_str());
            ASSERT_TRUE(reads1_file_exists);
            
            ifstream reads2_file_exists (reads2.c_str());
            ASSERT_TRUE(reads2_file_exists);
        }
      
        // Called after the last test in this test case.
        static void TearDownTestCase() {
            delete dataGen;
          
            // Remove reads files
            remove( reads1.c_str() );
            remove( reads2.c_str() );
        }
    
        virtual void SetUp() {
            
            // Create new SeqAlto params struct
            seqParams = getDefaultSeqParams();
            seqParams->index_id = INDEX_ID;
            //seqalto_params_t * seqParams = new seqalto_params_t();
            //parse_params(seqParams);
        }
    
        virtual void TearDown() {
            delete seqParams;
        }
        
        static testDataGenerator* dataGen;
        
        // Reads file names
        static string reads1, reads2;
        
        seqalto_params_t * seqParams;
    };
    
    // Shared Resources
    testDataGenerator* AlignTest::dataGen = NULL;
    string AlignTest::reads1;
    string AlignTest::reads2;
    
    /**
     * Tests that the get_reads_from_file function works correctly.
     */
    TEST_F(AlignTest, GetReadsFromFile) {

        // Set up reads files
        fast * readsFile1 = new fast();
        fast * readsFile2 = new fast();
        fast * readsFiles[2] = { readsFile1, readsFile2 };
        
        readsFiles[0]->open(reads1);
        readsFiles[1]->open(reads2);
    
        EXPECT_TRUE (readsFiles[0]->is_open());
        EXPECT_TRUE (readsFiles[1]->is_open());
            
        fast_t * lineBuffer[2] = { NULL, NULL };
        lineBuffer[0] = new fast_t[BUFFER_SIZE];
        lineBuffer[1] = new fast_t[BUFFER_SIZE];
            
        int numReads = get_reads_from_file(readsFiles, lineBuffer, BUFFER_SIZE, true);
    
        // Check that we have the correct number of reads
        ASSERT_EQ(dataGen->getNumReads(), numReads);
        
        fast_t line1, line2;
        
        for (int readIndex = 0; readIndex < numReads; readIndex++) {
    
            // Check forward read
            line1 = lineBuffer[0][readIndex];
            LOG (INFO) << "Read " << readIndex << " forward: '" << line1.name
              << "' = " << line1.seq << " (map quality = " << line1.qual << ").";
    
            EXPECT_EQ(line1.name, dataGen->getForwardReadName(readIndex));
            EXPECT_EQ(line1.seq, dataGen->getForwardReadSequence(readIndex));
            EXPECT_EQ(line1.qual, dataGen->getForwardReadQuality(readIndex));
    
            // Check reverse read
            line2 = lineBuffer[1][readIndex];   // paired read
            LOG (INFO) << "Read " << readIndex << " reverse: '" << line2.name
              << "' = " << line2.seq << " (map quality = " << line2.qual << ").";
    
            EXPECT_EQ(line2.name, dataGen->getReverseReadName(readIndex));
            EXPECT_EQ(line2.seq, dataGen->getReverseReadSequence(readIndex));
            EXPECT_EQ(line2.qual, dataGen->getReverseReadQuality(readIndex));
            
        }
        
    }

#ifdef ENABLE_ALIGN

    /*
     * Tests that the align_thread subfunction works correctly.
     */
    TEST_F(AlignTest, AlignThread) {
        
        // File setup
        fast * fastReadsFile;
        
        seqParams->reads_filename[0] = reads1;
        fastReadsFile = new fast();
        seqParams->reads_file[0] = fastReadsFile;
        seqParams->reads_file[0]->open(reads1);
        ASSERT_TRUE(seqParams->reads_file[0]->is_open());
        
        seqParams->reads_filename[1] = reads2;
        fastReadsFile = new fast();
        seqParams->reads_file[1] = fastReadsFile;
        seqParams->reads_file[1]->open(reads2);
        ASSERT_TRUE(seqParams->reads_file[1]->is_open());
        
        // Load index
        unsigned int options = LOAD_GENOME_FLAG | LOAD_INDEX_FLAG | LOAD_GENOME_START_FLAG |
		LOAD_GENOME_UPPERCASE;
      
        int errorNum = load_index_into_vars(FLAGS_test_index_file, seqParams, options);
        EXPECT_EQ(0, errorNum);
        
        LOG_IF (INFO, errorNum == 0) << "Genome load into memory was successful.";
        
        // Initialize libSW
        libsw_init(1, seqParams->enable_fpga_accel);
        set_genome ((const uint8_t *) seqParams->genome_c_str, seqParams->genome_length, ACTG_chars, 0);
        set_penalties_nw(seqParams->pen->get_m_score(), seqParams->pen->get_mm_pen(), seqParams->pen->get_gap_pen() + seqParams->pen->get_ext_pen(), 
                seqParams->pen->get_ext_pen(), seqParams->end_gap_dist);
        set_penalties_sw(seqParams->M_SW, seqParams->S_SW, seqParams->D_SW, seqParams->E_SW); 
        
        LOG (INFO) << "LibSW successfully initialized.";
        
        // Set up thread params
        thread_params_t * threadParams = (thread_params_t *) calloc(1, sizeof(thread_params_t));;
        threadParams[0].p = seqParams;
        threadParams[0].tid = 0;
        threadParams[0].pairing_stage_params.mean_template_size = FLAGS_test_mean_template_size;
        threadParams[0].kmer_locs_hist = (unsigned int *) calloc (65536, sizeof(unsigned int));
        
        LOG (INFO) << "Starting single aligner thread.";
        
        // Run the aligner
        align_thread((void*) threadParams);
        
        // Close the reads files
        seqParams->reads_file[0]->close();
        seqParams->reads_file[1]->close();
    
        if (seqParams->hybrid_mode) {
            seqParams->hybrid_out[0]->close();
            seqParams->hybrid_out[1]->close();
        }
        
        LOG (INFO) << "Aligner thread succeeded.";
        
    }
    
#endif

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
    LOG (INFO) << "Using genome index '" << FLAGS_test_index_file << "' as reference.";
    
    ::testing::InitGoogleTest(&argc, argv);
    
    int result = RUN_ALL_TESTS();
    
    LOG_IF (ERROR, result > 0) << "Unit testing on '" << TEST_TARGET << "' failed.";
    
    return result;
}
