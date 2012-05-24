/**
 * pairing_stage_test.cc
 * Tests the pairing_stage.cc module of SeqAlto aligner-fast.
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

#include "../src/constants.h"
#include "../src/pairing_stage.h"
#include "testUtils.h"

using namespace std;

// Constants
const string TEST_TARGET = "pairing_stage.cc";

// Directories
const char* const TMP_DIR_DEFAULT = "";

// Google Flags
DEFINE_string(test_tmp, TMP_DIR_DEFAULT, "tmp directory with sufficient space to store a small alignment.");
DEFINE_int32(test_mean_template_size, DEFAULT_MEAN_TEMPLATE_SIZE, "The mean template size to use for pairing.");


namespace {

    class PairingStageTest : public ::testing::Test {
      protected:
        
        // Called before the first test in this test case.
        static void SetUpTestCase() {
          
          // Check that tmp directory exists
          bool tmp_directory_exists = ( access( FLAGS_test_tmp.c_str(), 0 ) == 0 );
          ASSERT_TRUE (tmp_directory_exists);
        }
    
        // Called after the last test in this test case.
        static void TearDownTestCase() {
          
        }
    
        virtual void SetUp() {
          
          // Create new SeqAlto params struct
          seqParams = getDefaultSeqParams();
          //seqParams->index_id = FLAGS_index_id;
        
        }
    
        virtual void TearDown() {
          delete seqParams;
        }
        
        seqalto_params_t * seqParams;
        
    };
    
    
    TEST_F(PairingStageTest, DoPairing) {
        
        LOG (WARNING) << "Test not yet implemented.";
        /*
        // Set up thread params
        thread_params_t * threadParams;
        threadParams.p = seqalto_params_t;
        threadParams.tid = 0;
        threadParams.pairing_stage_params.mean_template_size = FLAGS_mean_template_size;
        threadParams.kmer_locs_hist = (unsigned int *) calloc (65536, sizeof(unsigned int));
        
        // Set up single-ended reads data
        single_ended_align_data_t * readAlign1, readAlign2;
        
        fast_t line1[2], line2[2];
        
        
        readAlign1 = new single_ended_align_data_t (line[0], aln_store[x], bad_read[x],
			curr_read_name, total_length[x], total_mega_length[x], line_seq_2bit[x], line_seq_2bit_ug[x], 
			line_len[x], barcode[x]);
        
        
        /* USAGE
        do_pairing(single_ended_align_data_t * read_data0, single_ended_align_data_t * read_data1, 
	vector<concordant_store_t> &concordant_pairs,
	vector<align_loc_pair_t> * temp_insert, stringstream * output_hybrid_ss, 
	pen_tester &pen, seqalto_params_t * p, thread_params_t * tp, MT_random &rand_gen, const int mean_template_size, int tid) {
        */ 
        /*
        pairing_stage_align_data_t * pairingStageResults;
        
        pairingStageResults = do_pairing (readAlign1, readAlign2, concordant_pairs, 
			    temp_insert, output_hybrid_ss, pen, p, tp, rand_gen, mean_template_size, tid);
        
        delete readAlign1;
        delete readAlign2; */
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
    
    ::testing::InitGoogleTest(&argc, argv);
    
    int result = RUN_ALL_TESTS();
    
    LOG_IF (ERROR, result > 0) << "Unit testing on '" << TEST_TARGET << "' failed.";
    
    return result;
}
