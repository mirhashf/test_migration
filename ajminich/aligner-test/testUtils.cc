/**
 * testUtils.cc
 * Defines several utilities common to test functionality.
 *
 * Written by AJ Minich, Jan 2012
 */

#include "testUtils.h"

using namespace std;

seqalto_params_t * getDefaultSeqParams() {
    
    // Set up seqalto_params_t object
    seqalto_params_t* seqParams = new seqalto_params_t();
    
    seqParams->hash_filename = "";

    seqParams->INDEX_MODE  = 1;

    seqParams->end_gap_dist = DEFAULT_END_GAP_DIST; // double penalty for gaps close to end
    seqParams->barcode_len = 0;

    seqParams->MAX_NUM_HASH = DEFAULT_MAX_NUM_HASH; // remove if bigger than this
    seqParams->PHASE_1_LENGTH = DEFAULT_PHASE_1_LENGTH; // 
    seqParams->PHASE_1_LENGTH_SW = DEFAULT_PHASE_1_LENGTH_SW; //
    seqParams->PHASE_2_LENGTH_SW = DEFAULT_PHASE_2_LENGTH_SW; //

    // Set up penalty tester
    int max_ext = DEFAULT_MAX_EXT;
    int m_score = DEFAULT_M_SCORE;
    int mm_pen  = DEFAULT_MM_PEN;
    int gap_pen = DEFAULT_GAP_PEN;
    int ext_pen = DEFAULT_EXT_PEN;
    int max_mm = -1;
    float mm_prob = DEFAULT_MM_PROB;
    int   max_gap  = -1;
    float gap_rate = DEFAULT_GAP_RATE;

    seqParams->pen = new pen_tester (mm_prob, max_mm, gap_rate, max_gap, max_ext, m_score, mm_pen, gap_pen, ext_pen);

    seqParams->M_SW = DEFAULT_M_SW;
    seqParams->S_SW = DEFAULT_S_SW;
    seqParams->D_SW = DEFAULT_D_SW;
    seqParams->E_SW = DEFAULT_E_SW;

    seqParams->SKIP_AMOUNT = 2;

    seqParams->MAX_HITS = 10;

    seqParams->LOOK_AHEAD_NUM = DEFAULT_LOOK_AHEAD_NUM; // so we can find the sub optimal hits
    seqParams->NUM_VAL = DEFAULT_NUM_VAL;
    seqParams->DENUM_VAL = 1;

    seqParams->LOOK_AHEAD_NUM_SW = DEFAULT_LOOK_AHEAD_NUM_SW;
    seqParams->NUM_VAL_SW = DEFAULT_NUM_VAL_SW;
    seqParams->DENUM_VAL_SW = DEFAULT_DENUM_VAL_SW;

    seqParams->NUM_THREADS = 1;

    seqParams->BUFFER_SIZE = DEFAULT_BUFFER_SIZE; // buffer this many reads 

    // paired end options
    seqParams->MAX_TEMPLATE_SIZE = DEFAULT_MAX_TEMPLATE_SIZE;
    seqParams->MEAN_TEMPLATE_SIZE = DEFAULT_MEAN_TEMPLATE_SIZE;
    seqParams->PHRED_PAIRING_PRIOR = DEFAULT_PHRED_PAIRING_PRIOR;
    seqParams->MIN_PERCENT_READ = DEFAULT_MIN_PERCENT_READ;

    // hybrid mode options
    seqParams->hybrid_mode = false;
    seqParams->mapq_reject = DEFAULT_MAPQ_REJECT;
    seqParams->hybrid_max_gap_reject = DEFAULT_HYBRID_MAX_GAP_REJECT; // 2 is selected for BWA
    seqParams->read_min_length = DEFAULT_READ_MIN_LENGTH;

    seqParams->ungapped   = false;
    seqParams->paired     = true;
    seqParams->no_pair_sw = false;

    seqParams->trim_reads  = false;
    seqParams->trim_qual   = DEFAULT_TRIM_QUAL;
    seqParams->qual_format = qual::SANGER;

    seqParams->output_RG = false;
    seqParams->output_PU = false;
    seqParams->output_PG = false;

    seqParams->cigar14 = false;

    seqParams->kmer_trim = DEFAULT_KMER_TRIM;

    seqParams->RG_ID = "";
    seqParams->PU    = "";

    seqParams->use_leveldb = false;

    seqParams->use_unknown_qual_bp_as_dont_care = false;
        
    seqParams->enable_fpga_accel = false;
    seqParams->enable_monitor_thread = false;
    seqParams->enable_server_thread = false;
    seqParams->disable_single_gap = false;
    seqParams->mix_ungapped_gapped = false;
    seqParams->enable_batch_pairing = false;
    seqParams->use_unknown_qual_bp_as_dont_care = false;
    seqParams->template_len_comp_method = 0;
    seqParams->index_id = 0;
    
    // zero progress variables
    seqParams->num_CS = 0;
    seqParams->num_NW = 0;
    seqParams->num_NW_SSE = 0;
    seqParams->num_NW_pair = 0;
    seqParams->num_SW = 0;
    seqParams->num_SW_SSE = 0;
    seqParams->num_MA = 0;
    seqParams->num_CS_outer = 0;
    seqParams->num_kmers_low_locs = 0;
    seqParams->num_kmers_gapped_low_locs = 0;
    seqParams->num_gapped = 0;

    seqParams->num_processed = 0;
    seqParams->num_aligned = 0;
    seqParams->num_done = 0;
    seqParams->sequence_counter = 1; // this is for leveldb
    seqParams->kmer_locs_hist = (unsigned int *) calloc (65536, sizeof(unsigned int));

    if (seqParams->paired) {
        seqParams->num_fragment = 2;
    } else {
        seqParams->num_fragment = 1;
    }

    // Set up threading mutexes
    pthread_mutex_t* read_lock_temp = new pthread_mutex_t();
    seqParams->read_lock = read_lock_temp;
    pthread_mutex_init(seqParams->read_lock, NULL);

    pthread_mutex_t* write_lock_temp = new pthread_mutex_t();
    seqParams->write_lock = write_lock_temp;
    pthread_mutex_init(seqParams->write_lock, NULL);

    pthread_mutex_t* mouse_lock_temp = new pthread_mutex_t();
    seqParams->mouse_lock = mouse_lock_temp;
    pthread_mutex_init(seqParams->mouse_lock, NULL);
    
    return seqParams;
}
