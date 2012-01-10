/* 
 * File:   main.cpp
 * Author: johnmu
 *
 * Created on August 25, 2011, 11:54 AM
 */


#include "main.h"

/*
 * 
 */

void write_line_buffer(vector<string> &line_buffer, MT_random &mt,params_t* p) {
    // generate some random numbers

    vector<uint32_t> rand_num; // these are the random locations for each read

    for (uint64_t i = 0; i < line_buffer.size(); i++) {
        rand_num.push_back((uint32_t) mt.genrand_int_range(1, 3200000000));
        
        //for(int i = 0;i<100;i++){
        //    uint64_t *mouse = new uint64_t();
        //    *mouse = mt.genrand64_int64();
        //    delete mouse;
        //}
        
    }



    // write to leveldb


    leveldb::WriteBatch batch;

    level_key_t key;


    vector<string>::iterator line_it = line_buffer.begin();
    vector<uint32_t>::iterator num_it = rand_num.begin();


    //timer.reset();
    pthread_mutex_lock(p->write_lock);
    while (line_it != line_buffer.end()) {

        key.loc = (*num_it - 1); // minus 1 for sorting 0 to bottom
        key.key = p->sequence_counter;


        leveldb::Slice s((char*) (&key), sizeof (level_key_t));


        batch.Put(s, *line_it);


        p->sequence_counter++;
        line_it++;
        num_it++;
    }
    pthread_mutex_unlock(p->write_lock);
    
    //cout << timer.elapsed_time();


    p->db->Write(leveldb::WriteOptions(), &batch);

    //cout << "\t" << timer.elapsed_time() << "\n";

    line_buffer.clear();
}


void* worker(void* params) {
    params_t* p = (params_t*) params;
    
    
    string line;
    MT_random mt;
    
    bool done = false;
    vector<string> line_buffer;

    while (!done) {

        pthread_mutex_lock(p->read_lock);
        while (!p->infile->eof() && line_buffer.size() < p->batch_size ) {
            getline(*(p->infile), line);

            trim2(line);
            if (line.length() == 0) {
                continue;
            }

            line_buffer.push_back(line);
            

        }
        pthread_mutex_unlock(p->read_lock);
        
        
        // write out last bit
        if (line_buffer.size() > 0) {
            write_line_buffer(line_buffer, mt, p);
        }else{
            done = true;
        }

    }

    return 0;
}


int main(int argc, char** argv) {

    params_t *p = new params_t();
    
    string usage_str = "ldbprof <db_path> <read_file> <write_buffer> <block_cache> <block_size> <batch_size> <num_threads>\n buffer and cache are in MB, block size in KB";
    string temp = "";
    string line = "";

    if (argc != 8) {
        cerr << usage_str << '\n';
        return 1;
    }

    chrLocComp cmp; // this should be in the heap.... :S

    string db_path = argv[1];
    string read_path = argv[2];

    temp = argv[3];
    uint64_t wb_size = strTo<uint64_t > (temp);

    temp = argv[4];
    uint64_t block_cache = strTo<uint64_t > (temp);

    temp = argv[5];
    uint64_t block_size = strTo<uint64_t > (temp);

    temp = argv[6];
    p->batch_size = strTo<uint64_t > (temp);
    
    temp = argv[7];
    uint64_t num_threads = strTo<uint64_t > (temp);

    

    leveldb::Options options;
    options.create_if_missing = true;
    options.comparator = &cmp;
    options.write_buffer_size = wb_size * 1048576;
    options.block_cache = leveldb::NewLRUCache(block_cache * 1048576); // 1.5GB cache
    options.block_size = block_size * 1024;

    leveldb::Status db_status = leveldb::DB::Open(options, db_path.c_str(), &(p->db));

    if (!db_status.ok()) cerr << db_status.ToString() << '\n';

    p->infile = new ifstream(read_path.c_str());

    if (!p->infile->is_open()) {
        cerr << "Could not open: " << read_path << '\n';
        return 1;
    }


    mu_timer all_timer;
    
    
    
    p->sequence_counter = 1;
    
    
    
    pthread_t* t_group = new pthread_t[num_threads]; 
    pthread_mutex_t* read_lock_temp = new pthread_mutex_t();
    p->read_lock = read_lock_temp;
    pthread_mutex_init(p->read_lock, NULL);

    pthread_mutex_t* write_lock_temp = new pthread_mutex_t();
    p->write_lock = write_lock_temp;
    pthread_mutex_init(p->write_lock, NULL);
    

    all_timer.reset();
    for (uint64_t k = 0; k < num_threads; k++) {
        pthread_create(&(t_group[k]), NULL, worker, (void*) p);
    }


    void* output;
    for (uint64_t k = 0; k < num_threads; k++) {

        pthread_join(t_group[k], &output);
    }

    cout << "Total time: " << all_timer.elapsed_time() << "s" << '\n';


    //////
    // Delete the new stuff

    delete(p->db);


    
    
    
    

    return 0;
}

