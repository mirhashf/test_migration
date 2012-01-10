//  Created by Amirhossein Kiani on 8/10/11.
//  Copyright 2011 Bina Technologies. All rights reserved.

#include "leveldb/db.h"
#include "leveldb/slice.h"
#include "leveldb/options.h"
#include "leveldb/comparator.h"
#include "leveldb/iterator.h"
#include "leveldb/write_batch.h"
#include "leveldb/status.h"
#include "leveldb/env.h"
#include "leveldb/table.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>

using namespace std;


class TwoPartComparator : public leveldb::Comparator {
public:
    // Three-way comparison function:
    //   if a < b: negative result
    //   if a > b: positive result
    //   else: zero result
    int Compare(const leveldb::Slice& a, const leveldb::Slice& b) const {
        uint64_t a1, b1;
        
        a1 = *((int*)a.data()); 
        b1 = *((int*)b.data());

        if (a1 < b1) return -1;
        if (a1 > b1) return +1;

        return 0;
    }
    
    // Ignore the following methods for now:
    virtual const char* Name() const { return "TwoPartComparator"; }
    void FindShortestSeparator(std::string*, const leveldb::Slice&) const { }
    void FindShortSuccessor(std::string*) const { }
    
};


// This class is specific to linux :(
class mu_timer{
private:
	struct timeval start_tv;
public:
    
    mu_timer(){
        gettimeofday(&start_tv, NULL);
    }
    
    void reset(){
        gettimeofday(&start_tv, NULL);
    }
    
    double elapsed_time(){
        struct timeval tv;
        gettimeofday(&tv, NULL);
        
        return (double)(tv.tv_sec - start_tv.tv_sec) + (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
    }
};


////////



void read(leveldb::DB* db){
    leveldb::Iterator* it = db->NewIterator(leveldb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        cout << it->key().ToString() << ": "  << it->value().ToString() << endl;
    }
    assert(it->status().ok());  // Check for any errors found during the scan
    delete it;
}


void atomic_write(leveldb::DB* db){
    leveldb::WriteBatch batch;
    srand (123456 );
    string value = "1000000;chrM:8901;5	0	chrM	8902	43	100M	*	0	0 AGCCCACTTCTTACCACAAGGCACACGTACACCCCTTATCCCCATACTAGTTATTATCGTAACCATCAGCCTACTCATACAACCAATAACCGTGGCCGTAIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:5	XM:i:5	XO:i:0	XG:i:0	X0:i:1	X1:i:0	XL:i:25	XR:i:0	XT:A:U";
    char key[50];
    int i = 0;
    mu_timer timer;
    while (i<5000){
        sprintf(key, "%d", rand());
        batch.Put(key, value);
        i++;
        //if(i%10000) printf("added %d\n", i);
    }
    printf("%f,", timer.elapsed_time());    
    leveldb::Status status  = db->Write(leveldb::WriteOptions(), &batch);
}


void ineficient_write(leveldb::DB* db){
    srand (123456 );
    string value = "1000000;chrM:8901;5	0	chrM	8902	43	100M	*	0	0 AGCCCACTTCTTACCACAAGGCACACGTACACCCCTTATCCCCATACTAGTTATTATCGTAACCATCAGCCTACTCATACAACCAATAACCGTGGCCGTAIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:5	XM:i:5	XO:i:0	XG:i:0	X0:i:1	X1:i:0	XL:i:25	XR:i:0	XT:A:U";
    char key[50];
    int i = 0;
    mu_timer timer;
    while (i<1000000){
        sprintf(key, "%d", rand());
        db->Put(leveldb::WriteOptions(), key,value);
        i++;
        //if(i%10000) printf("added %d\n", i);
    }
    printf("%f,",timer.elapsed_time());
}

int main(){
    TwoPartComparator cmp;
    leveldb::DB* db;
    leveldb::Options options;
    options.create_if_missing = true;
    options.comparator = &cmp;
    leveldb::Status status = leveldb::DB::Open(options, "/tmp/testdb-two", &db);
    mu_timer timer;
    if (!status.ok()) cout << status.ToString() << endl;
    atomic_write(db);
    delete(db);
    printf("%f\n",timer.elapsed_time());    
}






