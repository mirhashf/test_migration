/* 
 * File:   main.h
 * Author: johnmu
 *
 * Created on August 25, 2011, 12:53 PM
 */

#ifndef MAIN_H
#define	MAIN_H

#ifdef	__cplusplus
extern "C" {
#endif




#ifdef	__cplusplus
}
#endif


#include <cstdlib>

#include "leveldb/db.h"
#include "leveldb/slice.h"
#include "leveldb/options.h"
#include "leveldb/comparator.h"
#include "leveldb/iterator.h"
#include "leveldb/write_batch.h"
#include "leveldb/status.h"
#include "leveldb/env.h"
#include "leveldb/table.h"
#include "leveldb/cache.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <ctime>
#include <map>
#include <sstream>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <assert.h>
#include <ctime>
#include <cctype>
#include <climits>
#include <stdint.h>
#include <bitset>
#include <queue>
#include <deque>
#include <list>

#include <sys/time.h>
#include <pthread.h>


using namespace std;




struct params_t {

    leveldb::DB* db; // the db to store all the results

    ifstream *infile;
    
    uint64_t batch_size;
    uint64_t sequence_counter;
    
    pthread_mutex_t* read_lock;
    pthread_mutex_t* write_lock;

};


inline void ltrim(string& str) {
    string::size_type pos = 0;
    while (pos < str.size() && (isspace(str[pos]))) pos++;
    str.erase(0, pos);
}

inline void rtrim(string& str) {
    string::size_type pos = str.size();
    while (pos > 0 && (isspace(str[pos - 1]))) pos--;
    str.erase(pos);
}

inline void trim2(string& str) {
    ltrim(str);
    rtrim(str);
}




// This class is specific to linux :(

class mu_timer {
private:
    struct timeval start_tv;
public:

    mu_timer() {
        gettimeofday(&start_tv, NULL);
    }

    void reset() {
        gettimeofday(&start_tv, NULL);
    }

    double elapsed_time() {
        struct timeval tv;
        gettimeofday(&tv, NULL);

        return (double) (tv.tv_sec - start_tv.tv_sec) + (tv.tv_usec - start_tv.tv_usec) / 1000000.0;
    }
};

#pragma pack(push)  /* push current alignment to stack */
#pragma pack(4)     /* set alignment to 4 byte boundary */

struct level_key_t {
    uint32_t loc;
    uint64_t key;

    int comp(const level_key_t & y) const {

        if (loc < y.loc) {
            return -1;
        } else if (loc > y.loc) {
            return 1;
        }

        if (key < y.key) {
            return -1;
        } else if (key > y.key) {
            return 1;
        }

        return 0;
    }
};
#pragma pack(pop)   /* restore original alignment from stack */



class chrLocComp : public leveldb::Comparator {
public:
    // Three-way comparison function:
    //   if a < b: negative result
    //   if a > b: positive result
    //   else: zero result

    int Compare(const leveldb::Slice& a, const leveldb::Slice& b) const {
        level_key_t a1, b1;

        a1 = *((level_key_t*) a.data());
        b1 = *((level_key_t*) b.data());

        //cerr << "a1: loc: " << a1.loc << '\n';
        //cerr << "a1: key: " << a1.key << '\n';

        //cerr << "b1: loc: " << b1.loc << '\n';
        //cerr << "b1: key: " << b1.key << '\n';

        //cerr << "Compare: " << a1.comp(b1) << '\n';

        return a1.comp(b1);
    }

    // Ignore the following methods for now:

    const char* Name() const {
        return "chrLocComp";
    }

    void FindShortestSeparator(std::string*, const leveldb::Slice&) const {
    }

    void FindShortSuccessor(std::string*) const {
    }
};




// MT code taken from
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c
// Ported to C++ by John Mu

class MT_random {
private:
    static const uint64_t NN = 312;
    static const uint64_t MM = 156;
    static const uint64_t MATRIX_A = 0xB5026F5AA96619E9ULL;
    static const uint64_t UM = 0xFFFFFFFF80000000ULL; /* Most significant 33 bits */
    static const uint64_t LM = 0x7FFFFFFFULL; /* Least significant 31 bits */




    /* The array for the state vector */
    uint64_t mt[NN];
    /* mti==NN+1 means mt[NN] is not initialized */
    uint64_t mti;

    void init_all() {
        mti = NN + 1;
    }

public:

    MT_random() {

        init_all();
        //init_genrand64(time(NULL));
        init_genrand64(0);
    }

    MT_random(unsigned long long seed) {

        init_all();
        init_genrand64(seed);
    }

    /* initializes mt[NN] with a seed */
    void init_genrand64(unsigned long long seed) {
        mt[0] = seed;
        for (mti = 1; mti < NN; mti++)
            mt[mti] = (6364136223846793005ULL * (mt[mti - 1] ^ (mt[mti - 1] >> 62)) + mti);
    }

    /* initialize by an array with array-length */
    /* init_key is the array for initializing keys */

    /* key_length is its length */
    void init_by_array64(unsigned long long init_key[],
        unsigned long long key_length) {
        unsigned long long i, j, k;
        init_genrand64(19650218ULL);
        i = 1;
        j = 0;
        k = (NN > key_length ? NN : key_length);
        for (; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 3935559000370003845ULL))
                + init_key[j] + j; /* non linear */
            i++;
            j++;
            if (i >= NN) {
                mt[0] = mt[NN - 1];
                i = 1;
            }
            if (j >= key_length) j = 0;
        }
        for (k = NN - 1; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 62)) * 2862933555777941757ULL))
                - i; /* non linear */
            i++;
            if (i >= NN) {
                mt[0] = mt[NN - 1];
                i = 1;
            }
        }

        mt[0] = 1ULL << 63; /* MSB is 1; assuring non-zero initial array */
    }

    /* generates a random number on [0, 2^64-1]-interval */
    unsigned long long genrand64_int64(void) {
        uint64_t i;
        uint64_t x;
        static uint64_t mag01[2] = {0ULL, MATRIX_A};

        if (mti >= NN) { /* generate NN words at one time */

            /* if init_genrand64() has not been called, */
            /* a default initial seed is used     */
            if (mti == NN + 1)
                init_genrand64(5489ULL);

            for (i = 0; i < NN - MM; i++) {
                x = (mt[i] & UM) | (mt[i + 1] & LM);
                mt[i] = mt[i + MM] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
            }
            for (; i < NN - 1; i++) {
                x = (mt[i] & UM) | (mt[i + 1] & LM);
                mt[i] = mt[i + (MM - NN)] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];
            }
            x = (mt[NN - 1] & UM) | (mt[0] & LM);
            mt[NN - 1] = mt[MM - 1] ^ (x >> 1) ^ mag01[(int) (x & 1ULL)];

            mti = 0;
        }

        x = mt[mti++];

        x ^= (x >> 29) & 0x5555555555555555ULL;
        x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
        x ^= (x << 37) & 0xFFF7EEE000000000ULL;
        x ^= (x >> 43);

        return x;
    }

    /* generates a random number on [0, 2^63-1]-interval */
    long long genrand64_int63(void) {
        return (long long) (genrand64_int64() >> 1);
    }

    /* generates a random number on [0,1]-real-interval */
    double genrand64_real1(void) {
        return (genrand64_int64() >> 11) * (1.0 / 9007199254740991.0);
    }

    /* generates a random number on [0,1)-real-interval */
    double genrand64_real2(void) {
        return (genrand64_int64() >> 11) * (1.0 / 9007199254740992.0);
    }

    /* generates a random number on (0,1)-real-interval */
    double genrand64_real3(void) {
        return ((genrand64_int64() >> 12) + 0.5) * (1.0 / 4503599627370496.0);
    }

    // this is not the "correct" way... but good enough

    double genrand_norm(double mean, double sd) {
        return (genrand_norm() * sd) +mean;
    }

    double genrand_norm() {
        double sum = 0;

        for (int i = 0; i < 48; i++) {
            sum = sum + (double) (rand() + 1.0) / (RAND_MAX + 1.0);
        }

        // only for 48
        return (sum - 24) / 2;
    }

    double genrand_exp(double mean) {
        return -mean * log(1 - genrand64_real2());
    }

    int64_t genrand_binom(uint64_t n, double p) {

        if (p < 0 || p > 1) {
            return -1;
        }

        int64_t total = 0;
        for (uint64_t i = 0; i < n; i++) {
            total += (genrand64_real1() < p ? 1 : 0);
        }

        return total;
    }

    // generate random int inclusive of range [start,end]
    uint64_t genrand_int_range(uint64_t start, uint64_t end) {
        return (uint64_t) ((genrand64_real2()*(end - start + 1)) + start);
    }

    // bernouli p is true

    bool genrand_bern(double p) {
        return genrand64_real2() <= p;
    }

};



template <typename T> string toStr(T tmp) {
    ostringstream out;
    out << tmp;
    return out.str();
}

template <typename T> T strTo(string tmp) {
    T output;
    istringstream in(tmp);
    in >> output;
    return output;
}


#endif	/* MAIN_H */

