/**
 * testUtils.h
 * Declares several utilities common to test functionality.
 *
 * Written by AJ Minich, Jan 2012
 */

#ifndef TESTUTILS_H

#define TESTUTILS_H

#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdint.h>

#include "../src/general_utils.h"
#include "../src/constants.h"
#include "../src/pen_tester.h"
#include "../src/types.h"
#include "../src/align_store.h"

using namespace std;

seqalto_params_t * getDefaultSeqParams();


#endif