#!/usr/bin/env python

import sys, os
import re
import numpy

class TestPoint:
    num_threads = 0
    ref_len     = 0
    read_len    = 0

    result_list = [float]*0
    max_result  = 0
    min_result  = 0
    avg_result  = 0
    pruned_list = [float]*0
    pruned_max  = 0
    pruned_min  = 0
    pruned_avg  = 0

    def __init__(self, num_threads, ref_len, read_len):
        self.num_threads = num_threads
        self.ref_len     = ref_len
        self.read_len    = read_len

        self.result_list = [float]*0
        self.max_result  = 0
        self.min_result  = 0
        self.avg_result  = 0
        self.pruned_list = [float]*0
        self.pruned_max  = 0
        self.pruned_min  = 0
        self.pruned_avg  = 0

    def add_result(self, result):
        self.result_list.append(result)
        self.max_result = max(self.result_list)
        self.min_result = min(self.result_list)
        self.avg_result = sum(self.result_list)/float(len(self.result_list))

    def get_results(self):
        return [self.avg_result, self.min_result, self.max_result]

    def get_gcups(self):
        return [r * self.ref_len * self.read_len / 1e+09 for r in [self.avg_result, self.min_result, self.max_result]]

    def prune_results(self):
        std_dev = numpy.std(self.result_list)

        for r in self.result_list:
            if (abs(r - self.avg_result) < std_dev):
                self.pruned_list.append(r)

        self.pruned_max = max(self.pruned_list)
        self.pruned_min = min(self.pruned_list)
        self.pruned_avg = sum(self.pruned_list)/float(len(self.pruned_list))

    def get_pruned_results(self):
        return [self.pruned_avg, self.pruned_min, self.pruned_max]

    def get_pruned_gcups(self):
        return [r * self.ref_len * self.read_len / 1e+09 for r in [self.pruned_avg, self.pruned_min, self.pruned_max]]

def main():
    if len(sys.argv) == 1:
        print 'Need to specify a results file to process.'
        sys.exit(-1)

    results_filename = sys.argv[1]
    print 'Trying to read results file ' + results_filename + '.'

    if not os.path.isfile(results_filename):
        print 'Could not open file ' + results_filename + ' for reading.'
        sys.exit(-1)

    results_file = open(results_filename)

    results_list = [];

    re_header = re.compile("([0-9]+) threads, reflen ([0-9]+) readlen ([0-9]+)")
    re_result = re.compile("= (.*)$")

    while 1:
        line = results_file.readline()
        if line == "":
            break

        if re_header.search(line):
            redata = re_header.search(line).groups()
            result = TestPoint(int(redata[0]), int(redata[1]), int(redata[2]))
            results_list.append(result)
            last_result = len(results_list)-1
        elif re_result.search(line):
            redata = re_result.search(line).groups()
            results_list[last_result].add_result(float(redata[0]))

    print 'threads,ref_len,read_len,gcups,gcups_no-outliers'
    for r in results_list:
        r.prune_results()
        print str(r.num_threads) + ',' + str(r.ref_len) + ',' + str(r.read_len) + ',' + str(r.get_gcups()[0]) + ',' + str(r.get_pruned_gcups()[0])

    results_file.close()

if __name__ == '__main__':
    main()

