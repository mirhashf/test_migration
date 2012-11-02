#!/usr/bin/env python

import sys, os
import re

class TestPoint:
    num_threads = 0
    ref_len     = 0
    read_len    = 0
    result_list = [float]*0
    max_result  = 0
    min_result  = 0
    avg_result  = 0

    def __init__(self, num_threads, ref_len, read_len):
        self.num_threads = num_threads
        self.ref_len     = ref_len
        self.read_len    = read_len
        self.result_list = [float]*0
        self.max_result  = 0
        self.min_result  = 0
        self.avg_result  = 0

    def add_result(self, result):
        self.result_list.append(result);

    def results(self):
        return [min(self.result_list), max(self.result_list), sum(self.result_list)/float(len(self.result_list))]

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

    for r in results_list:
        print str(r.num_threads) + ' threads ref_len ' + str(r.ref_len) + ' read_len ' + str(r.read_len) + ' ' + str(r.results())

    results_file.close()

if __name__ == '__main__':
    main()

