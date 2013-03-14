#!/usr/bin/python

import sys, os
import re
import time
from datetime import datetime

# re.search('(\w{3}/[0-9]{2}/[0-9]{4}\W+[0-9]{2}:[0-9]{2}:[0-9]{2}\W+\w{2})', string)
# datetime.strptime(string, '%b/%d/%Y %I:%M:%S %p')

time_pattern = re.compile('(\w{3}/[0-9]{2}/[0-9]{4}\W+[0-9]{2}:[0-9]{2}:[0-9]{2}\W+\w{2})');

align_times = [];

while 1:
    line = sys.stdin.readline()
    if line == '':
        break

    if line.find('Getting status of job') != -1:
        jobid = re.search('Getting status of job (.*)...', line).group(1)
        print 'For job ' + jobid

        while 1:

            if line == '':
                break

            if line.find('Bina Aligner index') != -1:
                start_time_str = time_pattern.search(line).group(0)
                start_time = datetime.strptime(start_time_str, '%b/%d/%Y %I:%M:%S %p')
                align_times.append(start_time);

            line = sys.stdin.readline()

            if line.find('Completed alignment') != -1:
                align_time_str = time_pattern.search(line).group(0)
                align_time = datetime.strptime(align_time_str, '%b/%d/%Y %I:%M:%S %p')
                align_time_elapsed = align_time - align_times[-1]
                align_times.append(align_time);
                # print 'alignment done ' + align_time_str
                print 'alignment took %d seconds (%s)' % (align_time_elapsed.seconds, str(align_time_elapsed))

            if line.find('Indel realignment complete.') != -1:
                realign_time_str = time_pattern.search(line).group(0)
                realign_time = datetime.strptime(realign_time_str, '%b/%d/%Y %I:%M:%S %p')
                realign_time_elapsed = realign_time - align_time
                # print 'realignment done ' + realign_time_str
                print 'realignment took %d seconds (%s)' % (realign_time_elapsed.seconds, str(realign_time_elapsed))

            if line.find('Genotyping complete.') != -1:
                genotype_time_str = time_pattern.search(line).group(0)
                genotype_time = datetime.strptime(genotype_time_str, '%b/%d/%Y %I:%M:%S %p')
                genotype_time_elapsed = genotype_time - realign_time
                # print 'genotyping done ' + genotype_time_str
                print 'genotyping took %d seconds (%s)' % (genotype_time_elapsed.seconds, str(genotype_time_elapsed))

            if line.find('third-party structural variation detectors') != -1:
                svtools_time_str = time_pattern.search(line).group(0)
                svtools_time = datetime.strptime(svtools_time_str, '%b/%d/%Y %I:%M:%S %p')
                svtools_time_elapsed = svtools_time - genotype_time
                # print 'sv done ' + svtools_time_str
                print 'svtools took %d seconds (%s)' % (svtools_time_elapsed.seconds, str(svtools_time_elapsed))


            if line.find('Merged variant calling results') != -1:
                merge_time_str = time_pattern.search(line).group(0)
                merge_time = datetime.strptime(merge_time_str, '%b/%d/%Y %I:%M:%S %p')
                merge_time_elapsed = merge_time - svtools_time
                # print 'merging done ' + merge_time_str
                print 'merging took %d seconds (%s)' % (merge_time_elapsed.seconds, str(merge_time_elapsed))

            if line.find('All variant quality score recalibrations completed') != -1:
                vqsr_time_str = time_pattern.search(line).group(0)
                vqsr_time = datetime.strptime(vqsr_time_str, '%b/%d/%Y %I:%M:%S %p')
                vqsr_time_elapsed = vqsr_time - merge_time
                # print 'vqsr done ' + vqsr_time_str
                print 'vqsr took %d seconds (%s)' % (vqsr_time_elapsed.seconds, str(vqsr_time_elapsed))

            if line.find('Bina Structural Variation detector results written') != -1:
                bsv_time_str = time_pattern.search(line).group(0)
                bsv_time = datetime.strptime(bsv_time_str, '%b/%d/%Y %I:%M:%S %p')
                bsv_time_elapsed = bsv_time - vqsr_time
                # print 'bina sv done ' + bsv_time_str
                print 'bina sv took %d seconds (%s)' % (bsv_time_elapsed.seconds, str(bsv_time_elapsed))

            if line.find('Job dequeued.') != -1:
                end_time_str = time_pattern.search(line).group(0)
                end_time = datetime.strptime(end_time_str, '%b/%d/%Y %I:%M:%S %p')
                end_time_elapsed = end_time - start_time
                # print 'end time ' + end_time_str
                print 'total took %d seconds (%s)' % (end_time_elapsed.seconds, str(end_time_elapsed))
                break

        print ''


