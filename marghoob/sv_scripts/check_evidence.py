#!/usr/bin/env python

# This script check for potential evidence for paired-end and split-read methods given BAMs, golden data and events output

import os
import sys
import subprocess
import argparse
import logging

def get_coverage(coverage_file):
  coverage_file_fd = open(coverage_file, "r")


def valid_interval(a):
  return a[0] < a[1]

def interval_overlap(a, b):
  a_s, a_e = a
  b_s, b_e = b

  o_s, o_e = (max(a_s, b_s), min(a_e, b_e))

  return o_s, o_e

def get_overlapping_intervals(interval, intervals):
  overlapping_intervals = []
  for candidate in intervals:
    overlap = interval_overlap(interval, candidate)
    if valid_interval(overlap) and (overlap[1] - overlap[0]) >= 0.7 * (interval[1] - interval[0]) and (candidate[1] - candidate[0]) >= 0.7 * (interval[1] - interval[0]):
      overlapping_intervals.append(candidate)
  return overlapping_intervals

def get_overlap_events(events):
  overlap_events = []
  intervals = set()
  for event1 in events:
    overlap_found = False
    for event2 in events:
      if event1 == event2:
        continue
      if event1[1] < event2[1]:
        continue
      o_s, o_e = interval_overlap(event1[3], event2[3])
      if valid_interval((o_s, o_e)) and (o_s, o_e) not in intervals:
        intervals.add((o_s, o_e))
        if o_e - o_s > 100:
          overlap_events.append((event1[0], "overlap_" + event1[1] + "_" + event2[1], event1[2], (o_s, o_e)))
        overlap_found = True
    if not overlap_found:
      overlap_events.append(event1)
  return overlap_events
  

def load_deletion_events(bedfile):

  all_events = []
  for event in bedfile.readlines():
    fields = event.split()
    #(sample, tech, chromosome, (start, end)) = (fields[0], fields[1], fields[2], (int(fields[3]), int(fields[4])))
    (sample, tech, chromosome, (start, end)) = ("NA12878", "exp", fields[0], (int(fields[1]), int(fields[2])))

    if (options.chromosome and chromosome != options.chromosome) or sample != "NA12878": continue
    if (end - start) < options.min or (end - start) > options.max: continue
    all_events.append((sample, tech, chromosome, (start, end)))
    #print all_events[-1]

  if all_events:
    all_events.sort(key=lambda tup: (tup[3], tup[1]))
  bedfile.close()

  #print "%d deletion events loaded from %s" %(len(all_events), bedfile.name)
  return all_events

def load_gff_intervals(gff):
    all_events = []

    event_dict = {}
    event_dict["deletion"] = "Deletion"
    if gff:
        for line in gff.readlines():
            fields = line.split()
            chromosome, tool, event_type, start, end = fields[0:5]
            if event_type != event_dict[options.event_type] or chromosome != options.chromosome: continue
            start = int(start)
            end = int(end)
            if (end - start) < options.min or (end - start) > options.max: continue
            all_events.append((start, end))
        all_events.sort()
        gff.close()

        #print "Loaded %d events from %s" %(len(all_events), gff.name)
    return all_events

def load_cnvnator_intervals(cnvnator):
    all_events = []

    return all_events

def load_breakdancer_intervals(breakdancer):
    all_events = []

    event_dict = {}
    event_dict["deletion"] = "DEL"

    if breakdancer:
        for line in breakdancer.readlines():
            if line[0] == "#": continue
            fields = line.split('\t')
            score = int(fields[7])
            #if score < 10: continue
            chromosome, start, end, event_type = fields[0], fields[1], fields[4], fields[6]
            start, end = int(start), int(end)
            if (end - start) < options.min or (end - start) > options.max or chromosome != options.chromosome or event_type != event_dict[options.event_type]: continue
            all_events.append((start, end))
        all_events.sort()
        breakdancer.close()
        print "Loaded %d events from %s" %(len(all_events), breakdancer.name)
    return all_events

def load_lumpy_intervals(lumpy):
    all_events = []

    event_dict = {}
    event_dict["deletion"] = "TYPE:DELETION"

    if lumpy:
        for line in lumpy.readlines():
            fields = line.split('\t')
            chromosome, start, end, event_type = fields[0], fields[1], fields[5], fields[10]
            start, end = int(start), int(end)
            if (end - start) < options.min or (end - start) > options.max or chromosome != options.chromosome or event_type != event_dict[options.event_type]: continue
            all_events.append((start, end))
        all_events.sort()
        lumpy.close()
        print "Loaded %d events from %s" %(len(all_events), lumpy.name)
    return all_events

def load_intervals(intervals_file):
    all_events = []

    if intervals_file:
        for line in intervals_file.readlines():
            fields = line.split()
            start, end = int(fields[0]), int(fields[1])
            all_events.append((start, end))
        all_events.sort()
        intervals_file.close()
    return all_events    

def load_dgv_events(dgv):
    all_events = []

    if dgv:
        for line in dgv.readlines():
            fields = line.split()

            chromosome, start, end, size, gain, loss, year = fields

            if chromosome != options.chromosome or int(loss) == 0:
                continue

            #all_events.append(("dgv", "dgv", "chr1", (int(start), int(end))))
        dgv.close()
    return all_events

def get_overlapping_intervals(a, interval_list):
  overlap_list = []
  for interval in interval_list:
    if valid_interval(interval_overlap(a, interval)):
      overlap_list.append(interval)
  return overlap_list

parser = argparse.ArgumentParser(description='Process BAM to check for SV support')
parser.add_argument("-b", "--bed", help="Bedfile for deletion events", metavar="bed", required=True, type=file)
parser.add_argument("-f", "--bam", help="BAM file for analysis", metavar="bam", required=True)
parser.add_argument("-s", "--samtools", help="Samtools binary", metavar="Samtools binary", required=False, default="/home/marghoob/seqalto/third-party/native/bin/samtools")
parser.add_argument("--max", help="Max event size", metavar="MaxEventSize", required=False, type=int, default=10000000)
parser.add_argument("--min", help="Min event size", metavar="MinEventSize", required=False, type=int, default=500)
parser.add_argument("-w", "--window", type=int, help="Window size", metavar="Window size", default=1000)
parser.add_argument("--template_mean", help="Template length mean", type=int, metavar="template mean", required=False, default=320)
parser.add_argument("--template_sd", help="Template length s.d.", type=int, metavar="template s.d.", required=False, default=70)
parser.add_argument("--read_length", help="Read length", type=int, metavar="read length", required=False, default=100)
parser.add_argument("--coverage_file", help="Coverage file", metavar="FILE")
parser.add_argument("--min_overlap", type=float, help="Minimum percent overlap with the specified events", metavar="FLOAT", default=50.0)
parser.add_argument("--breakdancer_gff", help="Breakdancer GFF", metavar="BreakdancerGFF", type=file)
parser.add_argument("--breakdancer", help="Breakdancer Native", metavar="Breakdancer", type=file)
parser.add_argument("--cnvnator_gff", help="Cnvnator GFF", type=file)
parser.add_argument("--cnvnator", help="Cnvnator native", type=file)
parser.add_argument("--lumpy", help="Lumpy native output", metavar="Lumpy", type=file)
parser.add_argument("--pindel_gff", help="Pindel GFF", metavar="Pindel", type=file)
parser.add_argument("--pindel_intervals", help="Pindel intervals", type=file)
parser.add_argument("--dgv", help="DGV event file", metavar="DGV", type=file)
parser.add_argument("--chromosome", help="Chromosome", metavar="Chromosome")
parser.add_argument("--event_type", help="Event type", metavar="EventType", default="deletion")


options = parser.parse_args()

samtools = options.samtools

logfile = open("logs/check_evidence.log", "w")

total_events = 0
no_support_events = 0

techs = set()

breakdancer_intervals = load_gff_intervals(options.breakdancer_gff) + load_breakdancer_intervals(options.breakdancer)
pindel_intervals = load_gff_intervals(options.pindel_gff) + load_intervals(options.pindel_intervals)
lumpy_intervals = load_lumpy_intervals(options.lumpy)
cnvnator_intervals = load_gff_intervals(options.cnvnator_gff) + load_cnvnator_intervals(options.cnvnator)

logger = logging.getLogger(__name__)

def get_paired_end_support(sample, tech, interval, chromosome, breakdancer_intervals):
  start, end = interval
  event_size = end - start

  window_size = int(max(options.window, event_size * 0.1))
  left_start = max(1, start - window_size - options.template_mean - 3 * options.template_sd)
  left_end = start + window_size - options.template_mean + 3 * options.template_sd

  right_start = max(1, end - window_size + options.read_length - options.template_mean - 3 * options.template_sd)
  right_end = end + window_size + options.read_length - options.template_mean + 3 * options.template_sd

  support_file_left = "support/%d_%d.left.paired.sam" %(start, end)
  support_file_right = "support/%d_%d.right.paired.sam" %(start, end)
  support_fp_left = open(support_file_left, "w")
  support_fp_right = open(support_file_right, "w")
  
  samtools_cmd_left = [options.samtools, "view", "-F", "28", "-f", "32", "-q", "20", options.bam, "%s:%d-%d" %(chromosome, left_start, left_end)]
  awk_cmd_left = ["/usr/bin/awk", "-F", "\t", "{if (($4 + $9 > %d) && ($4 + $9 < %d) && ($9 > %d)) print $0}" %(right_start, right_end, (options.template_mean + event_size) * 0.8)]

  samtools_cmd_right = [options.samtools, "view", "-F", "44", "-f", "16", "-q", "20", options.bam, "%s:%d-%d" %(chromosome, right_start , right_end)]
  awk_cmd_right = ["/usr/bin/awk", "-F", "\t", "{if (($4 + %d + $9 > %d) && ($4 + %d + $9 < %d) && (-$9 > %d)) print $0}" %(options.read_length, left_start, options.read_length, left_end, (options.template_mean + event_size) * 0.8)]

  samtools_p_left = subprocess.Popen(samtools_cmd_left, stdout=subprocess.PIPE)
  samtools_p_right = subprocess.Popen(samtools_cmd_right, stdout=subprocess.PIPE)

  awk_p_left = subprocess.Popen(awk_cmd_left, stdin=samtools_p_left.stdout, stdout=support_fp_left)
  awk_p_right = subprocess.Popen(awk_cmd_right, stdin=samtools_p_right.stdout, stdout=support_fp_right)

  awk_ret_left = awk_p_left.wait()
  awk_ret_right = awk_p_right.wait()

  samtools_p_left.communicate()
  samtools_p_right.communicate()

  samtools_ret_left = samtools_p_left.wait()
  samtools_ret_right = samtools_p_right.wait()

  support_fp_left.close()
  support_fp_right.close()

  if samtools_ret_left or awk_ret_left or samtools_ret_right or awk_ret_right:
    sys.exit(1)

  support_fp_left = open(support_file_left, "r")
  support_fp_right = open(support_file_right, "r")

  support_size = len(support_fp_left.readlines()) + len(support_fp_right.readlines())
  overlapping_intervals = get_overlapping_intervals(interval, breakdancer_intervals)
  breakdancer_support_size = len(overlapping_intervals)

  #if (overlapping_intervals and not support_size) or (not overlapping_intervals and support_size):
  #  print tech, left_end, right_end, right_end - left_end, "paired-end support size = %d" %(support_size), get_overlapping_intervals(interval, breakdancer_events)
  support_fp_left.close()
  support_fp_right.close()

  return (support_size if support_size >= 2 else 0, breakdancer_support_size, overlapping_intervals)

def get_split_read_support(sample, tech, interval, pindel_intervals):
  start, end = interval
  window_size = int(max(options.window, (end - start) * 0.1))
  chromosome = options.chromosome

  left_start = max(1, start - window_size - options.template_mean - 3 * options.template_sd)
  left_end = start + window_size - options.template_mean + 3 * options.template_sd
  samtools_cmd_left = [options.samtools, "view", "-f", "8", "-F", "16", "-q", "20", options.bam, chromosome + ":" + str(left_start) + "-" + str(left_end)]

  right_start = max(1, end - window_size + options.read_length - options.template_mean - 3 * options.template_sd)
  right_end = end + window_size + options.read_length - options.template_mean + 3 * options.template_sd
  samtools_cmd_right = [options.samtools, "view", "-f", "24", "-q", "20", options.bam, chromosome + ":" + str(right_start) + "-" + str(right_end)]

  left_support_file = "support/" + str(start) + "_" + str(end) + ".left.sam"
  right_support_file = "support/" + str(start) + "_" + str(end) + ".right.sam"

  cmd_left_stdout = open(left_support_file, "w")
  cmd_right_stdout = open(right_support_file, "w")

  ret0 = subprocess.call(samtools_cmd_left, stdout=cmd_left_stdout, stderr=logfile)
  ret1 = subprocess.call(samtools_cmd_right, stdout=cmd_right_stdout, stderr=logfile)

  if ret0:
    print "Execution of", samtools_cmd_left, "failed with return code %d" %(ret0)
    sys.exit(1)

  if ret1:
    print "Execution of", samtools_cmd_right, "failed with return code %d" %(ret1)
    sys.exit(1)

  cmd_left_stdout.close()
  cmd_right_stdout.close()

  # now analyze support files
  left_support_fp = open(left_support_file, "r")
  left_differences = []
  for sam_record in left_support_fp.readlines():
    sam_fields = sam_record.split('\t')
    if int(sam_fields[4]) < 10:
      continue
    alignment_location = int(sam_fields[3])
    if alignment_location >= end: continue
    difference = start - alignment_location
    left_differences.append(difference)
  left_support_fp.close()

  right_support_fp = open(right_support_file, "r")
  right_differences = []
  for sam_record in right_support_fp.readlines():
    sam_fields = sam_record.split('\t')
    if int(sam_fields[4]) < 10:
      continue
    alignment_location = int(sam_fields[3])
    if alignment_location <= start: continue
    difference = alignment_location - end
    right_differences.append(difference)
  right_support_fp.close()

  support_size = len(left_differences) + len(right_differences)

  overlapping_intervals = get_overlapping_intervals(interval, pindel_intervals)
  pindel_support_size = len(overlapping_intervals)
  #print tech, left_end, right_end, right_end - left_end, "split-read support size = %d" %(support_size), left_differences, right_differences
  return support_size if support_size >= 2 else 0, pindel_support_size

def doc_support(sample, tech, interval, chromosome, doc_intervals):
    overlapping_intervals = get_overlapping_intervals(interval, doc_intervals)
    return overlapping_intervals

print "tech start end size paired_support_size tool_support_size"
for sample, tech, chromosome, (start, end) in load_deletion_events(options.bed) + load_dgv_events(options.dgv):

    if breakdancer_intervals:
        paired_support_size, breakdancer_support_size, breakdancer_overlapping_intervals = get_paired_end_support(sample, tech, (start, end), options.chromosome, breakdancer_intervals)
    #split_support_size, pindel_support_size = get_split_read_support(sample, tech, (start, end), pindel_intervals)

        if paired_support_size and not breakdancer_support_size:
            print "paired-breakdancer", tech, start, end, end - start, paired_support_size, breakdancer_support_size, breakdancer_overlapping_intervals

    if lumpy_intervals:
        paired_support_size, lumpy_support_size, lumpy_overlapping_intervals = get_paired_end_support(sample, tech, (start, end), options.chromosome, lumpy_intervals)
        if paired_support_size and not lumpy_support_size:
            print "paired-lumpy", tech, start, end, end - start, paired_support_size, lumpy_support_size, lumpy_overlapping_intervals

    if cnvnator_intervals:
        cnvnator_overlapping_intervals = doc_support(sample, tech, (start, end), options.chromosome, cnvnator_intervals)
        if not cnvnator_overlapping_intervals:
            print "doc-cnvnator", tech, start, end, end - start, len(cnvnator_overlapping_intervals), cnvnator_overlapping_intervals

    if pindel_intervals:
       split_support_size, pindel_support_size = get_split_read_support(sample, tech, (start, end), pindel_intervals)
       if split_support_size and not pindel_support_size:
           print "split-pindel", tech, start, end, end - start, split_support_size, pindel_support_size
