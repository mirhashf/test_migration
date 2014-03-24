#!/usr/bin/env python2.7

import sys
import os
import argparse
import subprocess
import pysam

parser = argparse.ArgumentParser("Convert genotyped BreakDancer output to VCF", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--sv_file", metavar="sv_file", help="SV file", required=False, default="-")
parser.add_argument("--reference", metavar="reference", help="Reference file", required=True)
parser.add_argument("--sample", metavar="Sample", help="Sample name", required=True)
parser.add_argument("--sort", action="store_true", help="Sort the input")

args = parser.parse_args()

input_handle = sys.stdin if args.sv_file == "-" else open(args.sv_file)

fasta_handle = pysam.Fastafile(args.reference)

def get_contigs(fai_filename):
  fai_file = open(fai_filename)
  contigs = {}
  contigs_order = {}
  linenum = 0
  for line in fai_file.readlines():
    line = line.strip()
    line_items = line.split("\t")
    name, length = line_items[0:2]
    name = name.split(" ")[0]
    contigs[name] = int(length)
    contigs_order[name] = linenum
    linenum += 1
  fai_file.close()
  return contigs, contigs_order

def line_to_tuple(line):
  line = line.strip()
  fields = line.split("\t")
  return tuple(fields[0:2]) + tuple([int(i) for i in fields[2:4]]) + tuple(fields[4:5]) + tuple([int(fields[5])])

contigs, contigs_order = get_contigs(args.reference + ".fai")
contig_names = contigs.keys()
contig_names.sort(key = lambda tup: contigs_order[tup])

contig_str = ""
for contig_name in contig_names:
  contig_str += "##contig=<ID=%s,length=%d>\n" % (contig_name, contigs[contig_name])

sys.stdout.write("""##fileformat=VCFv4.1
##reference=%s
##source=%s
##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">
##INFO=<ID=SVCNV,Number=0,Type=Flag,Description=\"Structural variation or copy number variation\">
##INFO=<ID=MEINFO,Number=4,Type=String,Description=\"Mobile element info of the form NAME,START,END,POLARITY\">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=NORMAL_COUNT,Number=1,Type=Integer,Description=\"Number of normal reads supporting reference\">
##INFO=<ID=NUM_READS,Number=1,Type=Integer,Description=\"Number of reads supporting event\">
##INFO=<ID=SCORE,Number=1,Type=Integer,Description=\"Score reported by tool\">
##INFO=<ID=SVTOOL,Number=1,Type=String,Description=\"Tool used to generate the SV\">
##INFO=<ID=SOURCES,Number=.,Type=String,Description=\"List of original raw SV calls as Toolname:Start:End:Size\">
##INFO=<ID=NUM_SVMETHODS,Number=1,Type=Integer,Description=\"Number of methods supporting the event\">
##INFO=<ID=VT,Number=1,Type=String,Description=\"indicates what type of variant the line represents\">
##INFO=<ID=SVMETHOD,Number=.,Type=String,Description=\"Type of approach used to detect SV: RP (read pair), RD (read depth), SR (split read), JT (junction) or AS (assembly)\">
##INFO=<ID=natorRD,Number=1,Type=Float,Description=\"CNVnator: Normalized RD\">
##INFO=<ID=natorP1,Number=1,Type=Float,Description=\"CNVnator: e-val by t-test\">
##INFO=<ID=natorP2,Number=1,Type=Float,Description=\"CNVnator: e-val by Gaussian tail\">
##INFO=<ID=natorP3,Number=1,Type=Float,Description=\"CNVnator: e-val by t-test (middle)\">
##INFO=<ID=natorP4,Number=1,Type=Float,Description=\"CNVnator: e-val by Gaussian tail (middle)\">
##INFO=<ID=natorQ0,Number=1,Type=Float,Description=\"CNVnator: Fraction of reads with 0 mapping quality\">
##FILTER=<ID=LowQual,Description=\"Low Quality\">
%s##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">
##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">
##ALT=<ID=DUP,Description=\"Duplication\">
##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">
##ALT=<ID=INS,Description=\"Insertion of novel sequence\">
##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">
##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">
##ALT=<ID=INV,Description=\"Inversion\">
##ALT=<ID=CNV,Description=\"Copy number variable region\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n""" % (args.reference, sys.argv[0], contig_str, args.sample))

records = map(line_to_tuple, input_handle.readlines())
#if args.sort:
#  records.sort(key = lambda tup: (contigs_order[tup[1]], tup[2], tup[3], tup[4]))


def process_merged_record(record):
  sv_type, chr1, pos1, pos2, sources, method_count = record

  sources_split = sources.split(";")
  tool_name_to_method = {"BreakSeq": "JT", "Breakdancer": "RP", "CNVnator": "RD", "Pindel": "SR", "HaplotypeCaller": "AS"}
  lists = {}
  lists["BreakSeq"] = []
  lists["Breakdancer"] = []
  lists["CNVnator"] = []
  lists["Pindel"] = []
  lists["HaplotypeCaller"] = []

  fix = 0
  if sv_type == "INS": fix = -1 
  for sources_split_fields in sources_split:
    if "Pindel" in sources_split_fields or "Breakdancer" in sources_split_fields:
      source_tool, source_start, source_end, source_size, normal_count, num_reads, source_gt = sources_split_fields.split("_")
      source_end = int(source_end) + fix
      lists[source_tool].append((int(source_start), int(source_end), abs(int(source_size)), source_gt, int(normal_count), int(num_reads)))
    elif "CNVnator" in sources_split_fields:
      source_tool, source_start, source_end, source_size, rd, p1, p2, p3, p4, q0, source_gt = sources_split_fields.split("_")
      lists[source_tool].append((int(source_start), int(source_end) + 1, int(source_size), source_gt, float(rd), float(p1), float(p2), float(p3), float(p4), float(q0)))
    elif "BreakSeq" in sources_split_fields:
      source_tool, source_start, source_end, source_size, filter_str, source_gt = sources_split_fields.split("_")
      source_end = int(source_end) + fix
      lists[source_tool].append((int(source_start), int(source_end), int(source_size), source_gt, filter_str))
    elif "HaplotypeCaller" in sources_split_fields:
      source_tool, source_start, source_end, source_size, source_gt = sources_split_fields.split("_")
      source_end = int(source_end) + fix
      lists[source_tool].append((int(source_start), int(source_end), int(source_size), source_gt))
    else: return None

  done = False
  size = None

  bd_normal_count = None
  bd_num_reads = None
  pindel_normal_count = None
  pindel_num_reads = None

  is_pass = False
  is_precise = False
  if sv_type in ["DEL", "INV", "INS", "DUP", "DUP:TANDEM"]:
    if len(lists["Pindel"]) == 1:
      source_start, source_end, source_size, source_gt, pindel_normal_count, pindel_num_reads = lists["Pindel"][0]
      if float(source_size) / (float(pos2) - float(pos1)) >= 0.5:
        pos1, pos2, gt, size = source_start, source_end, source_gt, source_size
        is_precise = True
        is_pass = True
      if sv_type == "INS":
        pos1, pos2, gt, size = source_start, source_end, source_gt, source_size
        is_precise = True
        is_pass = True
    elif len(lists["BreakSeq"]) == 1:
      source_start, source_end, source_size, source_gt = lists["BreakSeq"][0][0:4]
      if float(source_size) / (float(pos2) - float(pos1)) >= 0.5:
        pos1, pos2, gt, size = source_start, source_end, source_gt, source_size
        is_precise = True
        is_pass = True
      if sv_type == "INS":
        pos1, pos2, gt = source_start, source_end, source_gt
        is_precise = True
        size = None
        if lists["HaplotypeCaller"] == 1:
          size = lists["HaplotypeCaller"][0][2]
        is_pass = True
    elif len(lists["Breakdancer"]) == 1:
      source_start, source_end, source_size, source_gt, bd_normal_count, bd_num_reads = lists["Breakdancer"][0]
      if float(source_size) / (float(pos2) - float(pos1)) >= 0.5:
        pos1, pos2, gt, size = source_start, source_end, source_gt, source_size
        is_pass = True
    elif len(lists["CNVnator"]) == 1:
      source_start, source_end, source_size, source_gt = lists["CNVnator"][0][0:4]
      if float(source_size) / (float(pos2) - float(pos1)) >= 0.5:
        pos1, pos2, gt, size = source_start, source_end, source_gt, source_size
        is_pass = True
    if not is_pass:
      min_size = 10000
      for tool in ["Pindel", "BreakSeq", "Breakdancer", "CNVnator"]:
        if lists[tool] and len(lists[tool]) < min_size:
          size = lists[tool][0][2]
          gt = lists[tool][0][3]
          min_size = len(lists[tool])
      if sv_type in ["DEL", "INV"]:
        size = pos2 - pos1 - 1

  else: return None

  #if size is None: return None

  sources_list = []
  methods_list = []
  for tool in ["Breakdancer", "BreakSeq", "Pindel", "CNVnator", "HaplotypeCaller"]:
    if lists[tool]:
      sources_list += ["%s:%d:%d:%d" % (tool, item[0], item[1], item[2]) for item in lists[tool]]
      methods_list.append(tool)

  if sv_type == "INS":
    if lists["Pindel"] or lists["BreakSeq"]:
      pos1 = min([item[0] for item in lists["Pindel"] + lists["BreakSeq"]])
      pos2 = max([item[1] for item in lists["Pindel"] + lists["BreakSeq"]])

  if methods_list == ["HaplotypeCaller"]: return None

  sources_string = "" if len(sources_list) == 1 else (";SOURCES=" + ",".join(sources_list))

  tool = "MetaMergeSV"
  is_pass = is_pass and (len(methods_list) > 1)
  if len(sources_list) == 1: tool = sources_list[0].split(":")[0]
  if sv_type != "DEL":
    #if len(sources_list) == 1 and tool == "Breakdancer": is_pass = True
    if len(sources_list) == 1 and tool == "BreakSeq": is_pass = lists["BreakSeq"][0][4] == "PASS"
    if len(sources_list) == 1 and tool == "CNVnator" and sv_type == "DUP": is_pass = lists["CNVnator"][0][7] <= 1.0
    if sv_type == "INS" and len(lists["HaplotypeCaller"]) == 1: is_pass = True

  is_low_qual = False
  if len(sources_list) == 1 and tool == "BreakSeq": is_low_qual = lists["BreakSeq"][0][4] == "LowQual"

  methods_str = ",".join([tool_name_to_method[toolname] for toolname in methods_list])
  info = "VT=SV;SVTOOL=%s%s;NUM_SVMETHODS=%d;SVMETHOD=%s" % (tool, sources_string, len(methods_list), methods_str)
  if pindel_normal_count: info += ";NORMAL_COUNT=%d;NUM_READS=%d" % (pindel_normal_count, pindel_num_reads)
  if bd_normal_count: info += ";NORMAL_COUNT=%d;NUM_READS=%d" % (bd_normal_count, bd_num_reads)
  if not is_precise: info += ";IMPRECISE"

  pos1 = max(1, pos1)
  pos2 = max(pos1, pos2)
  
  return sv_type, chr1, pos1, pos2, size, gt, info, is_pass, is_low_qual 

processed_records = []
for record in records:
  processed_record = process_merged_record(record)
  if processed_record: processed_records.append(processed_record)

if args.sort:
  processed_records.sort(key = lambda tup: (contigs_order[tup[1]], tup[2], tup[3]))

for processed_record in processed_records:
  sv_type, chr1, pos1, pos2, size, gt, info, is_pass, is_low_qual = processed_record
    
  alt_allele = "<%s>" % (sv_type)

  ref_allele = fasta_handle.fetch(chr1, pos1-1, pos1)

  if sv_type in ["DEL"]: size = -size
  info = "END=%d;%s;SVTYPE=%s" % (pos2, info, sv_type)
  if size is not None: info += ";SVLEN=%d" % (size)
  #info = "TOOLNAME=%s;SVLEN=%d;SVTYPE=%s;END=%d;IMPRECISE;NORMAL_COUNT=%d;NUM_READS=%d;SCORE=%g" % (tool, size, sv_type, pos2, normal_read_count, num_reads, score)

  pass_str = "PASS" if is_pass else ("LowQual") # if is_low_qual else ".")
  sys.stdout.write("%s\t%d\t.\t%s\t%s\t.\t%s\t%s\tGT\t%s\n" % (chr1, pos1, ref_allele, alt_allele, pass_str, info, gt))

fasta_handle.close()
