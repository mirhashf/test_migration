#!/usr/bin/python

import json
import sys
import re
import glob
import os

if len(sys.argv) != 2:
  sys.stderr.write("./script <jobdir>\n")
  sys.exit(1)

jobdir = sys.argv[1]

total_reads = 0
total_unmapped = 0
total_uniquely_mapped = 0
total_multiply_mapped = 0
total_duplicate = 0
for alignment_metrics in glob.glob(os.path.join(jobdir, "metrics/*.alignment_metrics.json")):
  alignment_json = json.load(open(alignment_metrics, "r"))
  for key1 in alignment_json:
    for key2 in alignment_json[key1]:
      sub_dict = alignment_json[key1][key2]
      if "read_num" in sub_dict:
        total_reads += sub_dict["read_num"]
        if key1 == "*":
          total_unmapped += sub_dict["read_num"]
      if "multi_mapped_num" in sub_dict:
        total_multiply_mapped += sub_dict["multi_mapped_num"]
      if "duplicate_num" in sub_dict:
        total_duplicate += sub_dict["duplicate_num"]

total_uniquely_mapped = total_reads - total_unmapped - total_multiply_mapped
percent_uniquely_mapped = float(total_uniquely_mapped) * 100.0 / float(total_reads)
percent_unmapped = float(total_unmapped) * 100.0 / float(total_reads)
percent_multi_mapped = 100.0 - percent_uniquely_mapped - percent_unmapped

vcf_stats = {}
for vcf_metrics in glob.glob(os.path.join(jobdir, "metrics/*.vcf_metrics.json")):
  vcf_json = json.load(open(vcf_metrics, "r"))
  for chrom in vcf_json:
    for stat_type in vcf_json[chrom]:
      if stat_type not in vcf_stats:
        vcf_stats[stat_type] = {}
      if not vcf_json[chrom][stat_type]: continue
      for metric in vcf_json[chrom][stat_type]:
        if metric not in vcf_stats[stat_type]:
          vcf_stats[stat_type][metric] = 0
        vcf_stats[stat_type][metric] += float(vcf_json[chrom][stat_type][metric])

snp_ti_tv = vcf_stats["snpStats"]["tiCount"] / vcf_stats["snpStats"]["tvCount"]
snp_het_hom = vcf_stats["snpStats"]["hetCount"] / vcf_stats["snpStats"]["homCount"]
snp_indel = vcf_stats["snpStats"]["totalCount"] / vcf_stats["indelStats"]["totalCount"]
indel_het_hom = vcf_stats["indelStats"]["hetCount"] / vcf_stats["indelStats"]["homCount"]

# now extract fastqc stuff
fastq_metrics = {}
for fastqc_dir in glob.glob(os.path.join(jobdir, "fastqc/*/*")):
  if not os.path.isdir(fastqc_dir): continue
  fastq_metrics[fastqc_dir] = {"quality": 0.0}
  fastqc_file = open(os.path.join(fastqc_dir, "fastqc_data.txt"))

  state = "normal"
  for line in fastqc_file.readlines():
    line_t = line.strip()
    
    if line_t.find("%GC") == 0:
      fastq_metrics[fastqc_dir]["gc_content"] = float(line_t.split()[1])
    elif line_t.find("Total Sequences") == 0:
      fastq_metrics[fastqc_dir]["total_sequences"] = float(line_t.split()[2])
    elif line_t.find("#Total Duplicate Percentage") == 0:
      fastq_metrics[fastqc_dir]["duplication_percent"] = float(line_t.split()[3])
    elif line_t.find(">>Per sequence quality scores") == 0:
      state = "psqs"
    elif line_t.find(">>END_MODULE") == 0:
      state = "normal"
    elif state == "psqs":
      if line_t.find("#") < 0:
        fields = line_t.split()
        fastq_metrics[fastqc_dir]["quality"] += float(fields[1]) * float(fields[0])
  fastq_metrics[fastqc_dir]["quality"] /= fastq_metrics[fastqc_dir]["total_sequences"]

fastqc_metrics = {"gc_content": 0.0, "duplication_percent": 0.0, "quality": 0.0}
for fastqc_dir in fastq_metrics:
  weight = fastq_metrics[fastqc_dir]["total_sequences"] / float(total_reads)
  for metric in ["gc_content", "duplication_percent", "quality"]:
    fastqc_metrics[metric] += weight *  fastq_metrics[fastqc_dir][metric]

print "%g %g %g %g %g %g %g %g %g %g" % (percent_uniquely_mapped, percent_multi_mapped, percent_unmapped, snp_ti_tv, snp_het_hom, snp_indel, indel_het_hom, fastqc_metrics["quality"], fastqc_metrics["gc_content"], fastqc_metrics["duplication_percent"]), total_reads
