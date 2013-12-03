#!/usr/bin/env python2.7

import json
import argparse
import os

parser = argparse.ArgumentParser("Generate json for different coverages", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--prefix", help="Prefix for the FASTQs", required=True)
parser.add_argument("--nlanes", help="Number of lanes", required=True, type=int)
parser.add_argument("--coverage_per_lane", help="Coverage per lane", required=True, type=float)
parser.add_argument("--bedfile", help="Regions bedfile", required=False)

args = parser.parse_args()

# convert the prefix to a path on lake/river
prefix = os.path.realpath(args.prefix)
lake = "/net/kodiak/volumes/lake/shared/"
river = "/net/kodiak/volumes/river/shared/"
if prefix.startswith(lake):
  prefix = "lake:/" + os.path.relpath(prefix, lake)
elif prefix.startswith(river):
  prefix = "river:/" + os.path.relpath(prefix, river)

# generate the jsons for different coverages
datasets = {}
for nlanes in xrange(args.nlanes):
  dataset_name = "sim_%gx" % (args.coverage_per_lane * (nlanes + 1) * 2.0)
  dataset = {}
  alignment_groups = []
  for fastqs in xrange(1 + nlanes):
    alignment_group = { "library": "pairedend", "platform": "Illumina", "read_group": "l%d" % (fastqs), "sample": "NA12878"}
    alignment_group["fastq"] = [prefix + "%d.bwa.read1.fastq.gz" % (fastqs), prefix + "%d.bwa.read2.fastq.gz" % (fastqs)]
    alignment_groups.append(alignment_group)
  datasets[dataset_name] = {"alignment_groups": alignment_groups}
  if args.bedfile:
    bedfile = os.path.realpath(args.bedfile)
    if bedfile.startswith(lake): bedfile = "lake:/" + os.path.relpath(bedfile, lake)
    elif bedfile.startswith(river): bedfile = "river:/" + os.path.relpath(bedfile, river)
    datasets[dataset_name]["bedfile"] = bedfile

print json.dumps(datasets, indent=4)
