#!/usr/bin/env python2.7

import argparse
import json
import sys

parser = argparse.ArgumentParser("Submit WGS/WES jobs", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--dataset_file", metavar="dataset_file", help="Dataset file", required=True, type=file)

args = parser.parse_args()

datasets = json.load(args.dataset_file)

dataset_keys = datasets.keys()
dataset_keys.sort()
max_key_len = max([len(key) for key in dataset_keys])

fmt = "%%-%ds| %%s" % (max_key_len + 1)
print fmt % ("Dataset name", "Comment")
print "-"*120
for dataset_name in dataset_keys:
  dataset = datasets[dataset_name]
  comment = dataset["comment"] if "comment" in dataset else ""
  print fmt % (dataset_name, comment)
