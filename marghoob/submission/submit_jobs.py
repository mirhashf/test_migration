#!/usr/bin/env python2.7

import argparse
from client import Client
import json
import sys
import copy
import pwd, os

parser = argparse.ArgumentParser("Submit WGS/WES jobs", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--url", metavar="URL", help="URL for portal frontend", required=True)
parser.add_argument("--username", metavar="username", help="Username", default=pwd.getpwuid(os.getuid())[0])
parser.add_argument("--password", metavar="password", help="Password", default="b")
parser.add_argument("--binabox", metavar="binabox_id", help="Bina box id", type=int, default=1)
parser.add_argument("--datasets", metavar="dataset", help="Dataset file", required=True, type=file)
parser.add_argument("--output_dir", metavar="output-dir", help="Output directory (must be on river)", required=True)
parser.add_argument("--no_submit", action="store_true", help="Don't submit. Just print the job jsons.")
parser.add_argument("--enable_vqsr", action="store_true", help="Enable VQSR")
parser.add_argument("--enable_sv", action="store_true", help="Enable SV tools")
parser.add_argument("--enable_hc", action="store_true", help="Use HaplotypeCaller for genotyping")
parser.add_argument("--aligners", metavar="aligners", nargs='+', help="List of aligners (one of bina, bwa, bwamem)", required=False, default=["bwa"]) 
parser.add_argument("--gatk_versions", metavar="gatk-versions", nargs='+', help="List of gatk versions (one of 2.7-2 or 2.3-9)", required=False, default=["2.7-2"])
parser.add_argument("--dataset_names", metavar="dataset_names", nargs='+', help="List of dataset names in the dataset json", required=False)
parser.add_argument("--keep_sorted", action="store_true", help="Keep sorted BAMs")
parser.add_argument("--keep_realigned", action="store_true", help="Keep realigned BAMs")
parser.add_argument("--keep_recalibrated", action="store_true", help="Keep recalibrated BAMs")
parser.add_argument("--sample", help="Sample name", required=True)
parser.add_argument("--library", help="Library name", required=False, default="pairedend")
parser.add_argument("--platform", help="Platform name", required=False, default="Illumina")

sys.stderr.write("Command line\n")
sys.stderr.write("%s\n" % ' '.join(sys.argv))

args = parser.parse_args()
baseurl = args.url
username = args.username
password = args.password
binabox_id = args.binabox

datasets = json.load(args.datasets)
sample = args.sample
library = args.library
platform = args.platform

dataset_names=[]
if not args.dataset_names:
  dataset_names = datasets.keys()
else:
  for dataset_name in args.dataset_names:
    if dataset_name not in datasets:
      raise Exception("Dataset %s not found in %s" % (dataset_name, args.datasets.name))
  dataset_names = args.dataset_names

output_dir = os.path.realpath(args.output_dir)
river = "/net/kodiak/volumes/river/shared/"
if output_dir.startswith(river):
  output_prefix = "river:/" + os.path.relpath(output_dir, river)
else:
  output_prefix = output_dir
if not os.path.exists(output_dir): os.makedirs(output_dir)

jobs = []

gatk_jsons = {}
dirname = os.path.dirname(os.path.realpath(__file__))
for gatk_version in args.gatk_versions:
  gatk_json_file = open(dirname + "/gatk-" + gatk_version + ".json", "r")
  gatk_jsons[gatk_version] = json.load(gatk_json_file)

workflow_common = {
   "fasta": "lake:/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa",
   "dbsnp": "lake:/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf",
   "keep_sorted_bams": args.keep_sorted,
   "keep_realigned_bams": args.keep_realigned,
   "keep_recalibrated_bams": args.keep_recalibrated,
   "bplib": "lake:/users/marghoob/resources/bplib.fa",
   "pindel_args": {
       "-T": 4
   },
   "mounts": [
   {
       "name": "river",
       "address": "kodiak:/volumes/river/shared"
   },
   {
       "name": "lake",
       "address": "kodiak:/volumes/lake/shared"
   }
   ],
   "bwa_aln_args": {
       "-t": 4,
       "-q": 30
   },
   "bwa_sampe_args": {
       "-P": ""
   },
   "bwarunner_args": {
       "--nthreads": 18,
       "--splitter_blocksize": 65536,
       "--splitter_max_pending": 4
   },
   "thread_num": 24,
   "bina_aligner_args": {
       "-p": 32,
       "-trim": 30
   },
   "bwa_mem_args": {
       "-t": 32
   }
}

for gatk_version in args.gatk_versions:
  for aligner in args.aligners:
    for dataset_name in dataset_names:
      workflow = copy.deepcopy(workflow_common)
      if aligner == "bwa": workflow["run_bwa"] = True
      elif aligner == "bwamem": workflow["run_bwa_mem"] = True
      workflow["output_prefix"] = output_prefix
      workflow["run_haplotype_caller"] = args.enable_hc
      if "alignment_groups" not in datasets[dataset_name]:
        raise Exception("alignment_groups missing for dataset %s" % (dataset_name))
      for key in datasets[dataset_name]:
        if key == "comment": continue # skip comments
        workflow[key] = datasets[dataset_name][key]
      run_type = "wes" if "bedfile" in datasets[dataset_name] else "wgs" 
      tags = "%s,%s,%s,%s,%s,%s" % (dataset_name, run_type, aligner, "vqsr" if args.enable_vqsr else "no-vqsr", "sv" if (args.enable_sv and run_type == "wgs") else "no-sv", gatk_version)
      max_mem_gb = 85 if aligner == "bina" else 105
      if gatk_version == "2.7-2": workflow["run_reduce_reads"] = True
      workflow["sorter_args"] = {"-mark_duplicates": "", "-max_mem_gb": max_mem_gb}
      for key in gatk_jsons[gatk_version]:
        workflow[key] = gatk_jsons[gatk_version][key]
      if run_type == "wes":
        workflow["vqsr_snp_train_args"] = {"--maxGaussians": 4}
      workflow["run_vqsr_snp"] = args.enable_vqsr
      workflow["run_vqsr_indel"] = (run_type == "wgs") and args.enable_vqsr
      if args.enable_sv:
        for run_svtool in ["run_breakseq", "run_pindel", "run_breakdancer", "run_cnvnator"]:
          workflow[run_svtool] = (run_type == "wgs")
      workflow["worker_num"] = 4 if run_type == "wgs" else 1
      job = {
        "workflow": copy.deepcopy(workflow),
        "bina_box":{
            "id": binabox_id
        },
        "metadata": {
            "tags": tags,
            "project": "validation",
            "library": library,
            "pi": "mm",
            "sample": sample
        },
        "workflow_type": "wgs-workflow/current/bin/wgs-workflow.jar"
      }
      jobs.append(job)

c = Client(baseurl, username, password, False)
c.login()
for job in jobs:
  if not args.no_submit:
    submitted_job = c.post(job, "job_list")
    print submitted_job["id"]
  else:
    print json.dumps(job, indent=4)
c.logout()
