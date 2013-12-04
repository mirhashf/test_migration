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
parser.add_argument("--no_submit", action="store_true")
parser.add_argument("--enable_vqsr", action="store_true", help="Enable VQSR")
parser.add_argument("--enable_sv", action="store_true", help="Enable SV tools")

sys.stderr.write("Command line\n")
sys.stderr.write("%s\n" % ' '.join(sys.argv))

args = parser.parse_args()
baseurl = args.url
username = args.username
password = args.password
binabox_id = args.binabox

gatk_versions = ["2.3-9"] #, "2.7-2"]

datasets = json.load(args.datasets)
sample = 'NA12878'
library = 'pairedend'
platform = "Illumina"

output_dir = os.path.realpath(args.output_dir)
river = "/net/kodiak/volumes/river/shared/"
if output_dir.startswith(river):
  output_prefix = "river:/" + os.path.relpath(output_dir, river)
else:
  output_prefix = output_dir

jobs = []

gatk_jsons = {}
dirname = os.path.dirname(os.path.realpath(__file__))
for gatk_version in gatk_versions:
  gatk_json_file = open(dirname + "/gatk-" + gatk_version + ".json", "r")
  gatk_jsons[gatk_version] = json.load(gatk_json_file)

workflow_common = {
   "fasta": "lake:/users/marghoob/GATK-bundle-hg19/ucsc.hg19.fa",
   "dbsnp": "lake:/users/marghoob/GATK-bundle-hg19/dbsnp_137.hg19.vcf",
   "keep_sorted_bams": False,
   "keep_realigned_bams": False,
   "keep_recalibrated_bams": False,
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

for gatk_version in gatk_versions:
  for run_bwa in [True]: #, False]:
    for dataset_name in datasets:
      workflow = workflow_common
      workflow["output_prefix"] = output_prefix
      if "alignment_groups" not in datasets[dataset_name]:
        raise Exception("alignment_groups missing for dataset %s" % (dataset_name))
      for key in datasets[dataset_name]:
        workflow[key] = datasets[dataset_name][key]
      run_type = "wes" if "bedfile" in datasets[dataset_name] else "wgs" 
      tags = "%s,%s,%s,vqsr,%s" % (dataset_name, run_type, "bwa" if run_bwa else "bina", gatk_version)
      workflow["run_bwa"] = run_bwa
      max_mem_gb = 105 if run_bwa else 85
      if gatk_version == "2.7-2": workflow["run_reduce_reads"] = True
      workflow["sorter_args"] = {"-mark_duplicates": "", "-max_mem_gb": max_mem_gb}
      for key in gatk_jsons[gatk_version]:
        workflow[key] = gatk_jsons[gatk_version][key]
      run_vqsr_indel = True
      if run_type == "wes":
        workflow["vqsr_snp_train_args"] = {"--maxGaussians": 4}
        run_vqsr_indel = False
      if args.enable_vqsr:
        workflow["run_vqsr_snp"] = True
        workflow["run_vqsr_indel"] = run_vqsr_indel
      if args.enable_sv:
        for run_svtool in ["run_breakseq", "run_pindel", "run_breakdancer", "run_cnvnator"]:
          workflow[run_svtool] = (run_type == "wgs")
      workflow["worker_num"] = 1 if run_type == "wgs" else 1
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
  
  #print json.dumps(job, indent=True)
  if not args.no_submit:
    submitted_job = c.post(job, "job_list")
    print submitted_job["id"]
  else:
    print json.dumps(job, indent=4)
c.logout()
