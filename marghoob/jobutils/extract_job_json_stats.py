#!/usr/bin/python

import json
import sys
import re

if len(sys.argv) != 2:
  sys.stderr.write("./script <job_json>\n")
  sys.exit(1);

job_json = json.load(open(sys.argv[1], "r"))

workflow_json = job_json["config"]["meta"]["workflow"]
worker_num = workflow_json["worker_num"]
thread_num = workflow_json["thread_num"]

workers_json = job_json["workers"]
def get_time_taken(workers_json):
  all_tasks = [task for worker_json in workers_json for task in worker_json["tasks"]]
  re_pattern = re.compile("^W(\d)+\.INDEXGET.*")
  start_time = max([task["endTime"] for task in all_tasks if re_pattern.match(task["name"])])
  end_time = max([task["endTime"] for task in all_tasks])
  return (end_time - start_time)/1e3/3600

time_taken = get_time_taken(workers_json)
has_avx = ("haplotype_caller_args" in workflow_json) and ("-pairHMM" in workflow_json["haplotype_caller_args"])
is_gatk_3_1 = ("gatk_path" not in workflow_json) or (re.match(".*gatk-3\.1-1-.*", workflow_json["gatk_path"]) is not None) or (re.match(".*GenomeAnalysisTK-2014\.2-3\.1\.7-7-gb0ce0a5.*", workflow_json["gatk_path"]) is not None)
is_gatk_2_8 = ("gatk_path" in workflow_json) and ((re.match(".*gatk-2\.8.*", workflow_json["gatk_path"]) is not None) or (re.match(".*GenomeAnalysisTK-2014\.1-.*", workflow_json["gatk_path"]) is not None))
aligner = "bwa" if ("run_bwa" in workflow_json and workflow_json["run_bwa"]) else ("bwamem" if ("run_bwa_bem" in workflow_json and workflow_json["run_bwa_bem"]) else "bina")

gatk_version = "3.1" if is_gatk_3_1 else ("2.8" if is_gatk_2_8 else "unknown")
is_desktop = workflow_json["sorter_args"]["-max_mem_gb"] > 150
data_type = "wgs" if "bedfile" not in workflow_json else "wes"

print worker_num, thread_num, time_taken, has_avx, gatk_version, ("desktop" if is_desktop else "enterprise"), aligner, data_type, "dummy"
