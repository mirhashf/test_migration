#!/usr/bin/python

import json
import sys
import os
import argparse
import glob
import shutil

dir_list = sys.argv[1:]

job_dirs = [os.path.realpath(jobdirname) for dirname in dir_list for jobdirname in glob.glob("%s/*" % (dirname))]

mounts = ["/net/kodiak/volumes/river/shared/", "/net/kodiak/volumes/delta/shared/"]

def get_dot_bina_dir(jobdir):
  for mount in mounts:
    if jobdir.find(mount) == 0:
      return os.path.join(mount, ".bina/%s" % (os.path.basename(jobdir)))
  return None

def was_job_successful(job_json):
  total_tasks = 0
  failed_tasks = 0
  for worker_json in job_json["workers"]:
    tasks_json = worker_json["tasks"]
    for task_json in tasks_json:
      total_tasks = total_tasks + 1
      task_status = task_json["status"]
      if task_status != "FINISHED": failed_tasks = failed_tasks + 1
  return failed_tasks == 0

def completely_remove_job(jobdir):
  bina_job_dir = get_dot_bina_dir(jobdir)
  print "Will delete %s %s" % (bina_job_dir, jobdir)
  if os.path.isdir(bina_job_dir): shutil.rmtree(bina_job_dir, True)
  if os.path.isdir(jobdir): shutil.rmtree(jobdir, True)

for job_dir in job_dirs:
  if os.path.isdir(job_dir) and os.path.isfile(os.path.join(job_dir, "job.json")):
    with open(os.path.join(job_dir, "job.json")) as job_json_file:
      job_json = json.load(job_json_file)
    print job_dir, get_dot_bina_dir(job_dir), was_job_successful(job_json)
    if not was_job_successful(job_json):
      print job_dir, " will be removed"
      completely_remove_job(job_dir)
