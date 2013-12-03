#!/usr/bin/env python2.7

import argparse
from client import Client
import json
import sys
import copy
import fileinput
import signal, os
import time
import pwd

parser = argparse.ArgumentParser("Monitor WGS/WES jobs", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--url", metavar="URL", help="URL for portal frontend", required=True)
parser.add_argument("--username", metavar="username", help="Username", default=pwd.getpwuid(os.getuid())[0])
parser.add_argument("--password", metavar="password", help="Password", default="b")
parser.add_argument("--binabox", metavar="binabox_id", help="Bina box id", type=int, default=1)

args = parser.parse_args()

baseurl = args.url
username = args.username
password = args.password
binabox_id = args.binabox

job_ids = set()
finished_ids = set()
for line in sys.stdin:
  job_ids.add(line.strip())

print job_ids

c = Client(baseurl, username, password, False)
c.login()

def print_finished_status(finished_ids):
  box_finished_job_list = c.get("job_list/2")
  box_finished = set([(box_finished_job["id"], box_finished_job["status"]) for box_finished_job in box_finished_job_list])
  for job_id in job_ids:
    if (job_id, "Successful") in box_finished: print job_id, "Successful"
    elif (job_id, "Failed") in box_finished: print job_id, "Failed"
    else: print job_id, "Killed/Queued/Running"

def int_handler(signum, frame):
  print "Logging out"
  print finished_ids
  print_finished_status(finished_ids)
  c.logout()
  sys.exit(0)

signal.signal(signal.SIGINT, int_handler)

while job_ids:
  box_jobs = c.get("running_job_list/%s" % (binabox_id))
  box_jobs_ids = set([box_job["id"] for box_job in box_jobs])
  finished_ids |= job_ids - box_jobs_ids
  if not box_jobs_ids.intersection(job_ids): break
  print box_jobs_ids
  time.sleep(60)
print_finished_status(finished_ids)
c.logout()
