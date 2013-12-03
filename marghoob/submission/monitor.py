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
parser.add_argument("--monitor_interval", metavar="monitor-interval", help="Monitor interval in seconds", type=int, default=60)

args = parser.parse_args()

baseurl = args.url
username = args.username
password = args.password
binabox_id = args.binabox

running_job_ids = set()
finished_ids = set()
for line in sys.stdin:
  running_job_ids.add(line.strip())

c = Client(baseurl, username, password, False)
c.login()

def print_finished_status(finished_ids):
  box_finished_job_list = c.get("job_list/2")
  box_finished = {}
  for box_finished_job in box_finished_job_list:
    box_finished[box_finished_job["id"]] = box_finished_job["status"]
  for job_id in finished_ids:
    if job_id in box_finished:
      print job_id, box_finished[job_id]
    else:
      print job_id, "Killed/Queued/Running"

def int_handler(signum, frame):
  print_finished_status(finished_ids)
  c.logout()
  sys.exit(0)

while running_job_ids:
  box_jobs = c.get("running_job_list/%s" % (binabox_id))
  box_jobs_ids = set([box_job["id"] for box_job in box_jobs])
  finished_ids |= running_job_ids - box_jobs_ids
  running_job_ids &= box_jobs_ids
  if not running_job_ids: break
  sys.stderr.write("Jobs on box: %s\n" % (box_jobs_ids))
  time.sleep(args.monitor_interval)
print_finished_status(finished_ids)
c.logout()
