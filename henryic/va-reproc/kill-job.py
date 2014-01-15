#!/usr/bin/env python

# Simple script to submit a multisample job
from client import Client
import json
import sys
import ConfigParser

config = ConfigParser.ConfigParser();
config.read("config.cfg")

# Account credentials
baseurl = config.get("bina", "api-url")
username = config.get("bina", "username")
password = config.get("bina", "password")
binabox_id = config.get("bina", "box-id")

request = {
           "kill": True
}

if len(sys.argv) < 2:
    print "Usage: %s [JOB ID]" % sys.argv[0]
    exit(1)

# Main code
c = Client(baseurl, username, password, False)
c.login()
job = c.put(request, "running_job_list/%s/%s" % (binabox_id, sys.argv[1]))
print job
c.logout()