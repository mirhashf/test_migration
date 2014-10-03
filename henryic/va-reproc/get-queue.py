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

if len(sys.argv) < 2:
    print "Usage: %s [BINABOX ID]" % sys.argv[0]
    exit(1)

# Main code
c = Client(baseurl, username, password, False)
c.login()
job = c.get("running_job_list/%s" % (binabox_id))
print job
c.logout()