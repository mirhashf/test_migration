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

# Main code
c = Client(baseurl, username, password, False)
c.login()
job = c.get("mount")
print job
c.logout()