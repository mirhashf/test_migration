#!/usr/bin/env python
from jira import JIRA
import ConfigParser

config = ConfigParser.ConfigParser()
config.read("auth.cfg")
password = config.get("JIRA", "pass")
user = config.get("JIRA", "user")
host = config.get("JIRA", "host")
options = {'server': host}

jira = JIRA(options, basic_auth=(user, password))

# get all issues assigned to amir
issues = jira.search_issues('assignee=' + user)

print issues
