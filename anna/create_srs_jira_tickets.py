#!/usr/bin/env python
import json
import csv
from jira import JIRA
import ConfigParser
import sys
import getpass

# notes
# can also have TRQ or PRQ as a parameter in the command line

def get_file_from_cmd_line():
	if (len(sys.argv) != 2):
		sys.exit("Please inclue the csv file name of srs requirements as an argument.")

	file_name = sys.argv[1]
	return file_name



def read_file(csv_file):
    with open(csv_file, 'r') as f:
        data = [row for row in csv.reader(f.read().splitlines())]
    return data

def isNumber(req_id):
	try: 
		int(req_id)
		return True
	except ValueError:
		return False

def save_srs_data_to_ticket(srs_data, jira):

	for srs in srs_data:

		req_id = srs[0]
		if (not isNumber(req_id)): continue
		req_id = int(req_id)
		req_summary = srs[1]
		req_description = srs[2]
		req_component = srs[3]
		# the name that we give the issue
		# req_summary = str(req_id) + '-' + req_component

		# check if issue exists in jira
		issue = jira.search_issues('project= SBX and cf[12301]= ' +  str(req_id))

		if len(issue) > 0 :

			# need to compare current values to new values
			# version is hard coded
			update_issue(issue[0], req_summary, req_description, 12)

		else:

			new_issue = jira.create_issue(project='SBX', summary= req_summary,
                              description=req_description, issuetype={'name': 'Requirement'})

			# 12300 = HPALM Project name, 12301 = HP ALM Req ID, 12303 = HP ALM Version
			new_issue.update(fields={'customfield_12300': 'RSU_Development', 'customfield_12301': req_id, 'customfield_12303': 12})
	
def update_summary(issue, new_summary):
	current_summary = issue.fields.summary
	if not (current_summary == new_summary):
		print "updating summary"
		issue.update(fields={'summary': new_summary})
	else:
		print "keeping the same summary"

def update_description(issue, new_descript):
	current_description = issue.fields.description

	if not (current_description == new_descript):
		print "updating description "
		issue.update(fields={'description': new_descript})
	else:
		print "keeping the same description"

def update_version(issue, new_version):
	current_version = issue.fields.customfield_12303
	if not (current_version == new_version):
		print "updating version"
		issue.update(fields={'customfield_12303': new_version})

	else:
		print "keeping the same version"


def update_issue(issue, new_summary, new_descript, new_version):
	update_summary(issue, new_summary)
	update_description(issue, new_descript)
	update_version(issue, new_version)
	


def delete_all_issues(jira):
	issues_in_proj = jira.search_issues('project=SBX')
	for issue in issues_in_proj:
		print "deleting issue"
		issue.delete()

def login_to_jira():
	user = raw_input("JIRA Username: ")
	passwd = getpass.getpass()
	return user, passwd


def main():

	srs_data = read_file(get_file_from_cmd_line())

	# set up JIRA instance

	config = ConfigParser.ConfigParser()
	config.read("auth.cfg")
	user, password = login_to_jira()

	host = config.get("JIRA", "host")
	options = {'server': host}
	jira = JIRA(options, basic_auth=(user, password))


	save_srs_data_to_ticket(srs_data, jira)








if __name__ == "__main__":
    main()