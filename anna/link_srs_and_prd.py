import json
import csv
from jira import JIRA
import ConfigParser
import sys
import getpass

# pass in the trace matrix csv as the second argument



def parse_trace_file():
	trace_matrix = sys.argv[1]
	with open(trace_matrix, 'r') as f:
		# header = f.readline()

		data = [row for row in csv.reader(f.read().splitlines())]
		data = data[1:]
		
	return data


def link_issues(trace_matrix, jira):

	# data[1] refers to PRD, data[3] refers to SRS
	for link in trace_matrix:



		prd = link[0]
		srs = link[1]

		print link
		print prd
		print srs

		# need edge case if prd does not exist yet in jira
		# right now assume that it exists
		prd_issue = jira.search_issues('project= SBX and cf[12301]= ' +  str(prd))

		print prd_issue

		srs_issue = jira.search_issues('project= SBX and cf[12301]= ' +  str(srs))

		print srs_issue

		jira.create_issue_link("Relates", prd_issue[0].key, srs_issue[0].key, None)

		print "linked issue %s and %s" % (prd_issue, srs_issue)




def login_to_jira():
	user = raw_input("JIRA Username: ")
	passwd = getpass.getpass()
	return user, passwd


def main():

	trace_matrix = parse_trace_file()
	print trace_matrix

	config = ConfigParser.ConfigParser()
	config.read("auth.cfg")
	user, password = login_to_jira()

	host = config.get("JIRA", "host")
	options = {'server': host}
	jira = JIRA(options, basic_auth=(user, password))

	link_issues(trace_matrix, jira)





if __name__ == "__main__":
	main()