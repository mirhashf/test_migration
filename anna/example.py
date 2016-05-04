#!/usr/bin/env python
import json
import csv
from jira import JIRA
import ConfigParser

# breakdown of how to set requirement ids for user stories:
# 1. import csv file
# 2. create a dictionary from user stories to requirement ids
# 3. iterate through the keys of the dictionary (user stories)
# 4. save the list of requirement ids to the customfield_12100

def read_file(csv_file):
    with open(csv_file, 'r') as f:
        data = [row for row in csv.reader(f.read().splitlines())]
    return data

def get_user_stories_to_ids(csv_file):
	data = read_file(csv_file)
	user_stories_to_ids = dict()

	for row in data:
		# print row
		if row[1] in user_stories_to_ids:
			if (row[0] != ''):
				user_stories_to_ids[row[1]].append('RSUDev' + row[0])
		else:
			user_stories_to_ids[row[1]] = []
			if (row[0] != ''):
				user_stories_to_ids[row[1]].append('RSUDev' + row[0])


	return user_stories_to_ids
		# print row

	




def main():

	config = ConfigParser.ConfigParser()
	config.read("auth.cfg")
	password = config.get("JIRA", "pass")
	user = config.get("JIRA", "user")
	host = config.get("JIRA", "host")
	options = {'server': host}
	jira = JIRA(options, basic_auth=(user, password))

	user_stories_to_ids = get_user_stories_to_ids('../../../mappings.csv')

	print "user stories to ids"
	print user_stories_to_ids



	for user_story in user_stories_to_ids:

		print "before update"
		issue = jira.issue(user_story)
		print issue.fields.customfield_12100



		print "after update"
		issue.update(fields={'customfield_12100': user_stories_to_ids[user_story]})
		print issue.fields.customfield_12100








if __name__ == "__main__":
    main()