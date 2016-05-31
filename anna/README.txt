create_srs_jira_tickets.py takes another command line argument for the csv requirements file that should be an export from HP ALM.

If you run the script you will be prompted for your atlassian username and password.
The script creates new tickets for those items in the csv that have a Req ID that has not been seen before.
The script edits existing tickets if edits to the description have been made.

Need to do:
Pass argument for the project that we wish to add the tickets to (right now it adds ticketes to SBX)
See if there are changes to any other fields of an existing ticket
Potentially delete tickets if they have been deleted from HP ALM

