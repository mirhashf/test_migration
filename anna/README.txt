create_tickets_for_reqs.py takes a file name for the requirements (either a csv for SRSs or a csv for PRDs), then the argument of whether it's a PRD or SRS.

If you run the script you will be prompted for your atlassian username and password.
The script creates new tickets for those items in the csv that have a Req ID that has not been seen before.
The script edits existing tickets if edits to the description have been made.

Need to do:
Pass argument for the project that we wish to add the tickets to (right now it adds ticketes to SBX)
Potentially delete tickets if they have been deleted from HP ALM
