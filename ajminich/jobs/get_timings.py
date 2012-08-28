#!/usr/bin/env python
'''
Gets the timings for multiple jobs, as a CSV file.

Assumes the seqalto repository has been cloned into the same folder as the
sandbox repository.
'''

from datetime import datetime
import glob
import json
import os
import sys

# Import bina module from relative path
module_dir = os.path.abspath(os.path.join(os.path.realpath(__file__), '../../../../seqalto/loomis/sdk'))
sys.path.append(module_dir)

import bina

# Settings
status_pattern = "/home/ajminich/krakow-statuses/*.status"
output_file = "timings.csv"

header = [
    "Job", "Status", "Queue Time",
    "Started", "Alignment", "Realignment", "Genotyping", "Total"
    ]

if __name__ == '__main__':

    out = open(output_file, 'w')
    out.write(",".join(header) + "\n")

    for status_file in glob.glob(status_pattern):
        
        print "Now processing '" + status_file + "'."
        
        text = open(status_file, 'r').read()
        status_entry = bina.JobStatus(json.loads(text))
        job_entry = status_entry.get_job()
        
        logs = status_entry.get_logs()
        start = None
        alignment = None
        realignment = None
        genotyping = None
        total = None
        
        for log in logs:
            
            log_time = datetime.strptime(log.get_timestamp(), "%b/%d/%Y %I:%M:%S %p")
            
            if "Starting alignment on lane" in log.get_description():
                start = log_time
            elif "Completed alignment" in log.get_description():
                alignment = log_time
            elif "Indel realignment complete" in log.get_description():
                realignment = log_time
            elif "Genotyping complete" in log.get_description():
                genotyping = log_time
            elif "Bina job finished successfully" in log.get_description():
                total = log_time
                
        if not total:
            continue
        
        entries = [
            job_entry.get_description(),
            status_entry.get_status(),
            status_entry.get_queue_date_as_string(),
            datetime.strftime(start, "%b/%d/%Y %I:%M:%S %p"),
            str(alignment - start),
            str(realignment - alignment),
            str(genotyping - realignment),
            str(total - start)
            ]
        
        out.write(",".join(entries) + "\n")
        
    out.close()
    