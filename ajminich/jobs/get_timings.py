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
module_dir = os.path.abspath(os.path.join(os.path.realpath(__file__), '../../../../seqalto/loomis/client'))
sys.path.append(module_dir)

import bina

# Settings
status_pattern = "/home/ajminich/trex/stanford-statuses/*.status"
output_file = "timings.csv"

header = [
    "Job", "Status", "Queued",
    "Start Date", "Start Time",
    "Finish Date", "Finish Time",
    "Alignment", "Realignment", "Genotyping", "Elapsed"
    ]

if __name__ == '__main__':

    out = open(output_file, 'w')
    out.write(",".join(header) + "\n")

    for status_file in glob.glob(status_pattern):
        
        text = open(status_file, 'r').read()
        status_entry = bina.JobStatus(json.loads(text))
        job_entry = status_entry.get_job()
        
        if status_entry.get_status().lower() not in ["successful"]:
            continue
        
        logs = status_entry.get_logs()
        start = None
        alignment = None
        realignment = None
        genotyping = None
        total = None
        
        for log in logs:
            
            log_time = datetime.strptime(log.get_timestamp(), "%b/%d/%Y %I:%M:%S %p")
            
            if "Starting" in log.get_description() and not start:
                start = log_time
            elif "Completed alignment" in log.get_description():
                alignment = log_time
            elif "Indel realignment complete" in log.get_description():
                realignment = log_time
            elif "Genotyping complete" in log.get_description():
                genotyping = log_time
            elif "Bina job finished successfully" in log.get_description():
                total = log_time
                
        if not total or not start:
            continue
        
        alignment_time = str(alignment - start) if alignment else "N/A"
        realignment_time = str(realignment - alignment) if realignment else "N/A"
        genotyping_time = str(genotyping - realignment) if genotyping else "N/A"
        
        entries = [
            job_entry.get_description(),
            status_entry.get_status(),
            status_entry.get_queue_date_as_string(),
            datetime.strftime(start, "%m/%d/%Y"),
            datetime.strftime(start, "%I:%M:%S %p"),
            datetime.strftime(total, "%m/%d/%Y"),
            datetime.strftime(total, "%I:%M:%S %p"),
            alignment_time,
            realignment_time,
            genotyping_time,
            str(total - start)
            ]
        
        out.write(",".join(entries) + "\n")
        
    out.close()
    