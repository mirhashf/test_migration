#!/usr/bin/env python
'''
An example multi-lane job for Cuiping Pan at Stanford Genetics Center.

Created: 8/22/2012
Author: Bina Technologies, Inc.
'''

import os
import sys

# Import bina module from relative path
module_dir = os.path.abspath(os.path.join(os.path.realpath(__file__), '../..'))
sys.path.append(module_dir)

import bina

'''
        CONFIGURATION
'''

api_key = "gocardinal"
nodes = [
         "scg-bb-03-genstorage.sunet:8080",
         "scg-bb-05-genstorage.sunet:8080",
         "scg-bb-02-genstorage.sunet:8080",
         "scg-bb-04-genstorage.sunet:8080",
         "scg-bb-01-genstorage.sunet:8080"
         ]

data_dir = "bina://data/snyder"

# Set up the libraries
libraries = [
             ("blood", "A804NLABXX", 8),
             ("saliva", "A806WMABXX", 8),
             ("saliva", "A808HKABXX", 8)
             ]
sample = "snyder"
   
'''
        JOB SETUP
'''
   
# Create a new job
job = bina.Job()

# Set up the job
job.set_output_dir("bina://out")
job.set_description("Snyder 120x run with Bina Pipeline")
job.set_use_broad_gatk(False)

# Create alignment tasks for the Bina Aligner
for library in libraries:
    
    library_name = library[0]
    lane_prefix = library[1]
    num_lanes = library[2]
    
    print "Now queueing library '" + library_name + "' with " + str(num_lanes) + " and prefix '" + lane_prefix + "'." 
    
    for lane_index in range(1, num_lanes+1):
        
        aligner_job = bina.BinaAlignerJob(
            first_end = data_dir + "/" + lane_prefix + ".s_" + str(lane_index) + "_1.fq.gz",
            second_end = data_dir + "/" + lane_prefix + ".s_" + str(lane_index) + "_2.fq.gz",
            readgroup = lane_prefix + ".s_" + str(lane_index),
            library = library_name,
            sample = sample)
        aligner_job.set_trimming(0)
        
        # Set aligner template size calculation to automatic
        aligner_job.set_option("--template_len_comp_method", 2)
        
        # Use batched template size calculation to emulate BWA
        aligner_job.set_argument("--enable_batch_pairing", True)
    
        job.alignment.add_aligner_job(aligner_job)

# For extremely high coverage, we disable concurrent sorting to increase 
# the amount of sorting we can perform in main memory. 
job.alignment.set_disable_concurrent_sorting(True)

'''
        JOB SUBMISSION
'''

# Connect to Bina Box using API key and IP addresses of nodes
binabox = bina.BinaBox()
#binabox.connect(api_key, nodes)

# Submit the job
# job_id = binabox.run_job(job)
#print "Job submitted with ID " + job_id + "."

print str(job)

# Close the connection
binabox.close()
