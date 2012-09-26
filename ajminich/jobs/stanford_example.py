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

data_dir = "bina://data/public/snyder"

# Set up the libraries
libraries = [
             ("blood", "A804NLABXX", 8)
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

# Reference
job.reference.set_species("human")
job.reference.set_genome_build("CEUref")
job.reference.set_dbsnp_build("135_major")

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

# You can retain the sorted BAM files for analysis. However, this is
# generally not recommended due to limited space available on the system.
job.alignment.set_keep_sorted_bam(False)

'''
        REALIGNMENT
'''

# To configure the realignment:
# job.realignment.fast_realigner.set_argument(<argument>, <boolean>)
# job.realignment.fast_realigner.set_option(<option>, <value>)

'''
        GENOTYPING
'''

# To configure the genotyping:
# job.genotyping.fast_genotyper.set_argument(<argument>, <boolean>)
# job.genotyping.fast_genotyper.set_option(<option>, <value>)

'''
        STRUCTURAL VARIATION
'''

# Enable all structural variation tools
job.structural_variation.set_disable_bina_sv(False)
job.structural_variation.set_run_breakdancer(True)
job.structural_variation.set_run_breakseq(True)
job.structural_variation.set_run_cnvnator(True)
job.structural_variation.set_run_pindel(True)
job.structural_variation.pindel.set_use_breakdancer(True)

'''
        VQSR
'''

job.genotyping.set_perform_vqsr(True)

# Set up variant quality training resources.
# Assumes that the bina:// path is pointed to /srv/gs1/projects/snyder/cuiping.

hapmap_resource = bina.VariantRecalibrationResource("hapmap",
    "bina://data/referencefiles/hapmap/hapmap_3.3.hg19.majorchr.sites.vcf")
hapmap_resource.set_prior(15.0)
hapmap_resource.set_known(False)
hapmap_resource.set_training(True)
hapmap_resource.set_truth(True)

omni_resource = bina.VariantRecalibrationResource("omni",
    "bina://data/referencefiles/omni/1000G_omni2.5.hg19.majorchr.sites.vcf")
omni_resource.set_prior(12.0)
omni_resource.set_known(False)
omni_resource.set_training(True)
omni_resource.set_truth(False)

dbsnp_resource = bina.VariantRecalibrationResource("dbsnp",
    "bina://data/referencefiles/dbsnp/dbsnp_135.majorchr.hg19.vcf")
omni_resource.set_prior(6.0)
omni_resource.set_known(True)
omni_resource.set_training(False)
omni_resource.set_truth(False)

mills_resource = bina.VariantRecalibrationResource("mills",
    "bina://data/referencefiles/indel/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")
mills_resource.set_prior(12.0)
mills_resource.set_known(True)
mills_resource.set_training(True)
mills_resource.set_truth(True)

# SNP Recalibration
recal_operation = bina.VariantRecalOperation()
recal_operation.set_name("SNP recalibration")
recal_operation.set_variant_type("SNP")

recal_operation.variant_recalibrator.set_option("--use_annotation",
    "QD,HaplotypeScore,MQRankSum,ReadPosRankSum,HRun")
recal_operation.apply_recalibration.set_option("--ts_filter_level", 99.0)

recal_operation.add_resource(hapmap_resource)
recal_operation.add_resource(omni_resource)
recal_operation.add_resource(dbsnp_resource)

job.genotyping.add_recal_operation(recal_operation)

# Indel Recalibration based on 1000 Genomes index variants
recal_operation = bina.VariantRecalOperation()
recal_operation.set_name("indel recalibration")
recal_operation.set_variant_type("INDEL")

recal_operation.variant_recalibrator.set_option("--use_annotation",
    "QD,FS,HaplotypeScore,ReadPosRankSum")
recal_operation.apply_recalibration.set_option("--ts_filter_level", 99.0)

recal_operation.add_resource(mills_resource)

job.genotyping.add_recal_operation(recal_operation)

'''
        JOB SUBMISSION
'''

# Connect to Bina Box using API key and IP addresses of nodes
binabox = bina.BinaBox()
binabox.connect(api_key, nodes)

# Submit the job
job_id = binabox.run_job(job)
print "Job submitted with ID " + job_id + "."

# Close the connection
binabox.close()
