#!/usr/bin/env python
'''
Runs Snyder samples.

Created: 8/22/2012
Author: Bina Technologies, Inc.
'''

import os
import sys

# Import bina module from relative path
module_dir = os.path.abspath("/home/ajminich/seqalto/loomis/client")
sys.path.append(module_dir)

import bina

'''
        CONFIGURATION
'''

api_key = "gocardinal"
nodes = [
    "tehran-00:8080",
    "tehran-01:8080",
    "tehran-02:8080",
    "tehran-03:8080",
    "tehran-06:8080"
]

data_dir = "bina://data/snyder"

# Set up the libraries
libraries = [ ("blood", "A804NLABXX", 8) ]
sample = "snyder"

'''
        JOB SETUP
'''

# Create a new job
job = bina.Job()

# Set up the job
job.set_output_dir("bina://output/snyder")
job.set_description("Snyder Blood with Reference Calls")
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
            readgroup = "lane" + str(lane_index),
            library = library_name,
            sample = sample)
        aligner_job.set_trimming(0)

        # Snyder reads setting
        aligner_job.set_option("-m", 330)

        # Set aligner template size calculation to automatic
        #aligner_job.set_option("--template_len_comp_method", 2)

        # Use batched template size calculation to emulate BWA
        #aligner_job.set_argument("--enable_batch_pairing", True)

        aligner_job.set_option("-h", -1)
        aligner_job.set_argument("--nw_disable_match_at_ends", True)

        job.alignment.add_aligner_job(aligner_job)

job.alignment.set_disable_concurrent_sorting(False)
job.alignment.set_keep_sorted_bam(False)
job.set_keep_scratch_files(False)

job.genotyping.fast_genotyper.set_option("--output_mode", "EMIT_ALL_CONFIDENT_SITES")

'''
        VARIANT QUALITY SCORE RECALIBRATION
'''

job.genotyping.set_perform_vqsr(False)

hapmap_resource = bina.VariantRecalibrationResource("hapmap",
    "bina://data/variants/sites/hapmap_3.3.hg19.vcf")
hapmap_resource.set_prior(15.0)
hapmap_resource.set_known(False)
hapmap_resource.set_training(True)
hapmap_resource.set_truth(True)

omni_resource = bina.VariantRecalibrationResource("omni",
    "bina://data/variants/sites/1000G_omni2.5.hg19.vcf")
omni_resource.set_prior(12.0)
omni_resource.set_known(False)
omni_resource.set_training(True)
omni_resource.set_truth(False)

dbsnp_resource = bina.VariantRecalibrationResource("dbsnp",
    "bina://data/variants/sites/dbsnp_135.hg19.vcf")
omni_resource.set_prior(6.0)
omni_resource.set_known(True)
omni_resource.set_training(False)
omni_resource.set_truth(False)

# SNP Recalibration
recal_operation = bina.VariantRecalOperation()
recal_operation.set_name("SNP recalibration")
recal_operation.set_variant_type("SNP")

recal_operation.variant_recalibrator.set_option("--use_annotation",
    "QD,HaplotypeScore,MQRankSum,ReadPosRankSum,FS,MQ,DP")
recal_operation.apply_recalibration.set_option("--ts_filter_level", 99.0)

recal_operation.add_resource(hapmap_resource)
recal_operation.add_resource(omni_resource)
recal_operation.add_resource(dbsnp_resource)

job.genotyping.add_recal_operation(recal_operation)

# Indel Recalibration
recal_operation = bina.VariantRecalOperation()
recal_operation.set_name("indel recalibration")
recal_operation.set_variant_type("INDEL")

recal_operation.variant_recalibrator.set_option("--use_annotation",
    "QD,HaplotypeScore,FS")
recal_operation.apply_recalibration.set_option("--ts_filter_level", 99.0)

recal_operation.add_resource(hapmap_resource)
recal_operation.add_resource(omni_resource)
recal_operation.add_resource(dbsnp_resource)

job.genotyping.add_recal_operation(recal_operation)

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
        JOB SUBMISSION
'''

# Connect to Bina Box using API key and IP addresses of nodes
binabox = bina.BinaBox()
binabox.connect(api_key, nodes)

# Submit the job
job_id = binabox.run_job(job)
print "Job submitted with ID: " + job_id

# Close the connection
binabox.close()
