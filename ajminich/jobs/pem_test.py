#!/usr/bin/env python
'''
Runs a simulated dataset test.
'''

import os
import sys

# Import bina module from relative path
module_dir = '/home/ajminich/seqalto/loomis/client'
sys.path.append(module_dir)

import bina

'''
        CONFIGURATION
'''

api_key = "gocardinal"
servers = [ "t-rex:1380" ]
managers = [ "t-rex:1381" ]

dataset = "dsTest"
datasets_dir = "bina://data"

'''
        JOB SETUP
'''
   
# Create a new job
job = bina.Job()

# Set up the job
job.set_output_dir("bina://jobs/" + dataset)
job.set_description("Loomis PEM test on " + dataset)
job.set_use_broad_gatk(False)

# Reference
job.reference.set_species("human")
job.reference.set_genome_build("chr19_22")
job.reference.set_dbsnp_build("CEU-1409.15_21.chr")
#job.set_enable_local_run(True)

#job.structural_variation.bina_sv.set_perform_pem(False)

# Create alignment tasks for the Bina Aligner
aligner_job = bina.BinaAlignerJob(
            first_end = datasets_dir + "/" + dataset + "/" + dataset + "_1.fq",
            second_end = datasets_dir + "/" + dataset + "/" + dataset + "_2.fq",
            readgroup = dataset,
            library = dataset)
aligner_job.set_trimming(30)

job.alignment.add_aligner_job(aligner_job)

# For extremely high coverage, we disable concurrent sorting to increase 
# the amount of sorting we can perform in main memory. 
job.alignment.set_disable_concurrent_sorting(False)

# You can retain the sorted BAM files for analysis. However, this is
# generally not recommended due to limited space available on the system.
#job.alignment.set_keep_sorted_bam(True)
#job.set_keep_scratch_files(True)

'''
        REALIGNMENT
'''

# To configure the realignment:
# job.realignment.fast_realigner.set_argument(<argument>, <boolean>)
# job.realignment.fast_realigner.set_option(<option>, <value>)

'''
        GENOTYPING
'''

job.genotyping.set_disabled(True)

# To configure the genotyping:
# job.genotyping.fast_genotyper.set_argument(<argument>, <boolean>)
#job.genotyping.unified_genotyper.set_option("--output_mode", "EMIT_ALL_CONFIDENT_SITES")

'''
        VARIANT QUALITY SCORE RECALIBRATION
'''

job.genotyping.set_perform_vqsr(False)

hapmap_resource = bina.VariantRecalibrationResource("hapmap",
    "bina://data/genome/human/hg19/dbsnp/hapmap_3.3.hg19.vcf")
hapmap_resource.set_prior(15.0)
hapmap_resource.set_known(False)
hapmap_resource.set_training(True)
hapmap_resource.set_truth(True)

omni_resource = bina.VariantRecalibrationResource("omni",
    "bina://data/genome/human/hg19/dbsnp/1000G_omni2.5.hg19.vcf")
omni_resource.set_prior(12.0)
omni_resource.set_known(False)
omni_resource.set_training(True)
omni_resource.set_truth(False)

dbsnp_resource = bina.VariantRecalibrationResource("dbsnp",
    "bina://data/genome/human/hg19/dbsnp/dbsnp_135.hg19.vcf")
omni_resource.set_prior(6.0)
omni_resource.set_known(True)
omni_resource.set_training(False)
omni_resource.set_truth(False)

# SNP Recalibration
recal_operation = bina.VariantRecalOperation()
recal_operation.set_name("ds3 SNP recalibration")
recal_operation.set_variant_type("SNP")

recal_operation.variant_recalibrator.set_option("--maxGaussians", 2)
recal_operation.variant_recalibrator.set_option("--percentBadVariants", 0.8)
recal_operation.variant_recalibrator.set_option("--use_annotation", "QD,HaplotypeScore,MQRankSum,ReadPosRankSum,FS,MQ,DP")

recal_operation.apply_recalibration.set_option("--ts_filter_level", 99.0)

recal_operation.add_resource(hapmap_resource)
recal_operation.add_resource(omni_resource)
recal_operation.add_resource(dbsnp_resource)

job.genotyping.add_recal_operation(recal_operation)

# Indel Recalibration
recal_operation = bina.VariantRecalOperation()
recal_operation.set_name("ds3 indel recalibration")
recal_operation.set_variant_type("INDEL")

recal_operation.variant_recalibrator.set_option("--maxGaussians", 2)
recal_operation.variant_recalibrator.set_option("--percentBadVariants", 0.8)
recal_operation.variant_recalibrator.set_option("--use_annotation", "QD,HaplotypeScore,FS")

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

'''
job.structural_variation.set_run_breakdancer(True)
job.structural_variation.set_run_breakseq(True)
job.structural_variation.set_run_cnvnator(True)
job.structural_variation.set_run_pindel(True)
job.structural_variation.pindel.set_use_breakdancer(True)
'''

'''
        JOB SUBMISSION
'''

# Connect to Bina Box using API key and IP addresses of nodes
binabox = bina.BinaBox()
binabox.connect(api_key, servers)

# Submit the job
job_id = binabox.run_job(job)
print "Job submitted with ID: " + job_id

# Close the connection
binabox.close()
