#!/usr/bin/env python

# Simple script to submit a multisample job
import sys
from client import Client
import json
import ConfigParser

AVAIL_BOXES = [3, 6, 9]

if len(sys.argv) != 3:
    print 'Need to specify box id, sample name.'
    exit(0)
else:
    box_id = int(sys.argv[1])
    sample_name = sys.argv[2]

if box_id not in AVAIL_BOXES:
    print 'Invalid box specified ' + str(AVAIL_BOXES) + ': ' + str(box_id)
    exit(1)

if box_id in [6, 7]:
    hw_rev = 'A00'
elif box_id in [3, 9]:
    hw_rev = 'A01'
else:
    print 'The hardware rev mapper is confused'
    exit(1)

if hw_rev == 'A00':
    worker_num = 5
    max_mem = 60
    thread_num = 16
    bwa_threads = 14
elif hw_rev == 'A01':
    worker_num = 4
    max_mem = 85
    thread_num = 24
    bwa_threads = 18
else:
    print 'The hardware rev parameter mapper is confused'
    exit(1)

config = ConfigParser.ConfigParser();
config.read("config.cfg")

# Account credentials
baseurl = config.get("bina", "api-url")
username = config.get("bina", "username")
password = config.get("bina", "password")
binabox_id = config.get("bina", "box-id")

# define job (could use a for loop here)
job_obj = {
    "workflow": {
        "alignment_groups": [
            {
                "read_group": "RG0",
                "sample": sample_name,
                "library": "L0",
                "platform": "Illumina",
                "bam": "gsfs0:/VA/bam/" + sample_name + ".bam"
            }
        ],
        "fasta": "gsfs0:/referencefiles/ucsc.hg19.fa",
        "dbsnp": "gsfs0:/referencefiles/dbsnp_135.hg19.vcf",
        "output_prefix": "gsfs0:/output/VA/reproc-out/" + sample_name,
        "gatk_path": "gsfs0:/VA/gatk/gatk-2.7-4",
        "run_breakdancer": True,
        "run_cnvnator": True,
        "run_pindel": True,
        "run_breakseq": True,
        "bplib": "gsfs0:/referencefiles/bplib.fa",
        "keep_sorted_bams": False,
        "keep_realigned_bams": False,
        "keep_recalibrated_bams": True,
        "worker_num": worker_num,
        "run_bwa": True,
        "run_reduce_reads": True,
        "run_vqsr_snp": True,
        "run_vqsr_indel": True,
        "keep_scratch": False,
        "disable_variant_calling": False,
        "sorter_args": {
            "-mark_duplicates": "",
            "-max_mem_gb": max_mem
        },
        "mounts": [
            {
                "name": "gs1",
                "address": "scg-cnfs.sunet:/srv/gs1/projects/scg/Bina"
            },
            {
                "name": "gsfs0",
                "address": "scg-gs-cnfs.sunet:/srv/gsfs0/projects/gbsc/Bina"
            }
        ],
        "thread_num": thread_num,
        "bina_aligner_args": {
            "-p": 24,
            "-trim": 30
        },
        "genotyper_args": {
            "-stand_emit_conf": 0.1,
            "-A": [
                "AlleleBalance",
                "Coverage",
                "MappingQualityZero"
            ],
            "-baq": "CALCULATE_AS_NECESSARY",
            "--genotype_likelihoods_model": "BOTH",
            "--output_mode": "EMIT_ALL_CONFIDENT_SITES"
        },
        "reduce_reads_args": {
            "-dcov": 25
        },
        "bwa_mem_args": {
            "-t": 24
        },
        "bwa_aln_args": {
            "-t": 4,
            "-q": 30
        },
        "bwa_sampe_args": {
            "-P": ""
        },
        "bwarunner_args": {
            "--nthreads": bwa_threads,
            "--splitter_blocksize": 65536,
            "--splitter_max_pending": 4
        },
        "base_recalibrator_args": {
            "--disable_indel_quals": "",
            "-dcov": 250,
            "-cov": [
                "ReadGroupCovariate",
                "QualityScoreCovariate",
                "CycleCovariate",
                "ContextCovariate"
            ]
        },
        "count_covariates_args": {
            "-cov": [
                "ReadGroupCovariate",
                "QualityScoreCovariate",
                "CycleCovariate",
                "DinucCovariate"
            ]
        },
        "table_recal_args": {
            "-baq": "RECALCULATE",
            "--doNotWriteOriginalQuals": ""
        },
        "vqsr_snp_train_args": {
            "--numBadVariants": 3000
        },
        "vqsr_snp_apply_args": {
            "--ts_filter_level": 99.9
        },
        "vqsr_indel_train_args": {
            "--numBadVariants": 3000,
            "--maxGaussians": 4
        },
        "vqsr_indel_apply_args": {
            "--ts_filter_level": 99.9
        },
        "vqsr_groups": {
            "snp": {
                "annotations": [
                    "QD",
                    "MQRankSum",
                    "ReadPosRankSum",
                    "FS",
                    "DP",
                    "HaplotypeScore"
                ],
                "resources": [
                    "hapmap,VCF,known=false,training=true,truth=true,prior=15.0 gsfs0:referencefiles/hapmap_3.3.hg19.sites.vcf",
                    "omni,VCF,known=false,training=true,truth=true,prior=12.0   gsfs0:referencefiles/1000G_omni2.5.hg19.sites.vcf",
                    "dbsnp,VCF,known=true,training=false,truth=false,prior=2.0  gsfs0:referencefiles/dbsnp_135.hg19.vcf"
                ]
            },
            "indel": {
                "annotations": [
                    "DP",
                    "FS",
                    "MQRankSum",
                    "ReadPosRankSum"
                ],
                "resources": [
                    "mills,VCF,known=false,training=true,truth=true,prior=12.0 gsfs0:referencefiles/Mills_and_1000G_gold_standard.indels.hg19.sites.rightHeader.vcf",
                    "dbsnp,VCF,known=true,training=false,truth=false,prior=2.0 gsfs0:referencefiles/dbsnp_135.hg19.vcf"
                ]
            }
        },
        "pindel_args": {
            "-T": thread_num/2
        },
        "uid": 925,
        "gid": 1662
    },
    "bina_box": {
        "id": box_id
    },
    "metadata": {
        "tags": sample_name + ",bwa,gatk2.7-4,rread,refcall,cnv,bd,pd,bs",
        "project": "AAA Reprocess",
        "library": "L0",
        "pi": "Phil Tsao",
        "sample": sample_name
    },
    "workflow_type": "wgs-workflow/current/bin/wgs-workflow.jar"
}

# Main code
c = Client(baseurl, username, password, False)
c.login()
job = c.post(job_obj, "job_list")
print "Job submitted with ID: " + job["id"]
c.logout()
