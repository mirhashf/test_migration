{
    "bina_box": {
        "id": 1
    },
    "metadata": {
        "pi": "HENRY",
        "project": "BATCH SUBMIT",
        "tags": "ds9"
    },
    "workflow_type": "workflow-workgen/current/bin/workflow-workgen.jar",
    "workflow": {
        "system": {
            "workers": 1
        },
        "modules": {},
        "mounts": {
            "preload": "localhost:/mnt/pool/preload",
            "henryic": "localhost:/mnt/pool/home/henryic"
        },
        "output": {
            "uri": "mount:henryic/userdata/output"
        },
        "debug": {
            "keepScratch": "false"
        },
        "tool": {
            "name": "current/analyze/dna",
            "input": {
                "parameters": {
                    "reference": "mount:preload/ds9/b37.fa",
                    "dbsnp": "mount:preload/ds9/dbsnp_138.b37.vcf",
                    "breakpointLibrary": "mount:preload/ds9/bplib.fa"
                },
                "data": {
                    "captureRegionSize": "exome",
                    "capture": {
                        "include": {
                            "beds": [
                                "mount:preload/ds9/ds9.exome.bed"
                            ]
                        }
                    },
                    "pairs": [
                        {
                            "identifier": "<samples_list.readname>",
                            "sample": "<samples_list.samplename>",
                            "library": "<samples_list.library>",
                            "platform": "Illumina",
                            "reads": [
                                "<samples_list.end0_path>",
                                "<samples_list.end1_path>"
                            ]
                        }
                    ]
                },
                "configuration": {
                    "qc": {
                        "sequencing": {
                            "select": "ui-true",
                            "ui-true": {
                                "enabled": true
                            },
                            "ui-false": {
                                "enabled": false
                            }
                        },
                        "alignment": {
                            "select": "ui-true",
                            "ui-true": {
                                "enabled": true
                            },
                            "ui-false": {
                                "enabled": false
                            }
                        },
                        "variantCalling": {
                            "select": "ui-false",
                            "ui-true": {
                                "enabled": true
                            },
                            "ui-false": {
                                "enabled": false
                            }
                        }
                    },
                    "alignment": {
                        "output": false,
                        "unmapped": false,
                        "aligner": {
                            "select": "ui-bwamem",
                            "ui-none": {
                                "tool": null
                            },
                            "ui-bwamem": {
                                "tool": "bwamem"
                            },
                            "ui-bina": {
                                "tool": "bina",
                                "options": {
                                    "align": {
                                        "-trim": 30
                                    }
                                }
                            }
                        },
                        "sorter": {
                            "select": "ui-markduplicates",
                            "ui-markduplicates": {
                                "options": {
                                    "-mark_duplicates": true
                                }
                            },
                            "ui-base": {}
                        }
                    },
                    "reads": {
                        "realignment": {
                            "output": false
                        },
                        "recalibration": {
                            "output": false,
                            "options": {
                                "baseRecalibrator": {
                                    "--maximum_cycle_value": 500,
                                    "--sort_by_all_columns": "",
                                    "-cov": [
                                        "ReadGroupCovariate",
                                        "QualityScoreCovariate",
                                        "CycleCovariate",
                                        "ContextCovariate"
                                    ]
                                },
                                "printReads": {
                                    "-baq": "RECALCULATE"
                                }
                            }
                        }
                    },
                    "variants-simple": {
                        "call": {
                            "select": "ui-none",
                            "ui-hc": {
                                "tool": "haplotypecaller",
                                "enabled": true,
                                "options": {
                                    "-pairHMM": "VECTOR_LOGLESS_CACHING",
                                    "--annotation": [
                                        "AlleleBalance",
                                        "Coverage",
                                        "MappingQualityZero"
                                    ],
                                    "--standard_min_confidence_threshold_for_emitting": 0.1,
                                    "--variant_index_type": "LINEAR",
                                    "--variant_index_parameter": 128000
                                }
                            },
                            "ui-gvcf": {
                                "tool": "haplotypecaller",
                                "enabled": true,
                                "options": {
                                    "-pairHMM": "VECTOR_LOGLESS_CACHING",
                                    "--annotation": [
                                        "AlleleBalance",
                                        "Coverage",
                                        "MappingQualityZero"
                                    ],
                                    "--standard_min_confidence_threshold_for_emitting": 0.1,
                                    "--variant_index_type": "LINEAR",
                                    "--variant_index_parameter": 128000,
                                    "--emitRefConfidence": "GVCF"
                                }
                            },
                            "ui-none": {
                                "tool": "haplotypecaller",
                                "enabled": false
                            }
                        },
                        "vqsr": {
                            "select": "ui",
                            "ui": {
                                "indel": {
                                    "train": {
                                        "options": {
                                            "--badLodCutoff": 10,
                                            "--use_annotation": [
                                                "FS",
                                                "MQRankSum",
                                                "ReadPosRankSum"
                                            ],
                                            "--maxGaussians": 1
                                        },
                                        "resources": [
                                            {
                                                "name": "dbsnp",
                                                "known": true,
                                                "training": false,
                                                "truth": false,
                                                "prior": 2,
                                                "file": "mount:preload/ds9/dbsnp_138.b37.vcf"
                                            },
                                            {
                                                "name": "mills",
                                                "known": false,
                                                "training": true,
                                                "truth": true,
                                                "prior": 12,
                                                "file": "mount:preload/ds9/Mills_and_1000G_gold_standard.indels.b37.vcf"
                                            }
                                        ]
                                    },
                                    "apply": {
                                        "--ts_filter_level": 99.9
                                    },
                                    "enabled": false
                                },
                                "snp": {
                                    "train": {
                                        "options": {
                                            "--badLodCutoff": 10,
                                            "--use_annotation": [
                                                "QD",
                                                "FS",
                                                "MQRankSum",
                                                "ReadPosRankSum"
                                            ],
                                            "--maxGaussians": 1
                                        },
                                        "resources": [
                                            {
                                                "name": "dbsnp",
                                                "known": true,
                                                "training": false,
                                                "truth": false,
                                                "prior": 2,
                                                "file": "mount:preload/ds9/dbsnp_138.b37.vcf"
                                            },
                                            {
                                                "name": "hapmap",
                                                "known": false,
                                                "training": true,
                                                "truth": true,
                                                "prior": 15,
                                                "file": "mount:preload/ds9/hapmap_3.3.b37.vcf"
                                            },
                                            {
                                                "name": "omni",
                                                "known": false,
                                                "training": true,
                                                "truth": true,
                                                "prior": 12,
                                                "file": "mount:preload/ds9/1000G_omni2.5.b37.vcf"
                                            },
                                            {
                                                "name": "1000g",
                                                "known": false,
                                                "training": true,
                                                "truth": false,
                                                "prior": 10,
                                                "file": "mount:preload/ds9/1000G_phase1.snps.high_confidence.b37.vcf"
                                            }
                                        ]
                                    },
                                    "apply": {
                                        "--ts_filter_level": 99.9
                                    },
                                    "enabled": false
                                }
                            }
                        }
                    },
                    "variants-structural": {
                        "cnvnator": {
                            "select": "ui-false",
                            "ui-true": {
                                "enabled": true
                            },
                            "ui-false": {
                                "enabled": false
                            }
                        },
                        "breakdancer": {
                            "select": "ui-false",
                            "ui-true": {
                                "enabled": true,
                                "options": {
                                    "generate": {
                                        "-m": 2000000
                                    }
                                }
                            },
                            "ui-false": {
                                "enabled": false
                            }
                        },
                        "pindel": {
                            "select": "ui-false",
                            "ui-true": {
                                "enabled": true,
                                "options": {
                                    "generate": {
                                        "-T": 4
                                    }
                                }
                            },
                            "ui-false": {
                                "enabled": false
                            }
                        },
                        "breakseq": {
                            "select": "ui-false",
                            "ui-true": {
                                "enabled": true
                            },
                            "ui-false": {
                                "enabled": false
                            }
                        },
                        "metasv": {
                            "select": "ui-false",
                            "ui-true": {
                                "enabled": true
                            },
                            "ui-false": {
                                "enabled": false
                            }
                        }
                    }
                }
            }
        }
    }
}
