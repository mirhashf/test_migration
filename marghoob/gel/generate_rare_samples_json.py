#!/usr/bin/python

from glob import glob
import os
import sys
import pysam
import json

raw_path="/net/kodiak/volumes/delta/shared/prj/GELA20140221/rare/raw/"
delta_path="/net/kodiak/volumes/delta/shared/"

datasets = {}
#for sample_path in glob('%s/*/*/*' % (raw_path)):
for line in sys.stdin.readlines():
  sample_path = os.path.normpath(line.strip())
  sys.stderr.write("%s\n" % (sample_path))
  alignment_readgroups = []
  skip=False
  for bam in glob('%s/*.rg.bam' % (sample_path)):
    sys.stderr.write("Examining %s\n" % (bam))
    fields = bam.split(".")
    bam_rg = fields[1]
    bam_file = pysam.Samfile(bam, 'rb')
    readgroups = bam_file.header['RG']
    readgroup_found = False
    for readgroup in readgroups:
      if readgroup['ID'] == bam_rg:
        rg = bam_rg
        sample = readgroup['SM']
        library = readgroup['LB'] if 'LB' in readgroup else readgroup['PU']
        platform = readgroup['PL']
        readgroup_found = True
        break
    if not readgroup_found:
      sys.stderr.write("%s will need to be re-headered. it is missing readgroup %s\n" % (bam, bam_rg))
      skip = True
    alignment_readgroup = {"read_group": rg, "platform": platform, "sample": sample, "platform": platform, "library": library, "bam": "delta:/%s" % (os.path.relpath(bam, delta_path))}
    alignment_readgroups.append(alignment_readgroup)
    bam_file.close()
    if skip: break
  if skip: continue
  sample_path_split = sample_path.split("/")
  output_prefix = "/net/kodiak/volumes/delta/shared/prj/GELA20140221/rare/runs/single/%s/%s" % (sample_path_split[-3], sample_path_split[-2])
  if not os.path.exists(output_prefix): sys.stderr.write("%s non-existent\n" % (output_prefix))
  dataset = {"alignment_groups": alignment_readgroups, "output_prefix": "delta:/%s" % (os.path.relpath(output_prefix, delta_path))}
  datasets[sample_path_split[-2]] = dataset
  #print alignment_readgroups
  #print json.dumps(alignment_readgroups, indent=4)
  #sys.exit(0)

print json.dumps(datasets, indent=4)
