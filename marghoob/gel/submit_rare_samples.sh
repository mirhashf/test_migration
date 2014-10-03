#!/bin/bash

urls=
for port in 19001 19003 19005 19007 19011; do
  urls="$urls http://t-rex:$port"
done

PYTHONPATH=~/river/users/marghoob/git/seqalto/loomis2/python/client/:$PYTHONPATH ../submission/submit_jobs.py --urls $urls --datasets rare_samples.json --enable_vqsr --enable_sv --enable_hc --aligners bwamem --keep_recalibrated --project gel >> jobs.csv
