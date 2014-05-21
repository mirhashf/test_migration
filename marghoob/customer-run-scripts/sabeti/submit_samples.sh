#!/bin/bash

urls=
for port in 19003 ; do #19001 19005 19007 ; do
  urls="$urls http://t-rex:$port"
done

PYTHONPATH=~/river/users/marghoob/git/seqalto/loomis2/python/client/:$PYTHONPATH ../../submission/submit_jobs.py --urls $urls --datasets samples.json --enable_vqsr --enable_sv --enable_hc --keep_recalibrated --project sabeti >> jobs.csv
