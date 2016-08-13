#!/usr/bin/python

# Goes into the 20min job from the post request and digs out the workflow and turns off checkpointing
# Example:
# 	- Open network tab and copy body of the job saving request
#	- run `pbpaste | ./format_20min_job.py | pbcopy
#	- The properly formatted job should be in your clipboard now

import json
import sys

def getPredicate(gmsObject, predicate):
	for temp in gmsObject['properties']:
		if temp['predicate'] == predicate:
			return temp['object']

def getScript(gmsObject):
	return getPredicate(gmsObject, 'gms:code/workflow/script')

def getParameterJson(gmsObject):
	return getPredicate(gmsObject, 'gms:data/parameterset/json')
				
if len(sys.argv) == 2:
	job = getScript(json.load(open(sys.argv[1])))
elif len(sys.argv) == 1:
	job=getScript(json.load(sys.stdin))
else:
	print('get lost')
		
presets = job['workflow']['tool']['input']['presets']
for preset in presets:
	script = getScript(presets[preset])
	params = getParameterJson(script['workflow']['tool']['input']['configuration']['analyze']['options']['configuration'])
	for param in params:
		if param['option'] == 'params.runCheckpoints':
			param['value'] = False
	

print(json.dumps(job))
