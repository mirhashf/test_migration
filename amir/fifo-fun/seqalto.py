#!/usr/bin/python
import csv
import sys
import os
import subprocess
from threading import Thread
import re
from config import Config

def samtools(chr_fifo, chr_name):
	print "Starting SAMtools thread on: " + chr_fifo
	os.system("samtools view -ubS  " +  chr_fifo + " | " + " samtools sort /dev/stdin "+ chr_name + "_sorted")

def main():
	config_file = file(sys.argv[1])
	conf = Config(config_file)
	output_dir = conf.output_dir
	output_dict = {}
	merged_memory_header = ""
	header_written = False

	read_group_string = ""
	#prepare header readgroup tags
	for job in conf.jobs:
		read_group_string += "@RG\tID:"+ job.read_group + "\tLB:" + job.library + "\tSM:" + job.sample + "\tPL:" + job.platform +  "\n"

	for job in conf["jobs"]:
		aligner = subprocess.Popen(['seqalto','align', '-1', job["first_end"], '-2', job["second_end"], '-p', conf["num_threads"], '--rg', job["read_group"], '--template_len_comp_method', "2", "--trim" , conf.trim_q], stdout=subprocess.PIPE)
		for line in aligner.stdout:
			row = line.split("\t",3)
			if row[0][0] != '@' and row[2] != '*' and row[2] in conf['references']:
				if header_written == False:
					merged_memory_header += read_group_string
				header_written = True
				if row[2] not in output_dict:
					chr_fifo = output_dir + "/" + row[2]
					os.mkfifo(chr_fifo)
					t = Thread(target=samtools, args=(chr_fifo,row[2]))
					t.start()
					output_dict[row[2]] = open(output_dir + "/" + row[2], "w")
 	                                output_dict[row[2]].write(merged_memory_header)
				output_dict[row[2]].write(line)
			elif row[0][0] == '@' and header_written != True:
				if row[0][0:3]!= '@RG' and row[0][0:3]!= '@PG':
					merged_memory_header += line

	for fifo in output_dict.keys():
		os.remove(fifo)

if __name__ == "__main__":
	main()
