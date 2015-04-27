#!/usr/bin/env python
# coding=utf8

import os
import re
import sys
import random


# For an allowed mismatch rate during cutting = 0.1
FASTQFILE_PATH = "/data/fastq/all_pool_trimmed0.1"

pattern = re.compile('([NC])_(\d+)_(\d+)')
read_library = {'N':{},'C':{}}

FASTQFILE_ALL = os.listdir(FASTQFILE_PATH)
for i in range(0,len(FASTQFILE_ALL)):
	if FASTQFILE_ALL[i] == ".DS_Store":
		continue
	fragment = pattern.search(FASTQFILE_ALL[i]).group(1)
	individu = pattern.search(FASTQFILE_ALL[i]).group(2)
	fastq = open(FASTQFILE_PATH+"/"+FASTQFILE_ALL[i],'r')
	lines  = fastq.readlines()
	fastq.close()
	lines = map(str.strip, lines)
	if individu not in read_library[fragment]:
		read_library[fragment][individu] = []
	for i_line in range(1,len(lines),4):
		read_library[fragment][individu].append(lines[i_line])

# Sampling fonction from a coverage dict (with keys 'N' et 'C') 
def sampling(coverage_dict):
	read_sampling = []
	for fragment,coverage in coverage_dict.items():
		read_sampling_for_remainder = []
		read_number_by_sample = coverage/len(read_library[fragment])
		read_number_remainder = coverage % len(read_library[fragment])		
		for sample in read_library[fragment]:
			read_sampling += random.sample(read_library[fragment][sample],read_number_by_sample)
			read_sampling_for_remainder += random.sample(read_library[fragment][sample],1)
		read_sampling += random.sample(read_sampling_for_remainder,read_number_remainder)
	return read_sampling
