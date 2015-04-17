#!/usr/bin/env python
# coding=utf8

import os
import re
import sys

rate_list = [0.1,0.15,0.2]

for rate in rate_list:
	OUT1 = open("data/seq/N_454_"+str(rate)+".seq",'w')
	OUT2 = open("data/seq/C_454_"+str(rate)+".seq",'w')
	pattern = re.compile('([NC])_(\d+)_(\d+)')

	FASTQFILE_ALL = os.listdir("data/fastq/454/fastq_trimmed"+str(rate))
	for i in range(0,len(FASTQFILE_ALL)):
		if FASTQFILE_ALL[i] == ".DS_Store":
			continue
		fgmt = pattern.search(FASTQFILE_ALL[i]).group(1)
		fastq = open("data/fastq/454/fastq_trimmed"+str(rate)+"/"+FASTQFILE_ALL[i],'r')
		lignes  = fastq.readlines()
		fastq.close()
		if fgmt == "N":
			for i_lignes in range(1,len(lignes),4):
				OUT1.write(lignes[i_lignes])
		else:
			for i_lignes in range(1,len(lignes),4):
				OUT2.write(lignes[i_lignes])

	OUT1.close()
	OUT2.close()