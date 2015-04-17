#!/usr/bin/env python
# coding=utf8

import os
import re
import sys


# Parcours de l'ensemble de mes chemins --> A mettre en fonction!
OUT1 = open("/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/20150413/dbg/data/seq/NwP.seq",'w')
OUT2 = open("/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/20150413/dbg/data/seq/CwP.seq",'w')
pattern = re.compile('([NC])_(\d+)_(\d+)')

FASTQFILE_ALL = os.listdir("/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/20150413/primer_cutting/fastq/all_pool/fastq")
for i in range(0,len(FASTQFILE_ALL)):
	if FASTQFILE_ALL[i] == ".DS_Store":
		continue
	fgmt = pattern.search(FASTQFILE_ALL[i]).group(1)
	fastq = open("/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/20150413/primer_cutting/fastq/all_pool/fastq/"+FASTQFILE_ALL[i],'r')
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