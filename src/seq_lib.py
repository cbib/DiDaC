#!/usr/bin/env python
# coding=utf8
## random sampling in C.seq and N.seq files which are PacBio sequences files
import random

# Data loading and construction
reads_seq = {}

fragment_files = {'N': "data/seq/N_0.1.seq", 'C': "data/seq/C_0.1.seq"}
for fragment, path in fragment_files.items():
	reads_seq[fragment] = []
	with open(path, 'r') as sequence_file:
		for line in sequence_file:
			line = line.rstrip()
			reads_seq[fragment].append(line)

# Sampling fonction from a coverage dict (with keys 'N' et 'C') 
def sampling(cov_dict):
	readssampling = []
	for fragment, coverage in cov_dict.items():
		readssampling += random.sample(reads_seq[fragment], coverage)
	return readssampling

