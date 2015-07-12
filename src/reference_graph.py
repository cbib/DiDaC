#!/usr/bin/env python
# coding=utf8
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import networkx as nx

####################################################
####################################################
## Parameters loading 
# Data directories
data_directory = "data/ref"
VE_anno_file = data_directory + "/VE_anno.tab"
dbsnp_list_file = data_directory + "/dbsnp.tab"
dbsnp_anno_file = data_directory + "/dbsnp_anno.tab"
fasta_list = [data_directory + "/p53var1.fasta", data_directory + "/p53var3.fasta"] #, data_directory + "/p53var4.fasta"
end_3prime_utr = 202

####################################################
####################################################
## Loading SNP position and alt for each VE
# Return graph has argv : 
# node : ref_list = dict : key = ARNm reference ; value = start position of the kmer in the reference sequence

def ref_constructor(k):
	VE_dict = {}
	all_records = []
	Anno_hash = {}
	startPosition = {}
	dbg_ref = nx.DiGraph()
	IN_VELIST = open(VE_anno_file, 'r')
	lines = IN_VELIST.readlines()
	lines = map(str.strip, lines)
	# Record of each VE: name, ID, 3' & 5' primer position
	for l in lines:
		VE_dict[l.split("\t")[1]] = {'name': l.split("\t")[0], 'N_fgmt': l.split("\t")[2], 'C_fgmt': l.split("\t")[3]}
	IN_LISTdbSNP = open(dbsnp_list_file, 'r')
	# Record of SNP to include 
	dbSNPList = IN_LISTdbSNP.readlines()
	dbSNPList = map(str.strip, dbSNPList)
	# For each SNP, record of the position (same for the 3 VE)
	IN_ANNOdbSNP = open(dbsnp_anno_file, 'r')
	while 1:
		line_ = IN_ANNOdbSNP.readline()
		if line_ == '':
			break
		line_ = line_.rstrip()
		line_split = line_.split("\t")
		if line_split[0] in dbSNPList:
			if line_split[4] == "NM_000546":
				Anno_hash[line_split[0]] = {'pos': line_split[7], 'alt': Seq(line_split[1].split("|")[1], generic_dna)}
	# Record of all sequences to include in the Reference De Bruijn Graph 
	for f in fasta_list:
		for record in SeqIO.parse(f, "fasta", generic_dna):
			record.name = VE_dict[record.id]['name']
			# For each SNP, record of the 2 k-mers around it for all VE
			if record.name == "NM_000546.5":
				for rs in Anno_hash:
					alt_pos = int(Anno_hash[rs]['pos'])
					kmer_around = str(record.seq[alt_pos - k - 1:alt_pos - 1] + Anno_hash[rs]['alt'].complement() + record.seq[alt_pos:alt_pos + k])
					all_records.append(SeqRecord(Seq(kmer_around, generic_dna), id=rs, name=rs))
					startPosition[rs] = alt_pos - k - end_3prime_utr
			fragment = str(record.seq[int(VE_dict[record.id]['N_fgmt'].split(":")[0]) - 1:int(VE_dict[record.id]['N_fgmt'].split(":")[1]) - 1])
			all_records.append(SeqRecord(Seq(fragment, generic_dna), id=record.name + "_" + 'N', name=record.name))
			startPosition[record.name + "_" + 'N'] = int(VE_dict[record.id]['N_fgmt'].split(":")[0]) - end_3prime_utr
			fragment = str(record.seq[int(VE_dict[record.id]['C_fgmt'].split(":")[0]) - 1:int(VE_dict[record.id]['C_fgmt'].split(":")[1]) - 1])
			all_records.append(SeqRecord(Seq(fragment, generic_dna), id=record.name + "_" + 'C', name=record.name))
			startPosition[record.name + "_" + 'C'] = int(VE_dict[record.id]['C_fgmt'].split(":")[0]) - end_3prime_utr
	# Graph construction
	for i in range(0, len(all_records)):
		seq_s = str(all_records[i].seq)
		for i2 in range(0, len(seq_s) - k):
			curr_kmer = seq_s[(i2):(i2 + k)]
			next_kmer = seq_s[(i2 + 1):(i2 + 1 + k)]
			if next_kmer not in dbg_ref:
				dbg_ref.add_node(next_kmer, ref_list={all_records[i].name: i2 + startPosition[all_records[i].id]})
			else:
				dbg_ref.node[next_kmer]['ref_list'][all_records[i].name] = i2 + startPosition[all_records[i].id]	
			if curr_kmer in dbg_ref:
				dbg_ref.node[curr_kmer]['ref_list'][all_records[i].name] = i2 + startPosition[all_records[i].id]
				if dbg_ref[curr_kmer].get(next_kmer, 0) == 0:
					dbg_ref.add_edge(curr_kmer, next_kmer)
			else:
				dbg_ref.add_node(curr_kmer, ref_list={all_records[i].name: i2 + startPosition[all_records[i].id]})
				dbg_ref.add_edge(curr_kmer, next_kmer)
	return dbg_ref