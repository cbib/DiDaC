#!/usr/bin/env python
# coding=utf8
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import networkx as nx 
import re
from alteration import alteration as ALT

pattern = re.compile('([NC])_(\d+)_(\d+)')

class individu_graph:
	def __init__ (self,fastq_list,k):
		self.coverage = {}
		dbg=nx.DiGraph()
		for f in fastq_list: 
			fragment = pattern.search(f).group(1)
			comp = 0
			#reads_list = []
			for record_s in SeqIO.parse(f, "fastq", generic_dna):
				sequence = str(record_s.seq)
				comp += 1
				for i2 in range(0,len(sequence)-k):
					curr_kmer=sequence[(i2):(i2+k)]
					next_kmer=sequence[(i2+1):(i2+1+k)]
					if next_kmer not in dbg:
				 		dbg.add_node(next_kmer,reads_list_n=[fragment+"_"+str(comp)],fragment=[fragment])
					if curr_kmer in dbg:
						dbg.node[curr_kmer]['reads_list_n'].append(fragment+"_"+str(comp))
						if fragment not in dbg.node[curr_kmer]['fragment']:
							dbg.node[curr_kmer]['fragment'].append(fragment)
						if next_kmer not in dbg[curr_kmer]:
				  			dbg.add_edge(curr_kmer,next_kmer)
					else:
						dbg.add_node(curr_kmer,reads_list_n=[fragment+"_"+str(comp)],fragment=[fragment])
						dbg.add_edge(curr_kmer,next_kmer)
			self.coverage[fragment] = comp
		self.dbg = dbg
		self.alteration_list = []

	# Compute coverage from a node for the graph .dbg 
	def total_coverage_node(self,node):
		coverage_node = 0
		N,C = 0,0
		if len(self.dbg.node[node]['fragment']) == 2:
			coverage_node = self.coverage['N']+self.coverage['C']
		elif self.dbg.node[node]['fragment'] == "N":
			coverage_node = self.coverage['N']
		else:
			coverage_node = self.coverage['C']
		return(coverage_node)

	## Delete nodes of a graph G with count < coverage * X% 
	def graph_cleaned_init(self,X):
		self.dbgclean = self.dbg.copy()
		nodes_count_inf_seuil = []
		for n in self.dbg:
			total_coverage = self.total_coverage_node(n)
			if len(set(self.dbg.node[n]['reads_list_n'])) <= total_coverage*X/100:
				nodes_count_inf_seuil.append(n)
		self.dbgclean.remove_nodes_from(nodes_count_inf_seuil)

	# Removes edges in G_sample_test which are present in G_ref
	def graph_rmRefEdges_init(self,G2analyse,G_ref):
		self.dbg_refrm = G2analyse.copy()
		self.dbg_refrm.remove_edges_from(G_ref.edges())

# compute path differences dictionnary
	def alteration_list_init(self,G_ref,k):
		self.alteration_list = []
		# Only nodes in dbg_refrm & G_ref and with in degree > 0 for end nodes and out degree > 0 for start nodes  
		shared_nodes = list(set(self.dbg_refrm.nodes()) & set(G_ref.nodes()))
		out_d=self.dbg_refrm.out_degree()
		in_d=self.dbg_refrm.in_degree()
		shared_nodes_start = [x for x in shared_nodes if out_d[x]>0]
		shared_nodes_end = [x for x in shared_nodes if in_d[x]>0]
		G_ref_nodes_set = set(G_ref.nodes())
		for node_start in shared_nodes_start:
			for node_end in shared_nodes_end:
				for path in nx.all_simple_paths(self.dbg_refrm,node_start,node_end):
					# Just to print exeptions 
					if len(set(path) & G_ref_nodes_set) != 2:
						print (len(set(path) & G_ref_nodes_set))
						print (set(path) & G_ref_nodes_set)
						for path_ref in nx.all_simple_paths(G_ref,node_start,node_end):	
							print(path_ref)
							print(ALT.kmerpathToSeq(path_ref,k))
						print(path)
						print(ALT.kmerpathToSeq(path,k))
						continue

					# Read intersection of all nodes in the alt path for G_sample 
					read_set_pathAlt_G_sample = []
					for node in path:
				 		read_set_pathAlt_G_sample.append(set(self.dbg_refrm.node[node]['reads_list_n']))
					intersect_allnodes_pathAlt_G_sample = set.intersection(*read_set_pathAlt_G_sample)
					
					for path_ref in nx.all_simple_paths(G_ref,node_start,node_end):						
						self.alteration_list.append(ALT(path_ref,path,len(intersect_allnodes_pathAlt_G_sample),k))
						
