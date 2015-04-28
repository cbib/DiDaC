#!/usr/bin/env python
# coding=utf8
import re

from Bio import SeqIO
from Bio.Alphabet import generic_dna
import networkx as nx

from alteration import alteration as ALT
from helpers.logger import init_logger


pattern = re.compile('([NC])_(\d+)_(\d+)')
logger = init_logger("individu graph")


class IndividuGraph:
	def __init__(self, fastq_list, k):
		self.coverage = {}
		self.alteration_list = []
		self.dbg = nx.DiGraph()
		self.dbgclean = None

		for f in fastq_list:
			fragment = pattern.search(f).group(1)
			logger.info("Considering file %s for fragment %s", f, fragment)
			comp = 0
			for record_s in SeqIO.parse(f, "fastq", generic_dna):
				sequence = str(record_s.seq)
				comp += 1
				for i2 in range(0, len(sequence) - k):
					curr_kmer = sequence[(i2):(i2 + k)]
					next_kmer = sequence[(i2 + 1):(i2 + 1 + k)]
					if next_kmer not in self.dbg:
						self.dbg.add_node(next_kmer, read_list_n=set([fragment + "_" + str(comp)]), fragment=set([fragment]))
					if curr_kmer in self.dbg:
						self.dbg.node[curr_kmer]['read_list_n'].add(fragment + "_" + str(comp))
						# if fragment not in self.dbg.node[curr_kmer]['fragment']:
						self.dbg.node[curr_kmer]['fragment'].add(fragment)
						if next_kmer not in self.dbg[curr_kmer]:
							self.dbg.add_edge(curr_kmer, next_kmer)
					else:
						self.dbg.add_node(curr_kmer, read_list_n=set([fragment + "_" + str(comp)]), fragment=set([fragment]))
						self.dbg.add_edge(curr_kmer, next_kmer)
			self.coverage[fragment] = comp


	# Compute coverage from a node for the graph .dbg 
	def total_coverage_node(self, node):
		if len(self.dbg.node[node]['fragment']) == 2:
			coverage_node = self.coverage['N'] + self.coverage['C']
		elif self.dbg.node[node]['fragment'] == "N":
			coverage_node = self.coverage['N']
		else:
			coverage_node = self.coverage['C']
		return coverage_node

	## Delete nodes of a graph G with count < coverage * alpha %
	def graph_cleaned_init(self, alpha):
		self.dbgclean = self.dbg.copy()
		nodes_count_inf_seuil = []
		for n in self.dbg:
			total_coverage = self.total_coverage_node(n)
			if len(self.dbg.node[n]['read_list_n']) <= total_coverage * alpha / 100:
				nodes_count_inf_seuil.append(n)
		self.dbgclean.remove_nodes_from(nodes_count_inf_seuil)

	# Removes edges in G_sample_test which are present in G_ref
	def graph_rmRefEdges_init(self, G2analyse, G_ref):
		self.dbg_refrm = G2analyse.copy()
		self.dbg_refrm.remove_edges_from(G_ref.edges())

	# Init list of alterations dictionnary
	def alteration_list_init(self, G_ref, k):
		self.alteration_list = []
		# Only nodes in dbg_refrm & G_ref and with in degree > 0 for end nodes and out degree > 0 for start nodes  
		shared_nodes = list(set(self.dbg_refrm.nodes()) & set(G_ref.nodes()))
		out_d = self.dbg_refrm.out_degree()
		in_d = self.dbg_refrm.in_degree()
		shared_nodes_start = [x for x in shared_nodes if out_d[x] > 0]
		shared_nodes_end = [x for x in shared_nodes if in_d[x] > 0]
		G_ref_nodes_set = set(G_ref.nodes())
		for node_start in shared_nodes_start:
			for node_end in shared_nodes_end:
				for alternative_path in nx.all_simple_paths(self.dbg_refrm, node_start, node_end):
					if len(set(alternative_path) & G_ref_nodes_set) != 2:
						continue
					# Read intersection of all nodes in the alt path for G_sample 
					read_set_pathAlt_G_sample = []
					for node in alternative_path:
						read_set_pathAlt_G_sample.append(set(self.dbg_refrm.node[node]['read_list_n']))
					intersect_allnodes_pathAlt_G_sample = set.intersection(*read_set_pathAlt_G_sample)
					# Read intersection of all nodes in the reference path for G_sample 
					for reference_path in nx.all_simple_paths(G_ref, node_start, node_end):
						condition = 0
						read_set_pathRef_G_sample = []
						for node in reference_path:
							if node not in self.dbg:
								# print ("path de référence non représenté dans GDB individu")
								condition = 1
								intersect_allnodes_pathRef_G_sample = "0"
								break
							read_set_pathRef_G_sample.append(set(self.dbg.node[node]['read_list_n']))
						if condition == 0:
							intersect_allnodes_pathRef_G_sample = set.intersection(*read_set_pathRef_G_sample)
						self.alteration_list.append(ALT(reference_path, alternative_path, len(intersect_allnodes_pathRef_G_sample), len(intersect_allnodes_pathAlt_G_sample), k))


# logger.info("start building")
# foo = IndividuGraph(['data/fastq/all_pool_trimmed0.1/C_169_1.fastq', 'data/fastq/all_pool_trimmed0.1/N_169_1.fastq'], 20)
# foo.graph_cleaned_init(3)  # .dbgclean creation
# logger.info("end building")
#
# foo.dbg.nodes(data=True)[0:15]
#

