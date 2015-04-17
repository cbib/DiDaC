#!/usr/bin/env python
# coding=utf8
import networkx as nx
import seq_lib as SL


class randomreads_graph:
	def __init__(self, cov, k):
		self.cov = cov
		read_list = SL.sampling(self.cov)
		self.dbg = nx.DiGraph()
		for i_read in range(0, len(read_list)):
			for i_kmer in range(0, len(read_list[i_read]) - k):
				curr_kmer = read_list[i_read][(i_kmer):(i_kmer + k)]
				next_kmer = read_list[i_read][(i_kmer + 1):(i_kmer + 1 + k)]
				if next_kmer not in self.dbg:
					self.dbg.add_node(next_kmer, read_list_n=[i_read])
				if curr_kmer in self.dbg:
					self.dbg.node[curr_kmer]['read_list_n'].append(i_read)
					if next_kmer not in self.dbg[curr_kmer]:
						self.dbg.add_edge(curr_kmer, next_kmer)
				else:
					self.dbg.add_node(curr_kmer, read_list_n=[i_read])
					self.dbg.add_edge(curr_kmer, next_kmer)

	def check_path(self, path):
		condition = 0
		read_set_pathAlt_G_random = []
		if path[0] not in self.dbg:
			condition = 1
			intersect_allnodes_pat_G_random = ""
		else:
			for i_node in range(0, len(path) - 1):
				if path[i_node + 1] not in self.dbg[path[i_node]]:
					condition = 1
					intersect_allnodes_pat_G_random = ""
					break
				read_set_pathAlt_G_random.append(set(self.dbg.node[path[i_node]]['read_list_n']))
		if condition == 0:
			read_set_pathAlt_G_random.append(set(self.dbg.node[path[i_node + 1]]['read_list_n']))
			# Read intersection of all nodes in the alt path for random G
			intersect_allnodes_pat_G_random = set.intersection(*read_set_pathAlt_G_random)
		return len(intersect_allnodes_pat_G_random)