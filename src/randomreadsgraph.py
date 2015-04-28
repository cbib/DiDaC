#!/usr/bin/env python
# coding=utf8
import collections
import networkx as nx
from helpers.logger import init_logger
import seq_lib as SL


logger = init_logger("RANDGRAPH")


cached_kmers={}
cached_pairs={}
KMER_UID={}
curr_uid=0
last_sample=None

class RandomReadsGraph:
	def __init__(self, coverage_dict, k,restrict_to={}):
		global cached_kmers,curr_uid,last_sample
		self.coverage = coverage_dict
		self.kmer_map=collections.defaultdict(set)
		self.restrict_to=set(restrict_to)
		self.possible_pairs=set()
		read_list = SL.sampling(self.coverage)

		self.dbg = nx.DiGraph()
		# logger.info("Will process %d reads",len(read_list))
		for i_read in range(0, len(read_list)):
		# for i_read in range(0, len(read_list[0:1000])):

			this_read=read_list[i_read]
			# if this_read not in cached_kmers:
			# 	these_kmers=[this_read[i:i+k] for i in xrange(len(this_read)-k)]
			# 	for km in these_kmers:
			# 		if km not in KMER_UID:
			# 			KMER_UID[km]=curr_uid
			# 			curr_uid+=1
			# 	cached_kmers[this_read]=map(lambda x:KMER_UID[x],these_kmers)
			# 	cached_pairs[this_read]=[(f1,f2) for f1,f2 in zip(cached_kmers[this_read],cached_kmers[this_read][1:])]
			#
			# kkmers=cached_kmers[this_read]
			# kmers_pairs=cached_pairs[this_read]

			kkmers=[this_read[i:i+k] for i in xrange(len(this_read)-k) if this_read[i:i+k] in self.restrict_to]
			kmers_pairs=[(f1,f2) for f1,f2 in zip(kkmers,kkmers[1:])]

			for kmer in kkmers:
				self.kmer_map[kmer].add(i_read)

			self.possible_pairs.update(kmers_pairs)

			# self.dbg.add_path(kkmers)
			# for kmer in kkmers:
			# 	if 'read_list_n' not in self.dbg.node[kmer]:
			# 		self.dbg.node[kmer]['read_list_n']=set()
			# 	self.dbg.node[kmer]['read_list_n'].add(i_read)

			# for i_kmer in range(0, len(read_list[i_read]) - k):
			#
			# 	curr_kmer = read_list[i_read][(i_kmer):(i_kmer + k)]
			# 	next_kmer = read_list[i_read][(i_kmer + 1):(i_kmer + 1 + k)]
			#
			# 	if curr_kmer not in self.dbg:
			# 		self.dbg.add_node(curr_kmer, read_list_n={i_read})
			# 	else:
			# 		self.dbg.node[curr_kmer]['read_list_n'].add(i_read)
			#
			# 	if next_kmer not in self.dbg:
			# 		self.dbg.add_node(next_kmer, read_list_n={i_read})
			# 	else:
			# 		self.dbg.node[next_kmer]['read_list_n'].add(i_read)
			#
			# 	self.dbg.add_edge(curr_kmer, next_kmer)


	def build_read_set_for_path(self, a_path, verbose=False):
		# a_path=map(lambda x:KMER_UID[x],a_path)
		missing_kmers = set(a_path).difference(self.kmer_map)
		if len(missing_kmers):
			# logger.critical("Completely missing kmer (%d): %s", len(missing_kmers), missing_kmers)
			return set()

		# current_set = set(self.kmer_map[a_path[0]]['read_list_n'])
		current_set = set(self.kmer_map[a_path[0]])
		if verbose:
			print len(current_set)

		assert isinstance(current_set, set)

		for i, a_node in enumerate(a_path[1:]):
			# if not self.dbg.has_edge(a_path[i], a_node):
			# 	logger.critical("Missing edge between %s -> %s",a_path[i],a_node)

			if not (a_path[i],a_node) in self.possible_pairs:
				logger.critical("Missing edge between %s -> %s",a_path[i],a_node)

			# current_set.intersection_update(self.dbg.node[a_node]['read_list_n'])
			current_set.intersection_update(self.kmer_map[a_node])
			if verbose:
				print len(current_set)
		return current_set

	def check_path(self, reference_path, alternative_path):
		ref_path_read_set = self.build_read_set_for_path(reference_path)
		alt_path_read_set = self.build_read_set_for_path(alternative_path)
		if len(alt_path_read_set) == 0:
			logger.critical("Empty ALT read set for path %s", alternative_path)
		if len(ref_path_read_set) == 0:
			logger.critical("Empty REF read set for path %s", reference_path)
		if len(alt_path_read_set) + len(ref_path_read_set) == 0:
			ratio = float('+inf')
		else:
			ratio = float(len(alt_path_read_set)) / (len(alt_path_read_set) + len(ref_path_read_set))
		return ratio, len(ref_path_read_set), len(alt_path_read_set)
	#
	#
	# def check_path(self, reference_path, alternative_path):
	# 	condition = 0
	# 	read_set_pathalt_g_random = []
	# 	read_set_pathref_g_random = []
	# 	if alternative_path[0] not in self.dbg:
	# 		condition = 1
	# 		intersect_allnodes_alternative_path_g_random = set()
	# 		intersect_allnodes_reference_path_g_random = set()
	# 	else:
	# 		for i_node in range(0, len(alternative_path) - 1):
	# 			if alternative_path[i_node + 1] not in self.dbg[alternative_path[i_node]]:
	# 				condition = 1
	# 				intersect_allnodes_alternative_path_g_random = set()
	# 				break
	# 			read_set_pathalt_g_random.append(set(self.dbg.node[alternative_path[i_node]]['read_list_n']))
	# 		if condition == 0:
	# 			read_set_pathalt_g_random.append(set(self.dbg.node[alternative_path[i_node + 1]]['read_list_n']))
	# 			intersect_allnodes_alternative_path_g_random = set.intersection(*read_set_pathalt_g_random)
	# 		condition = 0
	# 		for i_node in range(0, len(reference_path) - 1):
	# 			if reference_path[i_node + 1] not in self.dbg[reference_path[i_node]]:
	# 				condition = 1
	# 				intersect_allnodes_reference_path_g_random = set()
	# 				break
	# 			read_set_pathref_g_random.append(set(self.dbg.node[reference_path[i_node]]['read_list_n']))
	# 		if condition == 0:
	# 			read_set_pathref_g_random.append(set(self.dbg.node[reference_path[i_node + 1]]['read_list_n']))
	# 			intersect_allnodes_reference_path_g_random = set.intersection(*read_set_pathref_g_random)
	# 	ratio = float(len(intersect_allnodes_alternative_path_g_random)) / (len(intersect_allnodes_alternative_path_g_random) + len(intersect_allnodes_reference_path_g_random))
	# 	return ratio, len(intersect_allnodes_reference_path_g_random), len(intersect_allnodes_alternative_path_g_random)

#
# #
# logger.info("will do")
an_alt_path = ['AAGCTCCCAGAATGCCAGAG', 'AGCTCCCAGAATGCCAGAGC', 'GCTCCCAGAATGCCAGAGCG', 'CTCCCAGAATGCCAGAGCGC', 'TCCCAGAATGCCAGAGCGCT', 'CCCAGAATGCCAGAGCGCTG',
			   'CCAGAATGCCAGAGCGCTGC', 'CAGAATGCCAGAGCGCTGCT', 'AGAATGCCAGAGCGCTGCTC', 'GAATGCCAGAGCGCTGCTCA', 'AATGCCAGAGCGCTGCTCAG', 'ATGCCAGAGCGCTGCTCAGA',
			   'TGCCAGAGCGCTGCTCAGAT', 'GCCAGAGCGCTGCTCAGATA', 'CCAGAGCGCTGCTCAGATAG', 'CAGAGCGCTGCTCAGATAGC', 'AGAGCGCTGCTCAGATAGCG', 'GAGCGCTGCTCAGATAGCGA']

# initial_seq="".join(an_alt_path)



a_ref_path = ['AAGCTCCCAGAATGCCAGAG', 'AGCTCCCAGAATGCCAGAGG', 'GCTCCCAGAATGCCAGAGGC', 'CTCCCAGAATGCCAGAGGCT', 'TCCCAGAATGCCAGAGGCTG', 'CCCAGAATGCCAGAGGCTGC',
			  'CCAGAATGCCAGAGGCTGCT', 'CAGAATGCCAGAGGCTGCTC', 'AGAATGCCAGAGGCTGCTCC', 'GAATGCCAGAGGCTGCTCCC', 'AATGCCAGAGGCTGCTCCCC', 'ATGCCAGAGGCTGCTCCCCG',
			  'TGCCAGAGGCTGCTCCCCGC', 'GCCAGAGGCTGCTCCCCGCG', 'CCAGAGGCTGCTCCCCGCGT', 'CAGAGGCTGCTCCCCGCGTG', 'AGAGGCTGCTCCCCGCGTGG', 'GAGGCTGCTCCCCGCGTGGC',
			  'AGGCTGCTCCCCGCGTGGCC', 'GGCTGCTCCCCGCGTGGCCC', 'GCTGCTCCCCGCGTGGCCCC', 'CTGCTCCCCGCGTGGCCCCT', 'TGCTCCCCGCGTGGCCCCTG', 'GCTCCCCGCGTGGCCCCTGC',
			  'CTCCCCGCGTGGCCCCTGCA', 'TCCCCGCGTGGCCCCTGCAC', 'CCCCGCGTGGCCCCTGCACC', 'CCCGCGTGGCCCCTGCACCA', 'CCGCGTGGCCCCTGCACCAG', 'CGCGTGGCCCCTGCACCAGC',
			  'GCGTGGCCCCTGCACCAGCA', 'CGTGGCCCCTGCACCAGCAG', 'GTGGCCCCTGCACCAGCAGC', 'TGGCCCCTGCACCAGCAGCT', 'GGCCCCTGCACCAGCAGCTC', 'GCCCCTGCACCAGCAGCTCC',
			  'CCCCTGCACCAGCAGCTCCT', 'CCCTGCACCAGCAGCTCCTA', 'CCTGCACCAGCAGCTCCTAC', 'CTGCACCAGCAGCTCCTACA', 'TGCACCAGCAGCTCCTACAC', 'GCACCAGCAGCTCCTACACC',
			  'CACCAGCAGCTCCTACACCG', 'ACCAGCAGCTCCTACACCGG', 'CCAGCAGCTCCTACACCGGC', 'CAGCAGCTCCTACACCGGCG', 'AGCAGCTCCTACACCGGCGG', 'GCAGCTCCTACACCGGCGGC',
			  'CAGCTCCTACACCGGCGGCC', 'AGCTCCTACACCGGCGGCCC', 'GCTCCTACACCGGCGGCCCC', 'CTCCTACACCGGCGGCCCCT', 'TCCTACACCGGCGGCCCCTG', 'CCTACACCGGCGGCCCCTGC',
			  'CTACACCGGCGGCCCCTGCA', 'TACACCGGCGGCCCCTGCAC', 'ACACCGGCGGCCCCTGCACC', 'CACCGGCGGCCCCTGCACCA', 'ACCGGCGGCCCCTGCACCAG', 'CCGGCGGCCCCTGCACCAGC',
			  'CGGCGGCCCCTGCACCAGCC', 'GGCGGCCCCTGCACCAGCCC', 'GCGGCCCCTGCACCAGCCCC', 'CGGCCCCTGCACCAGCCCCC', 'GGCCCCTGCACCAGCCCCCT', 'GCCCCTGCACCAGCCCCCTC',
			  'CCCCTGCACCAGCCCCCTCC', 'CCCTGCACCAGCCCCCTCCT', 'CCTGCACCAGCCCCCTCCTG', 'CTGCACCAGCCCCCTCCTGG', 'TGCACCAGCCCCCTCCTGGC', 'GCACCAGCCCCCTCCTGGCC',
			  'CACCAGCCCCCTCCTGGCCC', 'ACCAGCCCCCTCCTGGCCCC', 'CCAGCCCCCTCCTGGCCCCT', 'CAGCCCCCTCCTGGCCCCTG', 'AGCCCCCTCCTGGCCCCTGT', 'GCCCCCTCCTGGCCCCTGTC',
			  'CCCCCTCCTGGCCCCTGTCA', 'CCCCTCCTGGCCCCTGTCAT', 'CCCTCCTGGCCCCTGTCATC', 'CCTCCTGGCCCCTGTCATCT', 'CTCCTGGCCCCTGTCATCTT', 'TCCTGGCCCCTGTCATCTTC',
			  'CCTGGCCCCTGTCATCTTCT', 'CTGGCCCCTGTCATCTTCTG', 'TGGCCCCTGTCATCTTCTGT', 'GGCCCCTGTCATCTTCTGTC', 'GCCCCTGTCATCTTCTGTCC', 'CCCCTGTCATCTTCTGTCCC',
			  'CCCTGTCATCTTCTGTCCCT', 'CCTGTCATCTTCTGTCCCTT', 'CTGTCATCTTCTGTCCCTTC', 'TGTCATCTTCTGTCCCTTCC', 'GTCATCTTCTGTCCCTTCCC', 'TCATCTTCTGTCCCTTCCCA',
			  'CATCTTCTGTCCCTTCCCAG', 'ATCTTCTGTCCCTTCCCAGA', 'TCTTCTGTCCCTTCCCAGAA', 'CTTCTGTCCCTTCCCAGAAA', 'TTCTGTCCCTTCCCAGAAAA', 'TCTGTCCCTTCCCAGAAAAC',
			  'CTGTCCCTTCCCAGAAAACC', 'TGTCCCTTCCCAGAAAACCT', 'GTCCCTTCCCAGAAAACCTA', 'TCCCTTCCCAGAAAACCTAC', 'CCCTTCCCAGAAAACCTACC', 'CCTTCCCAGAAAACCTACCA',
			  'CTTCCCAGAAAACCTACCAG', 'TTCCCAGAAAACCTACCAGG', 'TCCCAGAAAACCTACCAGGG', 'CCCAGAAAACCTACCAGGGC', 'CCAGAAAACCTACCAGGGCA', 'CAGAAAACCTACCAGGGCAG',
			  'AGAAAACCTACCAGGGCAGC', 'GAAAACCTACCAGGGCAGCT', 'AAAACCTACCAGGGCAGCTA', 'AAACCTACCAGGGCAGCTAC', 'AACCTACCAGGGCAGCTACG', 'ACCTACCAGGGCAGCTACGG',
			  'CCTACCAGGGCAGCTACGGT', 'CTACCAGGGCAGCTACGGTT', 'TACCAGGGCAGCTACGGTTT', 'ACCAGGGCAGCTACGGTTTC', 'CCAGGGCAGCTACGGTTTCC', 'CAGGGCAGCTACGGTTTCCG',
			  'AGGGCAGCTACGGTTTCCGT', 'GGGCAGCTACGGTTTCCGTC', 'GGCAGCTACGGTTTCCGTCT', 'GCAGCTACGGTTTCCGTCTG', 'CAGCTACGGTTTCCGTCTGG', 'AGCTACGGTTTCCGTCTGGG',
			  'GCTACGGTTTCCGTCTGGGC', 'CTACGGTTTCCGTCTGGGCT', 'TACGGTTTCCGTCTGGGCTT', 'ACGGTTTCCGTCTGGGCTTC', 'CGGTTTCCGTCTGGGCTTCT', 'GGTTTCCGTCTGGGCTTCTT',
			  'GTTTCCGTCTGGGCTTCTTG', 'TTTCCGTCTGGGCTTCTTGC', 'TTCCGTCTGGGCTTCTTGCA', 'TCCGTCTGGGCTTCTTGCAT', 'CCGTCTGGGCTTCTTGCATT', 'CGTCTGGGCTTCTTGCATTC',
			  'GTCTGGGCTTCTTGCATTCT', 'TCTGGGCTTCTTGCATTCTG', 'CTGGGCTTCTTGCATTCTGG', 'TGGGCTTCTTGCATTCTGGG', 'GGGCTTCTTGCATTCTGGGA', 'GGCTTCTTGCATTCTGGGAC',
			  'GCTTCTTGCATTCTGGGACA', 'CTTCTTGCATTCTGGGACAG', 'TTCTTGCATTCTGGGACAGC', 'TCTTGCATTCTGGGACAGCC', 'CTTGCATTCTGGGACAGCCA', 'TTGCATTCTGGGACAGCCAA',
			  'TGCATTCTGGGACAGCCAAG', 'GCATTCTGGGACAGCCAAGT', 'CATTCTGGGACAGCCAAGTC', 'ATTCTGGGACAGCCAAGTCT', 'TTCTGGGACAGCCAAGTCTG', 'TCTGGGACAGCCAAGTCTGT',
			  'CTGGGACAGCCAAGTCTGTG', 'TGGGACAGCCAAGTCTGTGA', 'GGGACAGCCAAGTCTGTGAC', 'GGACAGCCAAGTCTGTGACT', 'GACAGCCAAGTCTGTGACTT', 'ACAGCCAAGTCTGTGACTTG',
			  'CAGCCAAGTCTGTGACTTGC', 'AGCCAAGTCTGTGACTTGCA', 'GCCAAGTCTGTGACTTGCAC', 'CCAAGTCTGTGACTTGCACG', 'CAAGTCTGTGACTTGCACGT', 'AAGTCTGTGACTTGCACGTA',
			  'AGTCTGTGACTTGCACGTAC', 'GTCTGTGACTTGCACGTACT', 'TCTGTGACTTGCACGTACTC', 'CTGTGACTTGCACGTACTCC', 'TGTGACTTGCACGTACTCCC', 'GTGACTTGCACGTACTCCCC',
			  'TGACTTGCACGTACTCCCCT', 'GACTTGCACGTACTCCCCTG', 'ACTTGCACGTACTCCCCTGT', 'CTTGCACGTACTCCCCTGTC', 'TTGCACGTACTCCCCTGTCC', 'TGCACGTACTCCCCTGTCCT',
			  'GCACGTACTCCCCTGTCCTC', 'CACGTACTCCCCTGTCCTCA', 'ACGTACTCCCCTGTCCTCAA', 'CGTACTCCCCTGTCCTCAAC', 'GTACTCCCCTGTCCTCAACA', 'TACTCCCCTGTCCTCAACAA',
			  'ACTCCCCTGTCCTCAACAAG', 'CTCCCCTGTCCTCAACAAGA', 'TCCCCTGTCCTCAACAAGAT', 'CCCCTGTCCTCAACAAGATG', 'CCCTGTCCTCAACAAGATGT', 'CCTGTCCTCAACAAGATGTT',
			  'CTGTCCTCAACAAGATGTTT', 'TGTCCTCAACAAGATGTTTT', 'GTCCTCAACAAGATGTTTTG', 'TCCTCAACAAGATGTTTTGC', 'CCTCAACAAGATGTTTTGCC', 'CTCAACAAGATGTTTTGCCA',
			  'TCAACAAGATGTTTTGCCAA', 'CAACAAGATGTTTTGCCAAC', 'AACAAGATGTTTTGCCAACT', 'ACAAGATGTTTTGCCAACTG', 'CAAGATGTTTTGCCAACTGG', 'AAGATGTTTTGCCAACTGGC',
			  'AGATGTTTTGCCAACTGGCC', 'GATGTTTTGCCAACTGGCCA', 'ATGTTTTGCCAACTGGCCAA', 'TGTTTTGCCAACTGGCCAAG', 'GTTTTGCCAACTGGCCAAGA', 'TTTTGCCAACTGGCCAAGAC',
			  'TTTGCCAACTGGCCAAGACC', 'TTGCCAACTGGCCAAGACCT', 'TGCCAACTGGCCAAGACCTG', 'GCCAACTGGCCAAGACCTGC', 'CCAACTGGCCAAGACCTGCC', 'CAACTGGCCAAGACCTGCCC',
			  'AACTGGCCAAGACCTGCCCT', 'ACTGGCCAAGACCTGCCCTG', 'CTGGCCAAGACCTGCCCTGT', 'TGGCCAAGACCTGCCCTGTG', 'GGCCAAGACCTGCCCTGTGC', 'GCCAAGACCTGCCCTGTGCA',
			  'CCAAGACCTGCCCTGTGCAG', 'CAAGACCTGCCCTGTGCAGC', 'AAGACCTGCCCTGTGCAGCT', 'AGACCTGCCCTGTGCAGCTG', 'GACCTGCCCTGTGCAGCTGT', 'ACCTGCCCTGTGCAGCTGTG',
			  'CCTGCCCTGTGCAGCTGTGG', 'CTGCCCTGTGCAGCTGTGGG', 'TGCCCTGTGCAGCTGTGGGT', 'GCCCTGTGCAGCTGTGGGTT', 'CCCTGTGCAGCTGTGGGTTG', 'CCTGTGCAGCTGTGGGTTGA',
			  'CTGTGCAGCTGTGGGTTGAT', 'TGTGCAGCTGTGGGTTGATT', 'GTGCAGCTGTGGGTTGATTC', 'TGCAGCTGTGGGTTGATTCC', 'GCAGCTGTGGGTTGATTCCA', 'CAGCTGTGGGTTGATTCCAC',
			  'AGCTGTGGGTTGATTCCACA', 'GCTGTGGGTTGATTCCACAC', 'CTGTGGGTTGATTCCACACC', 'TGTGGGTTGATTCCACACCC', 'GTGGGTTGATTCCACACCCC', 'TGGGTTGATTCCACACCCCC',
			  'GGGTTGATTCCACACCCCCG', 'GGTTGATTCCACACCCCCGC', 'GTTGATTCCACACCCCCGCC', 'TTGATTCCACACCCCCGCCC', 'TGATTCCACACCCCCGCCCG', 'GATTCCACACCCCCGCCCGG',
			  'ATTCCACACCCCCGCCCGGT', 'TTCCACACCCCCGCCCGGTA', 'TCCACACCCCCGCCCGGTAC', 'CCACACCCCCGCCCGGTACC', 'CACACCCCCGCCCGGTACCC', 'ACACCCCCGCCCGGTACCCG',
			  'CACCCCCGCCCGGTACCCGC', 'ACCCCCGCCCGGTACCCGCG', 'CCCCCGCCCGGTACCCGCGT', 'CCCCGCCCGGTACCCGCGTC', 'CCCGCCCGGTACCCGCGTCC', 'CCGCCCGGTACCCGCGTCCG',
			  'CGCCCGGTACCCGCGTCCGC', 'GCCCGGTACCCGCGTCCGCG', 'CCCGGTACCCGCGTCCGCGC', 'CCGGTACCCGCGTCCGCGCC', 'CGGTACCCGCGTCCGCGCCA', 'GGTACCCGCGTCCGCGCCAT',
			  'GTACCCGCGTCCGCGCCATG', 'TACCCGCGTCCGCGCCATGG', 'ACCCGCGTCCGCGCCATGGC', 'CCCGCGTCCGCGCCATGGCC', 'CCGCGTCCGCGCCATGGCCA', 'CGCGTCCGCGCCATGGCCAT',
			  'GCGTCCGCGCCATGGCCATC', 'CGTCCGCGCCATGGCCATCT', 'GTCCGCGCCATGGCCATCTA', 'TCCGCGCCATGGCCATCTAC', 'CCGCGCCATGGCCATCTACA', 'CGCGCCATGGCCATCTACAA',
			  'GCGCCATGGCCATCTACAAG', 'CGCCATGGCCATCTACAAGC', 'GCCATGGCCATCTACAAGCA', 'CCATGGCCATCTACAAGCAG', 'CATGGCCATCTACAAGCAGT', 'ATGGCCATCTACAAGCAGTC',
			  'TGGCCATCTACAAGCAGTCA', 'GGCCATCTACAAGCAGTCAC', 'GCCATCTACAAGCAGTCACA', 'CCATCTACAAGCAGTCACAG', 'CATCTACAAGCAGTCACAGC', 'ATCTACAAGCAGTCACAGCA',
			  'TCTACAAGCAGTCACAGCAC', 'CTACAAGCAGTCACAGCACA', 'TACAAGCAGTCACAGCACAT', 'ACAAGCAGTCACAGCACATG', 'CAAGCAGTCACAGCACATGA', 'AAGCAGTCACAGCACATGAC',
			  'AGCAGTCACAGCACATGACG', 'GCAGTCACAGCACATGACGG', 'CAGTCACAGCACATGACGGA', 'AGTCACAGCACATGACGGAG', 'GTCACAGCACATGACGGAGG', 'TCACAGCACATGACGGAGGT',
			  'CACAGCACATGACGGAGGTT', 'ACAGCACATGACGGAGGTTG', 'CAGCACATGACGGAGGTTGT', 'AGCACATGACGGAGGTTGTG', 'GCACATGACGGAGGTTGTGA', 'CACATGACGGAGGTTGTGAG',
			  'ACATGACGGAGGTTGTGAGG', 'CATGACGGAGGTTGTGAGGC', 'ATGACGGAGGTTGTGAGGCG', 'TGACGGAGGTTGTGAGGCGC', 'GACGGAGGTTGTGAGGCGCT', 'ACGGAGGTTGTGAGGCGCTG',
			  'CGGAGGTTGTGAGGCGCTGC', 'GGAGGTTGTGAGGCGCTGCC', 'GAGGTTGTGAGGCGCTGCCC', 'AGGTTGTGAGGCGCTGCCCC', 'GGTTGTGAGGCGCTGCCCCC', 'GTTGTGAGGCGCTGCCCCCA',
			  'TTGTGAGGCGCTGCCCCCAC', 'TGTGAGGCGCTGCCCCCACC', 'GTGAGGCGCTGCCCCCACCA', 'TGAGGCGCTGCCCCCACCAT', 'GAGGCGCTGCCCCCACCATG', 'AGGCGCTGCCCCCACCATGA',
			  'GGCGCTGCCCCCACCATGAG', 'GCGCTGCCCCCACCATGAGC', 'CGCTGCCCCCACCATGAGCG', 'GCTGCCCCCACCATGAGCGC', 'CTGCCCCCACCATGAGCGCT', 'TGCCCCCACCATGAGCGCTG',
			  'GCCCCCACCATGAGCGCTGC', 'CCCCCACCATGAGCGCTGCT', 'CCCCACCATGAGCGCTGCTC', 'CCCACCATGAGCGCTGCTCA', 'CCACCATGAGCGCTGCTCAG', 'CACCATGAGCGCTGCTCAGA',
			  'ACCATGAGCGCTGCTCAGAT', 'CCATGAGCGCTGCTCAGATA', 'CATGAGCGCTGCTCAGATAG', 'ATGAGCGCTGCTCAGATAGC', 'TGAGCGCTGCTCAGATAGCG', 'GAGCGCTGCTCAGATAGCGA']


# %timeit SL.sampling({"N": 120, "C": 120})
# %timeit a_rrg = RandomReadsGraph({"N": 1200, "C": 1200}, k=20,restrict_to=an_alt_path+a_ref_path)
# a_rrg.build_read_set_for_path(an_alt_path, verbose=True)
# len(a_ref_path)
# a_rrg.build_read_set_for_path(a_ref_path, verbose=True)
#
# a_rrg.check_path(a_ref_path, an_alt_path)
# a_rrg.check_path(a_ref_path, an_alt_path)
# logger.info("done")
