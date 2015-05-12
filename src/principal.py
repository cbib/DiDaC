#!/usr/bin/env python
# coding=utf8
## imports

import networkx as nx
from argparse import ArgumentParser
from helpers.helpers import time_iterator, get_or_create_dir
from helpers.logger import init_logger

### TEST ALIGNMENT
from Bio import pairwise2

### TEST ALIGNMENT

logger = init_logger('PACBIOP53')

## parameters

## imports
logger.info("Will import")
import reference_graph as RG
import visualization as VISU
from individugraph import IndividuGraph as IG
from randomreadsgraph import RandomReadsGraph as RRG

logger.info("Import finished")



def process_sample(kmer_length, sample_key=None, c_fastq_file=None, n_fastq_file=None, min_support_percentage=3, n_permutations=1000, destination_directory=".", export_gml=False):
	# global g_ref, fastq, g_test, g_test_visu, g_test_clean_visu, g_test_merged, g_test_clean_merged, g_test_merged_visu, g_test_clean_merged_visu, i, g_random, i_alteration, G_test_pval_clean, nodes_to_rm, G_test_pval_clean_merged, G_test_pval_clean_merged_visu

	# g_ref construction
	logger.info("Will build reference graph with k==%d", kmer_length)
	g_ref = RG.ref_constructor(kmer_length)


	# TODO: Test presence of cycles, abort if yes

	# G_ref_merge = VISU.merge_reference_graph(g_ref.copy())
	# G_ref_visu = VISU.reference_graph_visualization_formatting(g_ref.copy())
	# G_ref_merge_visu = VISU.reference_graph_merged_visualization_formatting(G_ref_merge.copy())
	# nx.write_gml(G_ref_visu,destination_directory+"/G_ref_visu"+str(k)+".gml")
	# nx.write_gml(G_ref_merge_visu,destination_directory+"/G_ref_merge_visu"+str(k)+".gml")
	# G_ind construction

	fastq = [c_fastq_file, n_fastq_file]
	fastq = [f for f in fastq if f]

	logger.info("Will build sample graph for %s with k==%d and minimum support (percentage) = %d", fastq, kmer_length, min_support_percentage)
	g_test = IG(fastq, kmer_length)
	g_test.graph_cleaned_init(min_support_percentage)  # .dbgclean creation
	g_test.graph_rmRefEdges_init(g_test.dbgclean, g_ref)  # .dbg_refrm creation

	# For visualisation
	graph_name = "G_%s_" % sample_key
	if export_gml:
		logger.info("Will save viz graph for %s with k==%d", fastq, kmer_length)
		get_or_create_dir(destination_directory)
		g_test_visu = VISU.individu_graph_visualization_formating(g_test.dbg.copy(), g_ref.copy())
		g_test_clean_visu = VISU.individu_graph_visualization_formating(g_test.dbgclean.copy(), g_ref.copy())
		cleaned_graph_name = graph_name + "clean%d_" % min_support_percentage
		nx.write_gml(g_test_visu, destination_directory + "/" + graph_name + str(kmer_length) + ".gml")
		nx.write_gml(g_test_clean_visu, destination_directory + "/" + cleaned_graph_name + str(kmer_length) + ".gml")

		# Graph merged
		logger.info("Will merge graph for %s with k==%d", fastq, kmer_length)
		g_test_merged = VISU.merge_individu_graph(g_test.dbg.copy(), g_ref.copy())
		g_test_merged_visu = VISU.individu_graph_merged_visualization_formating(g_test_merged.copy(), g_ref.copy())
		merged_graph_name = "G_%s_merged_" % sample_key
		nx.write_gml(g_test_merged_visu, destination_directory + "/" + merged_graph_name + str(kmer_length) + ".gml")

		g_test_clean_merged = VISU.merge_individu_graph(g_test.dbgclean.copy(), g_ref.copy())
		g_test_clean_merged_visu = VISU.individu_graph_merged_visualization_formating(g_test_clean_merged.copy(), g_ref.copy())

		merged_cleaned_graph_name = graph_name + "clean%d_merged_" % min_support_percentage

		nx.write_gml(g_test_clean_merged_visu, destination_directory + "/" + merged_cleaned_graph_name + str(kmer_length) + ".gml")

	### Permutation test ###

	g_test.alteration_list_init(g_ref, kmer_length)  # .alteration_list creation

	logger.info("Will create random graphs")
	all_possible_kmers=set()
	for an_alt in g_test.alteration_list:
		all_possible_kmers.update(an_alt.reference_path)
		all_possible_kmers.update(an_alt.alternative_path)

	for i, j in time_iterator(range(0, n_permutations), logger, msg_prefix="permuting"):
		g_random = RRG(g_test.coverage, kmer_length,restrict_to=all_possible_kmers)
		for i_alteration in range(0, len(g_test.alteration_list)):
			g_random_data = g_random.check_path(g_test.alteration_list[i_alteration].reference_path, g_test.alteration_list[i_alteration].alternative_path)
			g_test.alteration_list[i_alteration].random_ratio_list.append(g_random_data[0])
			g_test.alteration_list[i_alteration].random_reference_count_list.append(g_random_data[1])
			g_test.alteration_list[i_alteration].random_alternative_count_list.append(g_random_data[2])

	logger.info("Will generate p-values")
	for i_alteration in range(0, len(g_test.alteration_list)):
		g_test.alteration_list[i_alteration].pvalue_init()
		print "%s\t%s\t%s\t%f\t%f" % (
			sample_key,
			g_test.alteration_list[i_alteration].reference_sequence,
			g_test.alteration_list[i_alteration].alternative_sequence,
			g_test.alteration_list[i_alteration].pvalue_ratio,
			g_test.alteration_list[i_alteration].ratio_read_count
		)

	g_test.significant_alteration_list_init()
	# ICI

	print "Are significant:"
	for i_alteration in range(0, len(g_test.significant_alteration_list)):
		print "%d\t%s\t%s\t%s\t%f\t%f" % (
		i_alteration,
		sample_key,
		g_test.significant_alteration_list[i_alteration].reference_sequence,
		g_test.significant_alteration_list[i_alteration].alternative_sequence,
		g_test.significant_alteration_list[i_alteration].pvalue_ratio,
		g_test.significant_alteration_list[i_alteration].ratio_read_count
		)
	
	if len(g_test.significant_alteration_list) > 1:
		g_test.multiple_alternative_path_filter()


	print "After filtering:"
	for i_alteration in range(0, len(g_test.significant_alteration_list)):
		alignments = pairwise2.align.globalms(g_test.significant_alteration_list[i_alteration].reference_sequence, g_test.significant_alteration_list[i_alteration].alternative_sequence, 2, -3, -5, -2)
		print "%d\t%s\t%s\t%s\t%f\t%f" % (
		i_alteration,
		sample_key,
		g_test.significant_alteration_list[i_alteration].reference_sequence,
		g_test.significant_alteration_list[i_alteration].alternative_sequence,
		g_test.significant_alteration_list[i_alteration].pvalue_ratio,
		g_test.significant_alteration_list[i_alteration].ratio_read_count
		)
		print alignments
		if len(alignments) > 1:
			logger.critical("More than one alignements for %s vs %s", g_test.significant_alteration_list[i_alteration].reference_sequence,g_test.significant_alteration_list[i_alteration].alternative_sequence) 
			alignments = [alignments[0]]
		uncompact_cigar = ""
		compact_cigard = []
		for i_nucleotide in range(0,alignments[0][4]):
			if alignments[0][0][i_nucleotide] == alignments[0][1][i_nucleotide] :
				uncompact_cigar += "M"
			elif alignments[0][0][i_nucleotide] == "-":
				uncompact_cigar += "I"
			elif alignments[0][1][i_nucleotide] == "-":
				uncompact_cigar += "D"
			else:
				uncompact_cigar += "X"
		print uncompact_cigar
		operation = uncompact_cigar[0]
		count = 0
		for i_nucleotide in range(0,len(uncompact_cigar)):
			if uncompact_cigar[i_nucleotide] != operation:
				compact_cigard += [count,operation]
				operation = uncompact_cigar[i_nucleotide]
				count = 0
			count += 1
		compact_cigard += [count,operation]
		print compact_cigard
		if len(compact_cigard) == 6:
			# TO CONTINUE
			alteration_type = compact_cigard[3]
			if alteration_type == "X":
				# c.76A>T
				reference = g_test.significant_alteration_list[i_alteration].reference_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]
				alteration = g_test.significant_alteration_list[i_alteration].alternative_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]
				position = g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']['NM_000546.5']+compact_cigard[0]
				print "NM_000546.5:c.%d%s>%s"%(position,reference,alteration)
			elif alteration_type == "D":
				# c.76_78delACT
				reference = g_test.significant_alteration_list[i_alteration].reference_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]
				position = g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']['NM_000546.5']+compact_cigard[0]+1
				print "NM_000546.5:c.%d_%ddel%s"%(position,position+len(reference)-1,reference)
			else:
				# c.76_77insG
				alteration = g_test.significant_alteration_list[i_alteration].alternative_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]
				position = g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']['NM_000546.5']+compact_cigard[0]+1
				print "NM_000546.5:c.%d_%dins%s"%(position,position+1,alteration)






if __name__ == "__main__":
	parser = ArgumentParser()
	parser.add_argument('--kmer_length', help='Size of k-mer words', default=20, type=int, required=False)
	parser.add_argument('--cfastq', help='FASTQ file for C terminal reads', required=False, type=str)
	parser.add_argument('--nfastq', help='FASTQ file for N terminal reads', required=False, type=str)
	parser.add_argument('--samplekey', help='Unique sample key', default="", type=str, required=True)
	parser.add_argument("--npermutations", help="number of permutations / random samples to perform", default=10, type=int, required=False)
	parser.add_argument("--destdir", help="Output directory", default="output/gml", type=str, required=False)
	parser.add_argument("--export", help="Whether to export graphs to GML", action='store_true')

	# k = 20
	# c_fastq_file = "dbg/data/fatsq_test/C_158_1.fastq"
	# n_fastq_file = "dbg/data/fatsq_test/N_158_1.fastq"
	# sample_key = "158"
	# n_permutations = 1000
	args = parser.parse_args()

	process_sample(kmer_length=args.kmer_length, c_fastq_file=args.cfastq, n_fastq_file=args.nfastq, n_permutations=args.npermutations, sample_key=args.samplekey,
					destination_directory=args.destdir, export_gml=args.export)
