#!/usr/bin/env python
# coding=utf8
## imports

import sys
import os
import time
import networkx as nx
from argparse import ArgumentParser
import logging


## parameters

## imports
import reference_graph as RG
import visualization as VISU
from individu_graph import individu_graph as IG
from randomreads_graph import randomreads_graph as RRG

import logging

if "logger" not in globals():
	logger = logging.getLogger('PACBIOP53')
	logger.setLevel(logging.DEBUG)
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s', "%Y-%m-%d %H:%M:%S")
	ch.setFormatter(formatter)
	logger.addHandler(ch)

logger = logging.getLogger('PACBIOP53')


def process_sample(kmer_length, sample_key=None, c_fastq_file=None, n_fastq_file=None, min_support_percentage=3, n_permutations=1000, destination_directory="."):
	# global G_ref, fastq, G_test, G_test_visu, G_test_clean_visu, G_test_merged, G_test_clean_merged, G_test_merged_visu, G_test_clean_merged_visu, i, G_random, i_alteration, G_test_pval_clean, nodes_to_rm, G_test_pval_clean_merged, G_test_pval_clean_merged_visu

	# G_ref construction
	G_ref = RG.ref_constructor(kmer_length)


	# TODO: Test presence of cycles, abort if yes
	# TODO: Create output gml dir if needed

	# G_ref_merge = VISU.merge_reference_graph(G_ref.copy())
	# G_ref_visu = VISU.reference_graph_visualization_formatting(G_ref.copy())
	# G_ref_merge_visu = VISU.reference_graph_merged_visualization_formatting(G_ref_merge.copy())
	# nx.write_gml(G_ref_visu,destination_directory+"/G_ref_visu"+str(k)+".gml")
	# nx.write_gml(G_ref_merge_visu,destination_directory+"/G_ref_merge_visu"+str(k)+".gml")
	# G_ind construction

	fastq = [c_fastq_file, n_fastq_file]
	fastq = [f for f in fastq if f]

	G_test = IG(fastq, kmer_length)
	G_test.graph_cleaned_init(min_support_percentage)  # .dbgclean creation
	G_test.graph_rmRefEdges_init(G_test.dbgclean, G_ref)  # .dbg_refrm creation
	# G_test.dict_path_alternative_init(G_ref) # .path_alternative creation
	# For visualisation
	# Graph
	# G_test_visu = VISU.individu_graph_visualization_formating(G_test.dbg.copy(), G_ref.copy())
	# G_test_clean_visu = VISU.individu_graph_visualization_formating(G_test.dbgclean.copy(), G_ref.copy())

	# graph_name = "G_%s_" % sample_key
	# cleaned_graph_name = graph_name + "clean%d_" % min_support_percentage

	# # nx.write_gml(G_test_visu, destination_directory + "/" + graph_name + str(kmer_length) + ".gml")
	# # nx.write_gml(G_test_clean_visu, destination_directory + "/" + cleaned_graph_name + str(kmer_length) + ".gml")


	# # Graph merged
	# G_test_merged = VISU.merge_individu_graph(G_test.dbg.copy(), G_ref.copy())
	# G_test_clean_merged = VISU.merge_individu_graph(G_test.dbgclean.copy(), G_ref.copy())
	# G_test_merged_visu = VISU.individu_graph_merged_visualization_formating(G_test_merged.copy(), G_ref.copy())
	# G_test_clean_merged_visu = VISU.individu_graph_merged_visualization_formating(G_test_clean_merged.copy(), G_ref.copy())

	# merged_graph_name = "G_%s_merged_" % sample_key
	# merged_cleaned_graph_name = graph_name + "clean%d_merged_" % min_support_percentage

	# nx.write_gml(G_test_merged_visu, destination_directory + "/" + merged_graph_name + str(kmer_length) + ".gml")
	# nx.write_gml(G_test_clean_merged_visu, destination_directory + "/" + merged_cleaned_graph_name + str(kmer_length) + ".gml")

	### Permutation test ###

	G_test.alteration_list_init(G_ref, kmer_length)  # .alteration_list creation
	# print time.strftime('%d/%m/%y %H:%M', time.localtime())

	logger.info("Will create random graphs")
	for i in range(0, n_permutations):
		G_random = RRG(G_test.coverage, kmer_length)
		for i_alteration in range(0, len(G_test.alteration_list)):
			G_random_data = G_random.check_path(G_test.alteration_list[i_alteration].reference_path,G_test.alteration_list[i_alteration].alternative_path) 
			G_test.alteration_list[i_alteration].random_ratio_list.append(G_random_data[0])
			G_test.alteration_list[i_alteration].random_reference_count_list.append(G_random_data[1])
			G_test.alteration_list[i_alteration].random_alternative_count_list.append(G_random_data[2])

	logger.info("Will generate p-values")
	for i_alteration in range(0, len(G_test.alteration_list)):
		G_test.alteration_list[i_alteration].pvalue_init()
		print "%s\t%s\t%s\t%f\t%f" % (
		sample_key,
		G_test.alteration_list[i_alteration].reference_sequence,
		G_test.alteration_list[i_alteration].alternative_sequence,
		G_test.alteration_list[i_alteration].pvalue_ratio,
		G_test.alteration_list[i_alteration].ratio_read_count
		)


if __name__ == "__main__":
	parser = ArgumentParser()
	parser.add_argument('--kmer_length', help='Size of k-mer words', default=20, type=int, required=False)
	parser.add_argument('--cfastq', help='FASTQ file for C terminal reads', required=False, type=str)
	parser.add_argument('--nfastq', help='FASTQ file for N terminal reads', required=False, type=str)
	parser.add_argument('--samplekey', help='Unique sample key', default="", type=str, required=True)
	parser.add_argument("--npermutations", help="number of permutations / random samples to perform", default=10, type=int, required=False)
	parser.add_argument("--destdir", help="Output directory", default="output/gml", type=str, required=False)

	# k = 20
	# c_fastq_file = "dbg/data/fatsq_test/C_158_1.fastq"
	# n_fastq_file = "dbg/data/fatsq_test/N_158_1.fastq"
	# sample_key = "158"
	# n_permutations = 1000
	args = parser.parse_args()

	process_sample(kmer_length=args.kmer_length, c_fastq_file=args.cfastq, n_fastq_file=args.nfastq, n_permutations=args.npermutations, sample_key=args.samplekey,
				   destination_directory=args.destdir)
