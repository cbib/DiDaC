#!/usr/bin/env python
# coding=utf8
## imports

import networkx as nx
from argparse import ArgumentParser
from helpers.helpers import time_iterator, get_or_create_dir
from helpers.logger import init_logger
import sys

logger = init_logger('PACBIOP53')

## parameters

## imports
logger.info("Will import")
import reference_graph as RG
import visualization as VISU
from individugraph import IndividuGraph as IG
from randomreadsgraph import RandomReadsGraph as RRG
import forannotation as ANNO

logger.info("Import finished")



def process_sample(kmer_length, min_support_percentage,  n_permutations, sample_key=None, c_fastq_file=None, n_fastq_file=None, destination_directory=".", export_gml=False):

	# g_ref construction
	logger.info("Will build reference graph with k==%d", kmer_length)
	g_ref = RG.ref_constructor(kmer_length)

	# g_ind construction
	fastq = [c_fastq_file, n_fastq_file]
	fastq = [f for f in fastq if f]

	logger.info("Will build sample graph for %s with k==%d and minimum support (percentage) = %d", fastq, kmer_length, min_support_percentage)
	g_test = IG(fastq, kmer_length)
	g_test.graph_cleaned_init(min_support_percentage)  # .dbgclean creation


	# Is there cycles ?
	if list(nx.simple_cycles(g_test.dbgclean)):
		if kmer_length > 50:
			logger.info("There are always cycle(s) with k==50...exiting")
			sys.exit(0)
		# Check non depassement valeur limite de k 
		return process_sample(kmer_length=kmer_length+1,sample_key=sample_key,c_fastq_file=c_fastq_file,n_fastq_file=n_fastq_file, min_support_percentage=min_support_percentage, n_permutations=n_permutations, destination_directory=destination_directory, export_gml=export_gml)

	# Some prints for stats 
	dir_stat = get_or_create_dir("output/statistics") 
	# graph stat
	graph_stat_file = open(dir_stat+"/graph_stat_file"+sample_key+".tsv", 'w')
	graph_stat_file.write(
		"%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d"%(
		kmer_length,
		g_ref.size(),
		sample_key,
		g_test.coverage['C'],
		g_test.coverage['N'],
		g_test.dbg.size(),
		g_test.dbgclean.size(),
		g_test.dbg.in_degree().values().count(0),
		g_test.dbg.out_degree().values().count(0),
		g_test.dbgclean.in_degree().values().count(0),
		g_test.dbgclean.out_degree().values().count(0)
		))

	# kmer stat
	kmer_stat_file = open(dir_stat+"/kmer_stat_file"+sample_key+".tsv", 'w')
	for node_print in g_test.dbg.nodes():
		fragment_print = "".join(g_test.dbg.node[node_print]['fragment'])
		reads_print = len(g_test.dbg.node[node_print]['read_list_n'])
		kmer_stat_file.write(
			"%s\t%s\t%s\t%d\n"%(
			sample_key,
			node_print,
			fragment_print,
			reads_print,
			))

	g_test.graph_rmRefEdges_init(g_test.dbgclean, g_ref)  # .dbg_refrm creation

	# For visualisation
	graph_name = "G_%s_" % sample_key
	if export_gml:
		logger.info("Will save viz graph for %s with k==%d", fastq, kmer_length)
		get_or_create_dir(destination_directory)
		G_ref_merge = VISU.merge_reference_graph(g_ref.copy())
		G_ref_visu = VISU.reference_graph_visualization_formatting(g_ref.copy())
		G_ref_merge_visu = VISU.reference_graph_merged_visualization_formatting(G_ref_merge.copy())
		nx.write_gml(G_ref_visu,destination_directory+"/G_ref_visu"+str(kmer_length)+".gml")
		nx.write_gml(G_ref_merge_visu,destination_directory+"/G_ref_merge_visu"+str(kmer_length)+".gml")
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

	# .alteration_list creation
	g_test.alteration_list_init(g_ref, kmer_length,min_support_percentage)  

	### Permutation test ###
	logger.info("Will create random graphs")
	all_possible_kmers=set()
	for an_alt in g_test.alteration_list:
		all_possible_kmers.update(an_alt.reference_path)
		all_possible_kmers.update(an_alt.alternative_path)

	for i, j in time_iterator(range(0, n_permutations), logger, msg_prefix="permuting"):
		g_random = RRG(g_test.coverage, kmer_length,restrict_to=all_possible_kmers)
		for i_alteration in range(0, len(g_test.alteration_list)):
			g_random_data = g_random.check_path(g_test.alteration_list[i_alteration].reference_path, g_test.alteration_list[i_alteration].alternative_path, g_test.alteration_list[i_alteration].min_coverage)
			g_test.alteration_list[i_alteration].random_ratio_list.append(g_random_data[0])
			g_test.alteration_list[i_alteration].random_reference_count_list.append(g_random_data[1])
			g_test.alteration_list[i_alteration].random_alternative_count_list.append(g_random_data[2])

	logger.info("Will generate p-values")
	for i_alteration in range(0, len(g_test.alteration_list)):
		g_test.alteration_list[i_alteration].pvalue_init()

	g_test.significant_alteration_list_init()
	
	# If more than one significant alteration, check if they are not in "spike" (en épis)
	if len(g_test.significant_alteration_list) > 1:
	 	g_test.multiple_alternative_path_filter()

	## Stat 
	# graph stat
	alt_stat_file = open(dir_stat+"/alt_stat_file"+sample_key+".tsv", 'w')
	for i_alteration in range(0, len(g_test.significant_alteration_list)):
		if g_test.significant_alteration_list[i_alteration].pvalue_ratio <= 1:
			# print "%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%f\t%f" % (
			# alt_stat_file.write("%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%s" % (				
			alt_stat_file.write("%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%f\t%f\t%s" % (				
			i_alteration+1,
			sample_key,
			g_test.coverage['C'],
			g_test.coverage['N'],
			g_test.significant_alteration_list[i_alteration].reference_sequence,
			g_test.significant_alteration_list[i_alteration].alternative_sequence,
			g_test.significant_alteration_list[i_alteration].reference_read_count,
			g_test.significant_alteration_list[i_alteration].alternative_read_count,
			g_test.significant_alteration_list[i_alteration].ratio_read_count,
			g_test.significant_alteration_list[i_alteration].pvalue_ratio,
			# g_test.significant_alteration_list[i_alteration].zscore,
			"\t".join(map(str,g_test.significant_alteration_list[i_alteration].random_ratio_list))
			))

	### MICADo + ###
	ANNO.alteration_list_to_transcrit_mutation(g_test,g_ref)

if __name__ == "__main__":
	parser = ArgumentParser()
	parser.add_argument('--kmer_length', help='Size of k-mer words', default=20, type=int, required=False)
	parser.add_argument('--cfastq', help='FASTQ file for C terminal reads', required=False, type=str)
	parser.add_argument('--nfastq', help='FASTQ file for N terminal reads', required=False, type=str)
	parser.add_argument('--min_support_percentage', help='Minimum of read support percentage for node filter', default=3, type=int, required=False)
	parser.add_argument('--samplekey', help='Unique sample key', default="", type=str, required=True)
	parser.add_argument('--npermutations', help="number of permutations / random samples to perform", default=1000, type=int, required=False)
	parser.add_argument("--destdir", help="Output directory", default="output/gml", type=str, required=False)
	parser.add_argument("--export", help="Whether to export graphs to GML", action='store_true')

	args = parser.parse_args()

	process_sample(kmer_length=args.kmer_length, min_support_percentage=args.min_support_percentage, c_fastq_file=args.cfastq, n_fastq_file=args.nfastq, n_permutations=args.npermutations, sample_key=args.samplekey,
					destination_directory=args.destdir, export_gml=args.export)
