from helpers.helpers import time_iterator, get_or_create_dir
from helpers.logger import init_logger
from Bio import pairwise2
logger = init_logger('Annotation')

def alteration_list_to_transcrit_mutation(g_test,g_ref):
	for i_alteration in range(0, len(g_test.significant_alteration_list)):
		alignments = pairwise2.align.globalms(g_test.significant_alteration_list[i_alteration].reference_sequence, g_test.significant_alteration_list[i_alteration].alternative_sequence, 2, -3, -5, -2)
		# print alignments
		if len(alignments) > 1:
			logger.critical("More than one alignment for %s vs %s", g_test.significant_alteration_list[i_alteration].reference_sequence,g_test.significant_alteration_list[i_alteration].alternative_sequence) 
			alignments = [alignments[len(alignments)-1]]
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
		# print uncompact_cigar
		operation = uncompact_cigar[0]
		count = 0
		for i_nucleotide in range(0,len(uncompact_cigar)):
			if uncompact_cigar[i_nucleotide] != operation:
				compact_cigard += [count,operation]
				operation = uncompact_cigar[i_nucleotide]
				count = 0
			count += 1
		compact_cigard += [count,operation]
		# print compact_cigard
		if len(compact_cigard) == 6:
			# TO CONTINUE
			alteration_type = compact_cigard[3]
			# print g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']
			if len(g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']) == 1:
				splicing_variant = g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list'].keys()[0]
			elif "NM_000546.5" in g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']:
				splicing_variant = "NM_000546.5"
			elif "NM_001126114.2" in g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list']:
				splicing_variant = "NM_001126114.2"
			if alteration_type == "X":
				# n.76A>T
				reference = g_test.significant_alteration_list[i_alteration].reference_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]
				alteration = g_test.significant_alteration_list[i_alteration].alternative_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]
				position = g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list'][splicing_variant]+compact_cigard[0]
				print "%s:n.%d%s>%s"%(splicing_variant,position,reference,alteration)
			elif alteration_type == "D":
				# n.76_78delACT
				reference = g_test.significant_alteration_list[i_alteration].reference_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]
				position = g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list'][splicing_variant]+compact_cigard[0]
				print "%s:n.%d_%ddel%s"%(splicing_variant,position,position+len(reference)-1,reference)
			else:
				# n.76_77insG
				alteration = g_test.significant_alteration_list[i_alteration].alternative_sequence[compact_cigard[0]:compact_cigard[0]+compact_cigard[2]]
				position = g_ref.node[g_test.significant_alteration_list[i_alteration].reference_path[0]]['ref_list'][splicing_variant]+compact_cigard[0]
				print "%s:n.%d_%dins%s"%(splicing_variant,position,position+1,alteration)
