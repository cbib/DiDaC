#!/usr/bin/env python
# coding=utf8
## imports

import sys
import os
import time 
import networkx as nx 

## parameters

analyse_directory='/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/20150413/dbg/scripts' 
destination_directory='/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/20150413/dbg/outputs' 
os.chdir(analyse_directory)
k = 20

## imports 
import reference_graph as RG
import visualization as VISU
from individu_graph import individu_graph as IG
from randomreads_graph import randomreads_graph as RRG

# G_ref construction
G_ref = RG.ref_constructor(k)
# G_ref_merge = VISU.merge_reference_graph(G_ref.copy())
# G_ref_visu = VISU.reference_graph_visualization_formatting(G_ref.copy())
# G_ref_merge_visu = VISU.reference_graph_merged_visualization_formatting(G_ref_merge.copy())
# nx.write_gml(G_ref_visu,destination_directory+"/G_ref_visu"+str(k)+".gml")
# nx.write_gml(G_ref_merge_visu,destination_directory+"/G_ref_merge_visu"+str(k)+".gml")

# G_ind construction
# fastq = ["/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/data/fastq/fastq_pool0_trimmed/C_158_1.fastq","/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/data/fastq/fastq_pool0_trimmed/N_158_1.fastq"]
fastq = ["/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/20150413/dbg/data/fatsq_test/C_158_1.fastq","/Users/rudewicz/Documents/TP53_Pacbio/Analyses/GDB/20150413/dbg/data/fatsq_test/N_158_1.fastq"]
G_test = IG(fastq,k)
G_test.graph_cleaned_init(3) # .dbgclean creation
G_test.graph_rmRefEdges_init(G_test.dbgclean,G_ref) # .dbg_refrm creation
# G_test.dict_path_alternative_init(G_ref) # .path_alternative creation
# For visualisation
# Graph 
G_test_visu = VISU.individu_graph_visualization_formating(G_test.dbg.copy(),G_ref.copy())
G_test_clean_visu = VISU.individu_graph_visualization_formating(G_test.dbgclean.copy(),G_ref.copy())
nx.write_gml(G_test_visu,destination_directory+"/G_158_"+str(k)+".gml")
nx.write_gml(G_test_clean_visu,destination_directory+"/G_158_clean3_"+str(k)+".gml")
# Graph merged
G_test_merged = VISU.merge_individu_graph(G_test.dbg.copy(),G_ref.copy())
G_test_clean_merged = VISU.merge_individu_graph(G_test.dbgclean.copy(),G_ref.copy())
G_test_merged_visu = VISU.individu_graph_merged_visualization_formating(G_test_merged.copy(),G_ref.copy())
G_test_clean_merged_visu = VISU.individu_graph_merged_visualization_formating(G_test_clean_merged.copy(),G_ref.copy())
nx.write_gml(G_test_merged_visu,destination_directory+"/G_158_merged_"+str(k)+".gml")
nx.write_gml(G_test_clean_merged_visu,destination_directory+"/G_158_clean3_merged_"+str(k)+".gml")

### HERE ###

G_test.alteration_list_init(G_ref,k) # .alteration_list creation

print time.strftime('%d/%m/%y %H:%M',time.localtime())
for i in range(0,1000):
	G_random = RRG(G_test.coverage,k)
	for i_alteration in range(0,len(G_test.alteration_list)):
		G_test.alteration_list[i_alteration].random_count_list.append(G_random.check_path(G_test.alteration_list[i_alteration].alternative_path))	

print time.strftime('%d/%m/%y %H:%M',time.localtime())


for i_alteration in range(0,len(G_test.alteration_list)):
	# G_test.alteration_list[i_alteration].pvalue_init()
	print"%s vs %s pvalue: %f (%d)"%(G_test.alteration_list[i_alteration].reference_sequence,G_test.alteration_list[i_alteration].alternative_sequence,G_test.alteration_list[i_alteration].pvalue,G_test.alteration_list[i_alteration].read_count)





#### Test rm nodes in alternative_path where p.value < 0.001
G_test_pval_clean = G_test.dbgclean.copy()	
for i_alteration in range(0,len(G_test.alteration_list)):
	if G_test.alteration_list[i_alteration].pvalue > 0.001:
		print (G_test.alteration_list[i_alteration].pvalue)
		nodes_to_rm = G_test.alteration_list[i_alteration].alternative_path[1:len(G_test.alteration_list[i_alteration].alternative_path)-2]
		G_test_pval_clean.remove_nodes_from(nodes_to_rm)

# Visualization
G_test_pval_clean_merged = VISU.merge_individu_graph(G_test_pval_clean.copy(),G_ref.copy())
G_test_pval_clean_merged_visu = VISU.individu_graph_merged_visualization_formating(G_test_pval_clean_merged.copy(),G_ref.copy())
nx.write_gml(G_test_pval_clean_merged_visu,destination_directory+"/G_test_pval_clean_merged_visu_3_158_3"+str(k)+".gml")
