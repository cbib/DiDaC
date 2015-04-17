#!/usr/bin/env python
# coding=utf8
class alteration:

	def __init__(self,reference_path,alternative_path,read_count,k):
		self.reference_path = reference_path
		self.alternative_path = alternative_path
		self.read_count = read_count
		self.random_count_list = []

		self.reference_sequence = alteration.kmerpathToSeq(self.reference_path,k)
		self.alternative_sequence = alteration.kmerpathToSeq(self.alternative_path,k)
		
	def pvalue_init(self):
		self.pvalue = float(len([i for i in self.random_count_list if i >= self.read_count])) / len(self.random_count_list)

	@staticmethod
	def kmerpathToSeq (kmer_list,k):
		sequence = kmer_list[0]
		for i_kmer in range(1,len(kmer_list)):
			sequence += kmer_list[i_kmer][k-1]
		return(sequence)