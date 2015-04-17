rule some_samples:
	input: "data/all_pool_trimmed0.1/C_2_1.fastq"
	output: "output/alterations/sample_2_1.tsv"
	log: "output/logs/sample_2_1.txt"
	shell:
		"""
			python src/principal.py --samplekey 2 --cfastq {input} > {output} 2> {log}
		"""
