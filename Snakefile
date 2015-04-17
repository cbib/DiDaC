

rule process_sample:
	input: c_term="data/fastq/all_pool_trimmed0.1/C_{SAMPLE}.fastq", n_term="data/fastq/all_pool_trimmed0.1/N_{SAMPLE}.fastq"
	output: "output/alterations/{SAMPLE}.tsv"
	log: "output/logs/{SAMPLE}.txt"
	shell:
		"""
			python src/principal.py --samplekey {SAMPLE} --cfastq {input.c_term} --nfastq {input.n_term} --npermutations 1000 > {output} 2> {log}
		"""

rule some_samples:
	sam2test = open("data/samples2analyse.txt",'r').readlines()
	sam2test = map(str.strip, sam2test)
	input: expand("output/alterations/{SAMPLE}.tsv",SAMPLE=sam2test)





