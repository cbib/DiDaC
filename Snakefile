

rule process_sample_1:
	input: c_term="data/fastq/all_pool_trimmed0.1/C_{SAMPLE}.fastq", n_term="data/fastq/all_pool_trimmed0.1/N_{SAMPLE}.fastq"
	output: "output/alterations/sample_{SAMPLE}.tsv"
	log: "output/logs/sample_{SAMPLE}.txt"
	shell:
		"""
			python src/principal.py --samplekey {SAMPLE} --cfastq {input} > {output} 2> {log}
		"""

sam2test = open("data/samples2analyse.txt",'r').readlines()

rule some_samples:
	input: expand("output/alterations/sample_{id}_1.tsv",id=sam2test)

