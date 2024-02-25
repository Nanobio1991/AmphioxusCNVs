	'''
	Mask repeats in the genome using RepeatMasker
	'''
rule mask_genome_with_repeatmasker:
	input:
		amphioxus_genome = "data/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		masked_genome = "results/mask_genome_with_repeatmasker/Branchiostoma_lanceolatum.BraLan3_genome_masked.fa"
	log:
		err = "logs/mask_genome_with_repeatmasker/RepeatMasker.err",
		out = "logs/mask_genome_with_repeatmasker/RepeatMasker.out"
	conda:
		'envs/Finding_SDs.yaml'
	params:
		output_dir = "results/mask_genome_with_repeatmasker/",
		time = '10:00:00',
		threads = 2,
		mem = 20000,
		name = "RepeatMasker",
	shell:
		"RepeatMasker -s -xsmall -e ncbi {input.amphioxus_genome} -dir $(dirname {output.masked_genome}) >> {log.out} 2> {log.err}"


	'''
	Indexing genome with samtools
	'''
rule index_genome:
	input:
		masked_genome = rules.mask_genome_with_repeatmasker.output.masked_genome
	output:
		indexed_genome = "results/mask_genome_with_repeatmasker/Branchiostoma_lanceolatum.BraLan3_genome_masked.fa.fai"
	log:
		err = "logs/index_genome/index.err",
		out = "logs/index_genome/index.out"
	conda:
		"envs/Finding_SDs.yaml"
	params:
		time = '01:00:00',
		threads = 1,
		mem = 20000,
		name = "Samtools"
	shell:
		"samtools faidx {input.masked_genome} > {log.out} 2> {log.err}"


	'''
	Find SDs in the genome using BISER
	'''
rule finding_SDs:
	input:
		indexed_genome = rules.index_genome.output.indexed_genome
	output:
		SDs = "results/finding_SDs/Branchiostoma_lanceolatum.BraLan3_SDs.bedpe",
	log:
		err = "logs/finding_SDs/biser.err",
		out = "logs/finding_SDs/biser.out"
	conda:
		"envs/Finding_SDs.yaml"
	params:
		time = '10:00:00',
		threads = 4,
		mem = 50000,
		name = "Biser"		
	shell:
		"biser -o {output.SDs} -t {params.threads} {input.indexed_genome} > {log.out} 2> {log.err}"

