	'''
	Mask repeats in the genome using Tandem repeats Finder
	'''
rule mask_tandem_repeats_with_trf:
	input:
		amphioxus_genome = "data/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output :
		genome_trf = "results/mask_tandem_repeats_with_trf/Branchiostoma_lanceolatum.BraLan3_genome.fa.2.7.7.80.10.50.15.mask",
	log:
		err = "logs/mask_tandem_repeats_with_trf/trf.err",
		out = "logs/mask_tandem_repeats_with_trf/trf.out"
	conda:
		"envs/Finding_SDs.yaml"
	params:
		name = "Tandem_Repeat_Finder",
		time = '10:00:00',
		threads = 2,
		mem = 20000,
	shell:
		"""
		pwd=$(pwd)
		cd $(dirname {output.genome_trf})
		trf ${{pwd}}/{input.amphioxus_genome} 2 7 7 80 10 50 15 -l 10 -h -m > ${{pwd}}/{log.out} 2> ${{pwd}}/{log.err}
		cd ${{pwd}}
		"""

	'''
	Mask repeats in the genome using RepeatMasker
	'''
rule mask_genome_with_repeatmasker:
	input:
		genome_trf = rules.mask_tandem_repeats_with_trf.output.genome_trf
	output:
		masked_genome = "results/mask_genome_with_repeatmasker/Branchiostoma_lanceolatum.BraLan3_genome.fa.2.7.7.80.10.50.15.mask.masked"
	log:
		err = "logs/mask_genome_with_repeatmasker/RepeatMasker.err",
		out = "logs/mask_genome_with_repeatmasker/RepeatMasker.out"
	conda:
		"envs/Finding_SDs.yaml"
	params:
		output_dir = "results/mask_genome_with_repeatmasker/",
		time = '10:00:00',
		threads = 2,
		mem = 20000,
		name = "RepeatMasker",
	shell:
		"RepeatMasker -s -xsmall -e ncbi {input.genome_trf} -dir $(dirname {output.masked_genome}) > {log.out} 2> {log.err}"


	'''
	Replace lower case bases with N
	'''
rule replace_bases_with_N:
	input:
		masked_genome = rules.mask_genome_with_repeatmasker.output.masked_genome
	output:
		maskedN_genome = "results/masked_genome/Branchiostoma_lanceolatum.BraLan3_genome.fa.masked"
	log:
		err = "logs/replace/replace.err",
		out = "logs/replace/replace.out"
	conda:
		"envs/Finding_SDs.yaml"
	params:
		time = '01:00:00',
		threads = 1,
		mem = 20000,
		name = "Samtools"
	shell:
		"sed 's/[atgc]/N/g' {input.masked_genome} > {output.maskedN_genome}"


	'''
	Indexing genome with samtools
	'''
rule index_genome:
	input:
		maskedN_genome = rules.replace_bases_with_N.output.maskedN_genome
	output:
		indexed_genome = "results/masked_genome/Branchiostoma_lanceolatum.BraLan3_genome.fa.masked.fai"
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
		"samtools faidx {input.maskedN_genome} > {log.out} 2> {log.err}"


	'''
	Find SDs in the genome using BISER
	'''
rule finding_SDs:
	input:
		maskedN_genome = rules.replace_bases_with_N.output.maskedN_genome,
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
		"biser -o {output.SDs} -t {params.threads} --gc-heap 1G --hard {input.maskedN_genome} > {log.out} 2> {log.err}"
