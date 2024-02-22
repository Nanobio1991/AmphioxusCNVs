	'''
	Mask repeats in the genome using RepeatMasker
	'''
rule mask_genome_with_repeatmasker:
	input:
		amphioxus_genome = "data/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		masked_genome = "results/mask_genome_with_repeatmasker/Branchiostoma_lanceolatum.BraLan3_genome.fa.masked"
	log:
		err = "logs/mask_genome_with_repeatmasker/RepeatMasker.err",
		out = "logs/mask_genome_with_repeatmasker/RepeatMasker.out"
	conda:
		'envs/Finding_SDs.yaml'
	params:
		output_dir = "results/mask_genome_with_repeatmasker/",
		name = "RepeatMasker",
	shell:
		"RepeatMasker -s -xsmall -e ncbi {input.amphioxus_genome} -dir $(dirname {output.masked_genome}) >> {log.out} 2> {log.err}"


	'''
	Mask repeats in the genome using Tandem repeats Finder
	'''
rule mask_tandem_repeats_with_trf:
	input:
		masked_genome = rules.mask_genome_with_repeatmasker.output.masked_genome
	output :
		masked_genome_with_trf = "results/mask_tandem_repeats_with_trf/Branchiostoma_lanceolatum.BraLan3_masked_trf.fa" 
	log:
		err = "logs/mask_tandem_repeats_with_trf/trf.err",
		out = "logs/mask_tandem_repeats_with_trf/trf.out"
	conda:
		'envs/Finding_SDs.yaml'
	params:
		name = "Tandem_Repeat_Finder"
	shell:
		"trf {input.masked_genome} 2 7 7 80 10 50 500 -l 25 -h -ngs > {output.masked_genome_with_trf} >> {log.out} 2> {log.err}"


	'''
	Indexing genome with samtools
	'''
rule index_genome:
    input:
        masked_genome_with_trf = rules.mask_tandem_repeats_with_trf.output.masked_genome_with_trf
    output:
        indexed_genome = "results/mask_tandem_repeats_with_trf/Branchiostoma_lanceolatum.BraLan3_masked_trf.fa.fai"
#    log:
#		err = "logs/index_genome/index.err",
#		out = "logs/index_genome/index.out"
    conda:
        "envs/Finding_SDs.yaml"
    shell:
        "samtools faidx {input.masked_genome_with_trf}"


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
        threads = "4",  # Adjust
    shell:
        "biser -o {output.SDs} -t {params.threads} {input.indexed_genome} >> {log.out} 2> {log.err}"

