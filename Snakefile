	'''
	Mask repeats in the genome using Tandem repeats Finder
	'''
rule mask_tandem_repeats_with_trf:
	input:
		amphioxus_genome = "data/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		genome_trf = "results/mask_tandem_repeats_with_trf/Branchiostoma_lanceolatum.BraLan3_genome.fa.2.7.7.80.10.50.15.mask",
		dat_file = "results/mask_tandem_repeats_with_trf/Branchiostoma_lanceolatum.BraLan3_genome.fa.2.7.7.80.10.50.15.dat",
		trf_repeats = "results/mask_tandem_repeats_with_trf/repeats_trf.txt"
	log:
		err = "logs/mask_tandem_repeats_with_trf/trf.err",
		out = "logs/mask_tandem_repeats_with_trf/trf.out"
	conda:
		"envs/Finding_SDs.yaml"
	params:
		name = "Tandem_Repeat_Finder",
		time = '10:00:00',
		threads = 2,
		mem = 20000
	shell:
		"""
		set +euo pipefail
		pwd=$(pwd)
		cd $(dirname {output.genome_trf})
		trf ${{pwd}}/{input.amphioxus_genome} 2 7 7 80 10 50 15 -l 10 -h -m > ${{pwd}}/{log.out} 2> ${{pwd}}/{log.err}
		grep -v "^$" Branchiostoma_lanceolatum.BraLan3_genome.fa.2.7.7.80.10.50.15.dat | awk '/Sequence:/ {{chromosome=$2}} {{print chromosome, $1, $2}}' > repeats_info.txt
		tail -n +6 repeats_info.txt | grep -v "Sequence" | grep -v "Parameters" |grep -v "scaf" > repeats_trf.txt
		cd ${{pwd}}
		"""



	'''
	Mask repeats in the genome using RepeatMasker
	'''
rule mask_genome_with_repeatmasker:
	input:
		genome_trf = rules.mask_tandem_repeats_with_trf.output.genome_trf
	output:
		masked_genome = "results/mask_genome_with_repeatmasker/Branchiostoma_lanceolatum.BraLan3_genome.fa.2.7.7.80.10.50.15.mask.masked",
		repeatmasker_out = "results/mask_genome_with_repeatmasker/Branchiostoma_lanceolatum.BraLan3_genome.fa.2.7.7.80.10.50.15.mask.out",
		repeatmasker_repeats = "results/mask_genome_with_repeatmasker/repeats_repeatmasker.txt"
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
		"""
		RepeatMasker -s -xsmall -e ncbi {input.genome_trf} -dir $(dirname {output.masked_genome}) > {log.out} 2> {log.err}
		awk '{{if (NR > 1) print $5, $6, $7}}' results/mask_genome_with_repeatmasker/Branchiostoma_lanceolatum.BraLan3_genome.fa.2.7.7.80.10.50.15.mask.out | tail -n +3 | grep -v "scaf" > results/mask_genome_with_repeatmasker/repeats_repeatmasker.txt
		"""


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


	'''
	Making directories
	'''
rule make_directories:
	shell:
		"""
		mkdir -p results/filtered_sd_data
		mkdir -p results/plots
		"""


	'''
	Filtering with R 
	'''
rule filter_sd_for_bedtools:
	input:
		SDs = rules.finding_SDs.output.SDs
	output:
		sorted_bed="results/filtered_sd_data/sd_positions_sorted.bed",
		intra_bed="results/filtered_sd_data/sd_positions_intra.bed",
		inter_bed="results/filtered_sd_data/sd_positions_inter.bed",
		length_plot = "results/plots/Histogram_of_filtered_SD_length.pdf"
	conda:
		"envs/Finding_SDs.yaml"
	script:
		"scripts/Filtering_SDs.R"


	'''
	Merge with bedtools
	'''
rule merge_with_bedtools:
	input:
		sorted_bed = rules.filter_sd_for_bedtools.output.sorted_bed,
		intra_bed = rules.filter_sd_for_bedtools.output.intra_bed,
		inter_bed = rules.filter_sd_for_bedtools.output.inter_bed
	output:
		mixed_cases="results/filtered_sd_data/mixed_sd_cases.bed",
		pure_intra="results/filtered_sd_data/pure_intra_cases.bed",
		pure_inter="results/filtered_sd_data/pure_inter_cases.bed",
		merged_bed="results/filtered_sd_data/sd_positions_merged.bed",
	conda:
		"envs/Finding_SDs.yaml"
	params:
		sd_positions_inter_merged="results/filtered_sd_data/sd_positions_inter_merged.bed",
		sd_positions_intra_merged="results/filtered_sd_data/sd_positions_intra_merged.bed"
	shell:
		"""
		bedtools merge -i {input.intra_bed} > {params.sd_positions_intra_merged}
		bedtools merge -i {input.inter_bed} > {params.sd_positions_inter_merged}
		bedtools intersect -a {param.sd_positions_intra_merged} -b {param.sd_positions_inter_merged} > {output.mixed_cases}
		bedtools intersect -v -a {param.sd_positions_intra_merged} -b {param.sd_positions_inter_merged} > {output.pure_intra}
		bedtools intersect -v -a {param.sd_positions_inter_merged} -b {param.sd_positions_intra_merged} > {output.pure_inter}
		bedtools merge -i {input.sorted_bed} > {output.merged_bed}
		"""


	'''
	Plot with R 
	'''
rule plot_SDs:
	input:
		mixed_cases= rules.merge_with_bedtools.output.mixed_cases,
		pure_intra= rules.merge_with_bedtools.output.pure_intra,
		pure_inter= rules.merge_with_bedtools.output.pure_inter,
		merged_bed= rules.merge_with_bedtools.output.merged_bed,
		trf_repeats= rules.mask_tandem_repeats_with_trf.output.trf_repeats,
		repeatmasker_repeats= rules.mask_genome_with_repeatmasker.output.repeatmasker_repeats
	output:
		sd_bar_plot="results/plots/sd_bar_plot.pdf",
		sd_chr_plot1="results/plots/sd_chr_plot1.pdf",
		sd_chr_plot2="results/plots/sd_chr_plot2.pdf",
		sd_chr_plot3="results/plots/sd_chr_plot3.pdf"
	conda:
		"envs/Finding_SDs.yaml"
	script:
		"scripts/Plot_SDs.R"





configfile: "config.yaml"

rule run_all_samples:
    input:
        expand("results/BAM_Merging/{sample}_merged.bam", sample=config['samples'])
    output:
        "test"
    shell: 
        "echo test > {output}" 


rule Merge_BAM_Files_PerSample:
	'''
	Merge multiple BAM files for each sample into a single BAM file.
	'''
	input:
		bamFiles = expand("data/{{sample}}{combo}_sorted_markdup.bam", combo=config["combos"])
	output:
		mergedBAM = "results/BAM_Merging/{sample}_merged.bam"
	log:
		err = "logs/BAM_Merging/{sample}_merge.err",
		out = "logs/BAM_Merging/{sample}_merge.out"
	benchmark:
		"benchmarks/BAM_Merging/{sample}_merge.txt"
	conda:
		"envs/Detecting_CNVs.yaml"
	params:
		time = '02:00:00',
		name = "MergeBAM{sample}",
		threads = 4,
		mem = 16000
	shell:
		"""
		mkdir -p $(dirname {output.mergedBAM})
		samtools merge -@ {params.threads} {output.mergedBAM} {input.bamFiles} > {log.out} 2> {log.err}
		"""





	'''
	Cleaning the reference genome
	'''
rule reference_genome_clean:
	input:
		amphioxus_genome = "data/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		amphioxus_genome_cleaned = "data/Cleaned_Branchiostoma_lanceolatum.BraLan3_genome.fa",
		configfile_ref="results/CNVpytor/BraLan3_conf.py",
	log:
		err = "logs/clean_reference_genome/cleaning.err",
		out = "logs/clean_reference_genome/cleaning.out"
	conda:
		"envs/Detecting_CNVs.yaml",
	shell:
		"""
		cat {input.amphioxus_genome} | awk 'BEGIN {{p=1}} /^>/{{if ($0 ~ /scaf/) p=0; else p=1}} p' > {output.amphioxus_genome_cleaned}
		samtools faidx {output.amphioxus_genome_cleaned}
		echo 'from collections import OrderedDict' > {output.configfile_ref}
 		echo 'import_reference_genomes = {{' >> {output.configfile_ref}
		echo '    "Branchiostoma": {{' >> {output.configfile_ref}
		echo '        "name": "BraLan3",' >> {output.configfile_ref}
		echo '        "species": "Branchiostoma lanceolatum",' >> {output.configfile_ref}
		echo '        "chromosomes": OrderedDict([' >> {output.configfile_ref}
		echo '            ("chr1", (43860960, "A")),' >> {output.configfile_ref}
		echo '            ("chr2", (38510819, "A")),' >> {output.configfile_ref}
		echo '            ("chr3", (34610492, "A")),' >> {output.configfile_ref}
		echo '            ("chr4", (31719604, "A")),' >> {output.configfile_ref}
		echo '            ("chr5", (25701974, "A")),' >> {output.configfile_ref}
		echo '            ("chr6", (24533633, "A")),' >> {output.configfile_ref}
		echo '            ("chr7", (24230189, "A")),' >> {output.configfile_ref}
		echo '            ("chr8", (23752511, "A")),' >> {output.configfile_ref}
		echo '            ("chr9", (23231292, "A")),' >> {output.configfile_ref}
		echo '            ("chr10", (20381850, "A")),' >> {output.configfile_ref}
		echo '            ("chr11", (20367708, "A")),' >> {output.configfile_ref}
		echo '            ("chr12", (19917020, "A")),' >> {output.configfile_ref}
		echo '            ("chr13", (19776172, "A")),' >> {output.configfile_ref}
		echo '            ("chr14", (19709165, "A")),' >> {output.configfile_ref}
		echo '            ("chr15", (19381563, "A")),' >> {output.configfile_ref}
		echo '            ("chr16", (18823661, "A")),' >> {output.configfile_ref}
		echo '            ("chr17", (18214296, "A")),' >> {output.configfile_ref}
		echo '            ("chr18", (17113871, "A")),' >> {output.configfile_ref}
		echo '            ("chr19", (15322015, "A"))' >> {output.configfile_ref}
		echo '        ]),' >> {output.configfile_ref}
		echo '        "gc_file": "results/CNVpytor/BraLan3_gc.pytor"' >> {output.configfile_ref}
		echo '    }}' >> {output.configfile_ref}
		echo '}}' >> {output.configfile_ref}
		"""



configfile: "config.yaml"

rule run_all_samples_for_CNVs_pytor:
	input:
		expand("results/CNVpytor/{sample}_cnv_calls.tsv", sample=config['samples'])
	output:
		"cnvpytor"
	shell: 
		"echo cnvpytor > {output}" 

rule run_cnvpytor:
	'''
	Call CNVs with cnvpytor
	'''
	input:
		bam=rules.Merge_BAM_Files_PerSample.output.mergedBAM,
		configfile_ref=rules.reference_genome_clean.output.configfile_ref,
		amphioxus_genome_cleaned=rules.reference_genome_clean.output.amphioxus_genome_cleaned
	output:
		cnv_calls="results/CNVpytor/{sample}_cnv_calls.tsv"
	log:
		err="logs/CNVpytor/{sample}_cnvpytor.err",
		out="logs/CNVpytor/{sample}_cnvpytor.out"
	conda:
		"envs/Detecting_CNVs.yaml"
	params:
		pytor_file="results/CNVpytor/{sample}.pytor",
		gc_ref="results/CNVpytor/BraLan3_gc.pytor",
		bin_size=1000
	shell:
		"""
		cnvpytor -root {params.gc_ref} -gc {input.amphioxus_genome_cleaned} -make_gc_file
		cnvpytor -root {params.pytor_file} -rd {input.bam} -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 -conf input.configfile_ref} > {log.out} 2> {log.err} 
		cnvpytor -root {params.pytor_file} -his {params.bin_size} -conf {input.configfile_ref} > {log.out} 2> {log.err} 
		cnvpytor -root {params.pytor_file} -partition {params.bin_size} -conf {input.configfile_ref} > {log.out} 2> {log.err}
		cnvpytor -root {params.pytor_file} -call {params.bin_size} -conf {input.configfile_ref} > {log.out} 2> {log.err}
		"""
