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
		set +euo pipefail
		mkdir -p $(dirname {output.mergedBAM})
		samtools merge -@ {params.threads} {output.mergedBAM} {input.bamFiles} > {log.out} 2> {log.err}
		"""


	'''
	Split reference genome into chromosomes
	'''
rule split_reference_genome:
	input:
		ref_genome="data/Branchiostoma_lanceolatum.BraLan3_genome.fa"
	output:
		"data/chromosomes/.done"
	shell:
		"""
		mkdir -p data/chromosomes/
		cat {input.ref_genome} | awk '{if($1 ~ />/){if($1 ~ />chr/){out="data/chromosomes/"substr($1,2)".fa";x++}if(x>0){close(out)}}else{print $0 > out}}'
		touch {output}
		"""





configfile: "config.yaml"

rule run_all_samples_for_CNVs:
	input:
		expand("results/CNVnator/{sample}_cnv_calls.txt", sample=config['samples'])
	output:
		"cnvnator"
	shell: 
		"echo cnvnator > {output}" 

rule run_cnvnator:
	'''
	Call CNVs with cnvnator
	'''
	input:
		bam=rules.Merge_BAM_Files_PerSample.output.mergedBAM,
		chrom_list="data/chromosomes/chrom_list.txt"
	output:
		cnv_calls="results/CNVnator/{sample}_cnv_calls.txt"
	log:
		err="logs/CNVnator/{sample}_cnvnator.err",
		out="logs/CNVnator/{sample}_cnvnator.out"
	conda:
		"envs/Detecting_CNVs.yaml"
	params:
		root_file="results/CNVnator/{sample}.root",
		merged_root="results/CNVnator/all_samples.root",
		final_merged_root="results/CNVnator/merged_samples.root",
		bin_size=100,
		ref_genome_dir="data/chromosomes/"
	shell:
		"""
		cnvnator -root {params.root_file} -tree {input.bam} -chrom $(cat {input.chrom_list}) > {log.out} 2> {log.err} 
		cnvnator -root {params.root_file} -his {params.bin_size} -d {params.ref_genome_dir} -chrom $(cat {input.chrom_list}) > {log.out} 2> {log.err} 
		cnvnator -root {params.root_file} -stat {params.bin_size} -chrom $(cat {input.chrom_list}) > {log.out} 2> {log.err} 
		cnvnator -root {params.root_file} -partition {params.bin_size} -chrom $(cat {input.chrom_list}) > {log.out} 2> {log.err}
		cnvnator -root {params.root_file} -call {params.bin_size} -chrom $(cat {input.chrom_list}) > {output.cnv_calls} 2> {log.err} 
		"""



configfile: "config.yaml"

rule all_samples_transform_txt_into_bed:
	input:
		expand("results/filtering_cnv/cnv_{sample}.bed", sample=config['samples'])
	output:
		"rtrans"
	shell: 
		"echo rtrans > {output}" 

	'''
	Transform output cnvnator into bed files
	'''
rule transform_txt_into_bed:
	input:
		expand("results/CNVnator/{sample}_cnv_calls.txt", sample=config['samples'])
	output:
		expand("results/filtering_cnv/cnv_{sample}.bed", sample=config['samples'])
	conda:
		"envs/Detecting_CNVs.yaml"
	script:
		"scripts/Exploring_CNVs.R"




configfile: "config.yaml"

rule all_samples_merge_cnv:
	input:
		expand("results/filtering_cnv/intersect_{sample}.bed", sample=config['samples'])
	output:
		"merge"
	shell: 
		"echo merge > {output}" 

	'''
	Merge all CNV sample bed files with bedtools 
	'''
rule filtering_cnv:
	input:
		cnv_bed=rules.transform_txt_into_bed.output.cnv_bed
	output:
		adjusted_bed="results/filtering_cnv/cnv_{sample}_adjusted.bed.bed",
		intersection_bed="results/filtering_cnv/intersect_{sample}.bed"
	log:
		err="logs/filtering_cnv/{sample}_merging.err",
		out="logs/filtering_cnv/{sample}_merging.out"
	conda:
		"envs/Detecting_CNVs.yaml"
	params:
		multiinter="results/filtering_cnv/multiIntersectOutput.bed",
		all_samples_CNVs="results/filtering_cnv/paste_all_cnv.bed"
	shell:
		"""
		cat {input.cnv_bed} | sort -k1,2V | awk '{{if(end==$2){{st=$2+1;print $1" "st" "$3" "$4" "$5}}else{{print $1" "$2" "$3" "$4" "$5}}end=$3}}' > {output.adjusted_bed}
		bedtools multiinter -i results/filtering_cnv/cnv_F1D_adjusted.bed.bed results/filtering_cnv/cnv_F2D_adjusted.bed.bed results/filtering_cnv/cnv_F3D_adjusted.bed.bed results/filtering_cnv/cnv_F4D_adjusted.bed.bed results/filtering_cnv/cnv_F6D_adjusted.bed.bed results/filtering_cnv/cnv_F7D_adjusted.bed.bed results/filtering_cnv/cnv_F8D_adjusted.bed.bed results/filtering_cnv/cnv_F9D_adjusted.bed.bed results/filtering_cnv/cnv_F10D_adjusted.bed.bed results/filtering_cnv/cnv_M2D_adjusted.bed.bed results/filtering_cnv/cnv_M3D_adjusted.bed.bed results/filtering_cnv/cnv_M4D_adjusted.bed.bed results/filtering_cnv/cnv_M5D_adjusted.bed.bed results/filtering_cnv/cnv_M6D_adjusted.bed.bed results/filtering_cnv/cnv_M7D_adjusted.bed.bed results/filtering_cnv/cnv_M8D_adjusted.bed.bed results/filtering_cnv/cnv_M9D_adjusted.bed.bed results/filtering_cnv/cnv_RF1D_adjusted.bed.bed results/filtering_cnv/cnv_RF2D_adjusted.bed.bed results/filtering_cnv/cnv_RF3D_adjusted.bed.bed results/filtering_cnv/cnv_RF4D_adjusted.bed.bed results/filtering_cnv/cnv_RF5D_adjusted.bed.bed results/filtering_cnv/cnv_RF6D_adjusted.bed.bed results/filtering_cnv/cnv_RF7D_adjusted.bed.bed results/filtering_cnv/cnv_RF8D_adjusted.bed.bed results/filtering_cnv/cnv_RF10D_adjusted.bed.bed results/filtering_cnv/cnv_RM1D_adjusted.bed.bed results/filtering_cnv/cnv_RM2D_adjusted.bed.bed results/filtering_cnv/cnv_RM3D_adjusted.bed.bed results/filtering_cnv/cnv_RM4D_adjusted.bed.bed results/filtering_cnv/cnv_RM5D_adjusted.bed.bed results/filtering_cnv/cnv_RM6D_adjusted.bed.bed v results/filtering_cnv/cnv_RM7D_adjusted.bed.bed results/filtering_cnv/cnv_RM8D_adjusted.bed.bed results/filtering_cnv/cnv_RU5D_adjusted.bed.bed > {params.multiinter}
		bedtools intersect -a {params.multiinter} -b {output.adjusted_bed} -wao > {output.intersection_bed}
		paste <(cat intersect_F1D.bed | cut -f1,2,3,12,13) \
		<(cat intersect_F2D.bed | cut -f12,13) \
		<(cat intersect_F3D.bed | cut -f12,13) \
		<(cat intersect_F4D.bed | cut -f12,13) \
		<(cat intersect_F6D.bed | cut -f12,13) \
		<(cat intersect_F7D.bed | cut -f12,13) \
		<(cat intersect_F8D.bed | cut -f12,13) \
		<(cat intersect_F9D.bed | cut -f12,13) \
		<(cat intersect_F10D.bed | cut -f12,13) \
		<(cat intersect_M2D.bed | cut -f12,13) \
		<(cat intersect_M3D.bed | cut -f12,13) \
		<(cat intersect_M4D.bed | cut -f12,13) \
		<(cat intersect_M5D.bed | cut -f12,13) \
		<(cat intersect_M6D.bed | cut -f12,13) \
		<(cat intersect_M7D.bed | cut -f12,13) \
		<(cat intersect_M8D.bed | cut -f12,13) \
		<(cat intersect_M9D.bed | cut -f12,13) \
		<(cat intersect_RF1D.bed | cut -f12,13) \
		<(cat intersect_RF2D.bed | cut -f12,13) \
		<(cat intersect_RF3D.bed | cut -f12,13) \
		<(cat intersect_RF4D.bed | cut -f12,13) \
		<(cat intersect_RF5D.bed | cut -f12,13) \
		<(cat intersect_RF6D.bed | cut -f12,13) \
		<(cat intersect_RF7D.bed | cut -f12,13) \
		<(cat intersect_RF8D.bed | cut -f12,13) \
		<(cat intersect_RF10D.bed | cut -f12,13) \
		<(cat intersect_RM1D.bed | cut -f12,13) \
		<(cat intersect_RM2D.bed | cut -f12,13) \
		<(cat intersect_RM3D.bed | cut -f12,13) \
		<(cat intersect_RM4D.bed | cut -f12,13) \
		<(cat intersect_RM5D.bed | cut -f12,13) \
		<(cat intersect_RM6D.bed | cut -f12,13) \
		<(cat intersect_RM7D.bed | cut -f12,13) \
		<(cat intersect_RM8D.bed | cut -f12,13) \
		<(cat intersect_RU5D.bed | cut -f12,13) > {params.all_samples_CNVs}
		"""









