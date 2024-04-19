# Set working directory
setwd("//wsl.localhost/Ubuntu/home/nanobio/AmphioxusCNVs/AmphioxusCNVs")

# Libraries 
library(dplyr)
library(tidyr)

# Ensure numbers are not printed in scientific notation
options(scipen=999)


# Load CNV files and set column names directly
cnv_F1D <- read.table('results/CNVnator/F1D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_F2D <- read.table('results/CNVnator/F2D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_F3D <- read.table('results/CNVnator/F3D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_F4D <- read.table('results/CNVnator/F4D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_F6D <- read.table('results/CNVnator/F6D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_F7D <- read.table('results/CNVnator/F7D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_F8D <- read.table('results/CNVnator/F8D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_F9D <- read.table('results/CNVnator/F9D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_F10D <- read.table('results/CNVnator/F10D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_M2D <- read.table('results/CNVnator/M2D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_M3D <- read.table('results/CNVnator/M3D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_M4D <- read.table('results/CNVnator/M4D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_M5D <- read.table('results/CNVnator/M5D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_M6D <- read.table('results/CNVnator/M6D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_M7D <- read.table('results/CNVnator/M7D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_M8D <- read.table('results/CNVnator/M8D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_M9D <- read.table('results/CNVnator/M9D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RF1D <- read.table('results/CNVnator/RF1D_cnv_calls.txt', header = FALSE, sep = "\t",
                      col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RF2D <- read.table('results/CNVnator/RF2D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RF3D <- read.table('results/CNVnator/RF3D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RF4D <- read.table('results/CNVnator/RF4D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RF5D <- read.table('results/CNVnator/RF5D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RF6D <- read.table('results/CNVnator/RF6D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RF7D <- read.table('results/CNVnator/RF7D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RF8D <- read.table('results/CNVnator/RF8D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RF10D <- read.table('results/CNVnator/RF10D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RM1D <- read.table('results/CNVnator/RM1D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RM2D <- read.table('results/CNVnator/RM2D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RM3D <- read.table('results/CNVnator/RM3D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RM4D <- read.table('results/CNVnator/RM4D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RM5D <- read.table('results/CNVnator/RM5D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RM6D <- read.table('results/CNVnator/RM6D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RM7D <- read.table('results/CNVnator/RM7D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RM8D <- read.table('results/CNVnator/RM8D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
cnv_RU5D <- read.table('results/CNVnator/RU5D_cnv_calls.txt', header = FALSE, sep = "\t",
                       col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth", "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))


# Create a list of data frames for iteration
sample_data <- list(
  cnv_F1D = cnv_F1D, cnv_F2D = cnv_F2D, cnv_F10D = cnv_F10D,
  cnv_F3D = cnv_F3D, cnv_F4D = cnv_F4D, cnv_F6D = cnv_F6D,
  cnv_F7D = cnv_F7D, cnv_F8D = cnv_F8D, cnv_F9D = cnv_F9D,
  cnv_M2D = cnv_M2D, cnv_M3D = cnv_M3D, cnv_M4D = cnv_M4D,
  cnv_M5D = cnv_M5D, cnv_M6D = cnv_M6D, cnv_M7D = cnv_M7D,
  cnv_M8D = cnv_M8D, cnv_M9D = cnv_M9D, cnv_RF1D = cnv_RF1D,
  cnv_RF2D = cnv_RF2D, cnv_RF3D = cnv_RF3D, cnv_RF4D = cnv_RF4D,
  cnv_RF5D = cnv_RF5D, cnv_RF6D = cnv_RF6D, cnv_RF7D = cnv_RF7D,
  cnv_RF8D = cnv_RF8D, cnv_RF10D = cnv_RF10D, cnv_RM1D = cnv_RM1D,
  cnv_RM2D = cnv_RM2D, cnv_RM3D = cnv_RM3D, cnv_RM4D = cnv_RM4D,
  cnv_RM5D = cnv_RM5D, cnv_RM6D = cnv_RM6D, cnv_RM7D = cnv_RM7D,
  cnv_RM8D = cnv_RM8D, cnv_RU5D = cnv_RU5D
)

# Define output filenames with the directory path
output_files <- c(
  "results/filtering_cnv/cnv_F1D.bed", "results/filtering_cnv/cnv_F2D.bed", "results/filtering_cnv/cnv_F10D.bed",
  "results/filtering_cnv/cnv_F3D.bed", "results/filtering_cnv/cnv_F4D.bed", "results/filtering_cnv/cnv_F6D.bed",
  "results/filtering_cnv/cnv_F7D.bed", "results/filtering_cnv/cnv_F8D.bed", "results/filtering_cnv/cnv_F9D.bed",
  "results/filtering_cnv/cnv_M2D.bed", "results/filtering_cnv/cnv_M3D.bed", "results/filtering_cnv/cnv_M4D.bed",
  "results/filtering_cnv/cnv_M5D.bed", "results/filtering_cnv/cnv_M6D.bed", "results/filtering_cnv/cnv_M7D.bed",
  "results/filtering_cnv/cnv_M8D.bed", "results/filtering_cnv/cnv_M9D.bed", "results/filtering_cnv/cnv_RF1D.bed",
  "results/filtering_cnv/cnv_RF2D.bed", "results/filtering_cnv/cnv_RF3D.bed", "results/filtering_cnv/cnv_RF4D.bed",
  "results/filtering_cnv/cnv_RF5D.bed", "results/filtering_cnv/cnv_RF6D.bed", "results/filtering_cnv/cnv_RF7D.bed",
  "results/filtering_cnv/cnv_RF8D.bed", "results/filtering_cnv/cnv_RF10D.bed", "results/filtering_cnv/cnv_RM1D.bed",
  "results/filtering_cnv/cnv_RM2D.bed", "results/filtering_cnv/cnv_RM3D.bed", "results/filtering_cnv/cnv_RM4D.bed",
  "results/filtering_cnv/cnv_RM5D.bed", "results/filtering_cnv/cnv_RM6D.bed", "results/filtering_cnv/cnv_RM7D.bed",
  "results/filtering_cnv/cnv_RM8D.bed", "results/filtering_cnv/cnv_RU5D.bed"
)




# Loop through each data frame in the list
for (i in names(sample_data)) {
  cnv_data <- sample_data[[i]]
  
  # Process data to convert and format fields, and sort correctly
  cnv_data <- cnv_data %>%
    select(CNV_type, Coordinates, CNV_size, Normalized_read_deapth) %>%
    separate(Coordinates, into = c("Chr", "Range"), sep = ":") %>%
    separate(Range, into = c("Start", "End"), sep = "-") %>%
    mutate(
      Start = as.numeric(format(as.numeric(Start) - 1, scientific = FALSE)),  # Adjust for BED 0-based start
      End = as.numeric(format(as.numeric(End), scientific = FALSE)),
      ChrNum = as.numeric(gsub("chr", "", Chr))  # Extract numeric part for sorting
    ) %>%
    arrange(ChrNum, Start) %>%
    select(Chr, Start, End, CNV_type, Normalized_read_deapth)  # Ensure correct column order for BED
  
  # Save as BED file directly in the specified directory
  write.table(cnv_data, output_files[which(names(sample_data) == i)], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}


  
## MULTIINTER WITH BEDTOOLS 


bedtools multiinter -i /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/cnv_F1D_adjusted.bed /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/cnv_F2D_adjusted.bed /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/cnv_F10D_adjusted.bed > /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/multiIntersectOutput.bed

bedtools intersect -a /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/multiIntersectOutput.bed -b /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/cnv_F1D.bed -wao > /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/intersect_F1D.bed

bedtools intersect -a /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/multiIntersectOutput.bed -b /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/cnv_F2D.bed -wao > /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/intersect_F2D.bed

bedtools intersect -a /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/multiIntersectOutput.bed -b /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/cnv_F10D.bed -wao > /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/intersect_F10D.bed
bedtools intersect -a /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/multiIntersectOutput.bed -b /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/cnv_F1D_adjusted.bed -wao > /home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/filtering_cnv/intersect_F1D.bed

cat cnv_F2D.bed | sort -k1,2V | awk '{if(end==$2){st=$2+1;print $1"\t"st"\t"$3"\t"$4"\t"$5}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5}end=$3}' > cnv_F2D_adjusted.bed

paste <(cat intersect_F1D.bed | cut -f1,2,3,12,13) <(cat intersect_F2D.bed | cut -f12,13)  <(cat intersect_F10D.bed | cut -f12,13) > paste_all_cnv.bed




'''
	Merge all CNV sample bed files with bedtools 
	'''
rule merge_cnv:
  input:
  cnv_bed=rules.transform_txt_into_bed.output.cnv_bed
output:
  adjusted_bed="results/filtering_cnv/cnv_{sample}_adjusted.bed.bed",
multiinter="results/filtering_cnv/multiIntersectOutput.bed",
intersection_bed="results/filtering_cnv/intersect_{sample}.bed",
all_samples_CNVs="results/filtering_cnv/paste_all_cnv.bed"
conda:
  "envs/Detecting_CNVs.yaml"
shell:
  """
		cat {input.cnv_bed} | sort -k1,2V | awk '{if(end==$2){st=$2+1;print $1"\t"st"\t"$3"\t"$4"\t"$5}else{print $1"\t"$2"\t"$3"\t"$4"\t"$5}end=$3}' > {output.adjusted_bed}
	bedtools multiinter -i results/filtering_cnv/cnv_F1D_adjusted.bed results/filtering_cnv/cnv_F2D_adjusted.bed results/filtering_cnv/cnv_F10D_adjusted.bed > {output.multiinter}
		bedtools intersect -a {output.multiinter} -b {output.adjusted_bed} -wao > {output.intersection_bed}
	paste <(cat intersect_F1D.bed | cut -f1,2,3,12,13) <(cat intersect_F2D.bed | cut -f12,13)  <(cat intersect_F10D.bed | cut -f12,13) > {output.all_samples_CNVs}
		"""



