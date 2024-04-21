# Set working directory
#setwd("//wsl.localhost/Ubuntu/home/nanobio/AmphioxusCNVs/AmphioxusCNVs")

# Libraries
library(dplyr)
library(tidyr)

# Ensure numbers are not printed in scientific notation
options(scipen=999)

# Function to load data safely
safe_read_table <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  read.table(file_path, header = FALSE, sep = "\t",
             col.names = c("CNV_type", "Coordinates", "CNV_size", "Normalized_read_deapth",
                           "Evalue1", "Evalue2", "Evalue3", "Evalue4", "q0"))
}

# Define file names and paths
sample_names <- c("F1D", "F2D", "F3D", "F4D", "F6D", "F7D", "F8D", "F9D", "F10D",
                  "M2D", "M3D", "M4D", "M5D", "M6D", "M7D", "M8D", "M9D",
                  "RF1D", "RF2D", "RF3D", "RF4D", "RF5D", "RF6D", "RF7D", "RF8D", "RF10D",
                  "RM1D", "RM2D", "RM3D", "RM4D", "RM5D", "RM6D", "RM7D", "RM8D", "RU5D")
sample_data <- list()

# Load data for each sample
for (sample in sample_names) {
  file_path <- paste0('results/CNVnator/', sample, '_cnv_calls.txt')
  sample_data[[sample]] <- safe_read_table(file_path)
}

# Define output filenames with the directory path
output_files <- paste0("results/filtering_cnv/cnv_", sample_names, ".bed")

# Process and save data for each sample
for (i in seq_along(sample_data)) {
  cnv_data <- sample_data[[i]]
  
  # Data processing
  cnv_data <- cnv_data %>%
    select(CNV_type, Coordinates, CNV_size, Normalized_read_deapth) %>%
    separate(Coordinates, into = c("Chr", "Range"), sep = ":") %>%
    separate(Range, into = c("Start", "End"), sep = "-") %>%
    mutate(Start = as.integer(as.numeric(Start) - 1),  # Adjust for BED 0-based start
           End = as.integer(as.numeric(End)),
           ChrNum = as.numeric(gsub("chr", "", Chr))) %>%
    arrange(ChrNum, Start) %>%
    select(Chr, Start, End, CNV_type, Normalized_read_deapth)
  
  # Save as BED file
  write.table(cnv_data, output_files[i], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
