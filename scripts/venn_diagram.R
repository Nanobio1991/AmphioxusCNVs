install.packages("openxlsx")

library(openxlsx)
library(dplyr)
library(tidyr)



multiinter <- read.table('results/plots/venn/multiintersect_exons.bed', header = FALSE, sep = "\t")


multiinter <- multiinter %>%
  select(Chr = V1, st = V2, end = V3, num = V4, list = V5, cnvs = V6, sds = V7, exons = V8)


multiinter <- multiinter %>%
  mutate(size = end-st)



# Convert numeric 0/1 to logical TRUE/FALSE for Venn diagram
multiinter <- multiinter %>%
  mutate(
    cnvs = as.logical(cnvs),
    exons = as.logical(exons),
    sds = as.logical(sds)
  )

# Expand the dataframe by replicating rows according to the size
multiinter_expanded <- multiinter %>%
  uncount(size, .id = "replicate")  




summary_counts <- multiinter_expanded %>%
  summarise(
    Pure_CNVs = sum(cnvs == TRUE & sds == FALSE & exons == FALSE),
    Pure_SDs = sum(cnvs == FALSE & sds == TRUE & exons == FALSE),
    Pure_Exons = sum(cnvs == FALSE & sds == FALSE & exons == TRUE),
    CNVs_and_SDs = sum(cnvs == TRUE & sds == TRUE & exons == FALSE),
    SDs_and_Exons = sum(cnvs == FALSE & sds == TRUE & exons == TRUE),
    CNVs_and_Exons = sum(cnvs == TRUE & sds == FALSE & exons == TRUE),
    CNVs_SDs_and_Exons = sum(cnvs == TRUE & sds == TRUE & exons == TRUE)
  )




print(summary_counts)

# Write the dataframe to an Excel file
write.xlsx(summary_counts, file = "results/plots/venn/counts_for_venn_diagram.xlsx")




###############################################################################
############ PLOT GRAPH FIXED AND NON-FIXED SDS ###############################
###############################################################################

## CREATE DATA FRAME FOR THE LENGTH OF CHROMOSOMES ################

length_chr <- read.table('data/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt', header = FALSE, sep = "\t")
head(length_chr)

# Ensure lenDF is prepared correctly
length_chr <- length_chr %>%
  rename(Chr = V1, Length = V2)

length_chr <- length_chr %>%
  filter(!grepl("^>scaf", Chr))

length_chr$Chr <- sub(">", "", length_chr$Chr)  # Remove ">" from chromosome names



############################################
fixed_sds <- read.table('results/plots/venn/fixed_sds.bed', header = FALSE, sep = "\t")
non_fixed_sds <- read.table('results/plots/venn/non_fixed_sds.bed', header = FALSE, sep = "\t")
##############################################

fixed_sds <- fixed_sds %>%
  rename(Chr = V1, st = V2, end = V3)

non_fixed_sds <- non_fixed_sds %>%
  rename(Chr = V1, st = V2, end = V3)

plot_all_chr_rows_len <- function(lenDF, regDF, color){
  height <- 0.6
  maxY <- length(lenDF$Chr) * (height + 0.2) # Adjust spacing
  maxChrLength <- max(lenDF$Length) # Find the maximum chromosome length
  
  # Adjust plot to reflect chromosome lengths properly
  plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, maxY), xlim=c(0, maxChrLength), col=NA)
  
  for(c in 1:length(lenDF$Chr)){
    chrName <- lenDF$Chr[c]
    print(chrName)
    pos <- maxY - (c * (height + 0.2)) # Adjust for spacing
    chromosomeLength <- lenDF$Length[c]
    
    # Filter regDF for the current chromosome
    currentRegDF <- regDF[regDF$Chr == chrName,]
    sts <- currentRegDF$st / chromosomeLength * chromosomeLength
    ends <- currentRegDF$end / chromosomeLength * chromosomeLength
    
    # Draw the chromosome and the SDs
    draw_chr_row(height, pos, sts, ends, color)
    polygon(c(0, chromosomeLength, chromosomeLength, 0), c(pos-height/2, pos-height/2, pos+height/2, pos+height/2), col=NA, border="black")
  }
  
  # Create descriptive labels for each chromosome
  chrLabels <- lenDF$Chr
  
  # Add y-axis with descriptive labels
  axis(2, at = seq(0, maxY - (height + 0.2), by = (height + 0.2)), labels=rev(chrLabels), las=1, cex.axis=0.7)
}

draw_chr_row <- function(h, p, starts, ends, color){
  for(r in 1:length(starts)){
    polygon(c(starts[r], ends[r], ends[r], starts[r]), c(p-h/2, p-h/2, p+h/2, p+h/2), col=color, border=NA)
  }
}


pdf("results/plots/fixed_sds.pdf")

plot_all_chr_rows_len(length_chr, fixed_sds, "blue")

dev.off()



pdf("results/plots/non_fixed_sds.pdf")

plot_all_chr_rows_len(length_chr, non_fixed_sds, "blue")

dev.off()






# Load necessary libraries
library(dplyr)


# Function to perform KS test within a chromosome
ks_test_within_chromosome <- function(df, chrom_length) {
  midpoints <- (df$st + df$end) / 2
  normalized_midpoints <- midpoints / chrom_length
  ks_test <- ks.test(normalized_midpoints, "punif")
  return(ks_test)
}

# Apply the KS test to each chromosome
ks_results <- sd_positions_merged %>%
  group_by(Chr) %>%
  do({
    chrom_length <- length_chr$Length[length_chr$Chr == unique(.$Chr)]
    ks_test <- ks_test_within_chromosome(., chrom_length)
    data.frame(Chr = unique(.$Chr), D = ks_test$statistic, p.value = ks_test$p.value)
  })

# Print KS test results
ks_results


fixed_sds <- read.table('results/plots/venn/true_sds.bed', header = FALSE, sep = "\t")


non_fixed_sds <- read.table('results/plots/venn/nonfixed_sds.bed', header = FALSE, sep = "\t")
fixed_sd <- fixed_sd %>%
  rename(Chr = V1, st = V2, end = V3)


# Apply the KS test to each chromosome
ks_results <- fixed_sd %>%
  group_by(Chr) %>%
  do({
    chrom_length <- length_chr$Length[length_chr$Chr == unique(.$Chr)]
    ks_test <- ks_test_within_chromosome(., chrom_length)
    data.frame(Chr = unique(.$Chr), D = ks_test$statistic, p.value = ks_test$p.value)
  })

# Print KS test results
ks_results

ks_results <- ks_results %>%
  rename(D_fixed = D, p.value_fixed = p.value)

# Convert to a dataframe
fixed_sds <- as.data.frame(ks_results)

# Print the dataframe
print(fixed_sds)




