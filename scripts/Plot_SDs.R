setwd("//wsl.localhost/Ubuntu/home/nanobio/AmphioxusCNVs/AmphioxusCNVs/")


#Libraries 
library(dplyr)
library(ggplot2)
library(tidyr)


###############################################################################
##########################  MERGED BAR PLOT  ##################################
#############################################################################

#upload the intra data in R
sd_intra <- read.table('results/filtered_sd_data/pure_intra_cases.bed', header = FALSE, sep = "\t")

sd_intra <- sd_intra %>%
  mutate(V4 = V3 - V2)

#upload the inter data in R
sd_inter <- read.table('results/filtered_sd_data/pure_inter_cases.bed', header = FALSE, sep = "\t")

sd_inter <- sd_inter %>%
  mutate(V4 = V3 - V2)

#upload mixed data in R
sd_mixed <- read.table('results/filtered_sd_data/mixed_sd_cases.bed', header = FALSE, sep = "\t")
sd_mixed <- sd_mixed %>%
  mutate(V4 = V3 - V2)

#################################################
#Prepare the data for plotting 

# Combine all steps into a single pipeline
total_length <- bind_rows(
  sd_intra %>% mutate(Type = "Intra"),
  sd_inter %>% mutate(Type = "Inter"),
  sd_mixed %>% mutate(Type = "Intersection")
) %>%
  group_by(V1, Type) %>%
  summarise(Total_Length = sum(V4, na.rm = TRUE), .groups = "drop") %>%
  ungroup() %>%
  mutate(V1 = factor(V1, levels = unique(V1))) %>%
  arrange(V1, Type)



# Put it in numeric order
numeric_part <- as.numeric(sub("chr", "", total_length$V1))
ordered_indices <- order(numeric_part)
total_length <- total_length[ordered_indices, ]
total_length$V1 <- factor(total_length$V1, levels = unique(total_length$V1))


pdf("results/plots/Length_of_SD_type_per_chromosome.pdf", width = 10, height = 8)

# Plot
ggplot(total_length, aes(x = V1, y = Total_Length, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Intra" = "blue", "Intersection" = "green", "Inter" = "orange")) +
  labs(title = "Length of SD type per Chromosome",
       x = "Chromosome",
       y = "Total Length",
       fill = "SD Category") +
  theme_minimal()

dev.off()



###############################################################################
############ PLOT GRAPH PER CHROMOSMES ########################################
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
sd_positions_merged <- read.table('results/filtered_sd_data/sd_positions_merged.bed', header = FALSE, sep = "\t")
##############################################




sd_positions_merged <- sd_positions_merged %>%
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


pdf("results/plots/sd_chr_plot1.pdf.pdf")

plot_all_chr_rows_len(length_chr, sd_positions_merged, "blue")



dev.off()







#################################################################################
############## CORRECTED CHROMOSOME PLOT ########################################
#################################################################################


# Combine into one dataframe and add a 'Type' column
combined_sd_df <- bind_rows(
  sd_intra %>% mutate(Type = "Intra"),
  sd_inter %>% mutate(Type = "Inter"),
  sd_mixed %>% mutate(Type = "Mixed")
)

# Now you can use your combined_sd_df with the function

combined_sd_df <- combined_sd_df %>%
  rename(Chr = V1, st = V2, end = V3)



#Create function and plot 
plot_all_chr_rows_len <- function(lenDF, regDF, plot_title = "Chromosome Visualization"){
  height <- 0.6
  maxY <- length(lenDF$Chr) * (height + 0.2) # Adjust spacing
  maxChrLength <- max(lenDF$Length) # Find the maximum chromosome length
  
  # Adjust plot to reflect chromosome lengths properly
  plot(c(1:10), c(1:10), axes=F, xlab="", ylab="", ylim=c(0, maxY), xlim=c(0, maxChrLength), col=NA)
  
  for(c in 1:length(lenDF$Chr)){
    chrName <- lenDF$Chr[c]
    pos <- maxY - (c * (height + 0.2)) # Adjust for spacing
    chromosomeLength <- lenDF$Length[c]
    
    # Filter regDF for the current chromosome
    currentRegDF <- regDF[regDF$Chr == chrName,]
    
    # Draw the chromosome and the SDs
    draw_chr_row(height, pos, chromosomeLength, currentRegDF)
    polygon(c(0, chromosomeLength, chromosomeLength, 0), c(pos-height/2, pos-height/2, pos+height/2, pos+height/2), col=NA, border="black")
  }
  
  # Create descriptive labels for each chromosome
  chrLabels <- lenDF$Chr
  
  # Add y-axis with descriptive labels
  axis(2, at = seq(0, maxY - (height + 0.2), by = (height + 0.2)), labels=rev(chrLabels), las=1, cex.axis=0.7)
  
  # Adding a title to the plot
  title(main = plot_title, col.main = "black", font.main = 1)
}


draw_chr_row <- function(h, p, chromosomeLength, regDF){
  for(r in 1:nrow(regDF)){
    sts <- regDF$st[r] / chromosomeLength * chromosomeLength
    ends <- regDF$end[r] / chromosomeLength * chromosomeLength
    type <- regDF$Type[r]
    
    # Assign color based on the type
    color <- ifelse(type == "Intra", "blue", ifelse(type == "Inter", "orange", "chartreuse3"))
    
    polygon(c(sts, ends, ends, sts), c(p-h/2, p-h/2, p+h/2, p+h/2), col=color, border=NA)
  }
}


pdf("results/plots/Types_of_SDs_across_chromosomes.pdf",width = 15, height = 8)


plot_all_chr_rows_len(length_chr, combined_sd_df, "Types of Segmental Duplications across chromosomes")
# Add a legend to the plot

legend("bottomright",            # Position of the legend
       legend = c("Intra", "Inter", "Intersection"),  # Labels for the legend
       fill = c("blue", "orange", "chartreuse3"), # Colors corresponding to the labels
       title = "SD Types",    # Title for the legend
       cex = 0.8,             # Character expansion factor for the legend text
       bg = 'white')          # Background color of the legend box


dev.off()



######################################################################################


######################################################################################
# COMPLETE THE GRAPH WITH THE REPEAT REGIONS 
######################################################################
#ON BASH : 

# EXTRACT POSITION FILES FOR TRF MASKED REGIONS

#grep -v "^$" Branchiostoma_lanceolatum.BraLan3_genome.fa.masked.2.7.7.80.10.50.15.dat | awk '/Sequence:/ {chromosome=$2} {print chromosome, $1, $2}' > repeats_info.txt
#tail -n +6 repeats_info.txt | grep -v "Sequence" | grep -v "Parameters" |grep -v "scaf" > repeats_trf.txt

repeats_trf <- read.table('results/mask_tandem_repeats_with_trf/repeats_trf.txt', header = FALSE, sep = " ")

#EXTRACT POSITIONS FOR THE REPEATMASKER MASKED REGIONS

#awk '{if (NR > 1) print $5, $6, $7}' Branchiostoma_lanceolatum.BraLan3_genome.fa.masked.2.7.7.80.10.50.15.mask.out | tail -n +3 | grep -v "scaf" > repeats_repeatmasker.txt

#when uploading data that doesn't have collumns just put sep " " 
repeats_repeatmasker <- read.table('results/mask_genome_with_repeatmasker/repeats_repeatmasker.txt', header = FALSE, sep = " ")




# Combine into one dataframe and add a 'Type' column
combined_sd_df2 <- bind_rows(
  sd_intra %>% mutate(Type = "Intra"),
  sd_inter %>% mutate(Type = "Inter"),
  sd_mixed %>% mutate(Type = "Mixed"),
  repeats_repeatmasker %>% mutate(Type = "repeatmasker"),
  repeats_trf %>% mutate(Type = "trf")
)

# use combined_sd_df with the function
combined_sd_df2 <- combined_sd_df2 %>%
  rename(Chr = V1, st = V2, end = V3)


draw_chr_row <- function(h, p, chromosomeLength, regDF){
  for(r in 1:nrow(regDF)){
    sts <- regDF$st[r] / chromosomeLength * chromosomeLength
    ends <- regDF$end[r] / chromosomeLength * chromosomeLength
    type <- regDF$Type[r]
    
    # Assign color based on the type
    color <- ifelse(type == "Intra", "blue", 
                    ifelse(type == "Inter", "blue", 
                           ifelse(type == "Mixed", "blue",
                                  ifelse(type == "trf", "red", "red"))))
    
    polygon(c(sts, ends, ends, sts), c(p-h/2, p-h/2, p+h/2, p+h/2), col=color, border=NA)
  }
}

pdf("results/plots/SD_and_repeats_across_chromosomes.pdf",width = 15, height = 8)

plot_all_chr_rows_len(length_chr, combined_sd_df2, "Segmental duplications and repeats across chromosomes")
# Add a legend to the plot

legend("bottomright",            # Position of the legend
       legend = c("SDs", "Repeats"),  # Labels for the legend
       fill = c("blue", "red"), # Colors corresponding to the labels
       title = "SD Types",    # Title for the legend
       cex = 0.8,             # Character expansion factor for the legend text
       bg = 'white')          # Background color of the legend box

dev.off()




#################################################################################
############## LENGTH DISTRIBUTION PLOT #########################################
#################################################################################




# Add a category column to each dataframe
sd_intra <- sd_intra %>% mutate(Category = 'Intra')
sd_inter <- sd_inter %>% mutate(Category = 'Inter')
sd_mixed <- sd_mixed %>% mutate(Category = 'Intersection')

# Combine the dataframes into one
combined_sd <- bind_rows(sd_intra, sd_inter, sd_mixed)



pdf("results/plots/length_distribution_SDs.pdf",width = 15, height = 8)

# Create the density plot with lines for each classification
ggplot(combined_sd, aes(x = V4, color = Category)) +
  geom_density(size = 0.7) +  # Increase line thickness
  scale_x_continuous(limits = c(0, 7500), breaks = seq(0, 20000, by = 1000)) +  # Set x-axis limits and breaks
  scale_color_manual(values = c('Intra' = 'blue', 'Inter' = 'orange', 'Intersection' = 'green')) +
  labs(title = "Length Distribution of Segmental Duplications", x = "Length of SDs (bp)", y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())

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




