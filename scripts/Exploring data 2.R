# Set Working directory 
setwd("//wsl.localhost/Ubuntu/home/nanobio/finding_SDs_3")

#Load environment 
load("Amphioxus_SDs.RData")


#save.image("Amphioxus_SDs.RData")


#Libraries 
library(dplyr)
library(ggplot2)
library(tidyr)


# Load the data
sd_data <- read.table('Branchiostoma_lanceolatum.BraLan3_SDs.bedpe', header = FALSE, sep = "\t")


###############################################################################
############################ ADDING COLLUMNS ##################################
###############################################################################

###Adding a collumn with the mean length of each SDs
sd_data <- sd_data %>%
  mutate(mean_length = ((V3 - V2 + 1) + (V6 - V5 + 1)) / 2)


###Adding collumns for similarity 

# Function to calculate total alignment span and total number of matches
calculate_totals <- function(cigar) {
  # Extract numbers and letters (operations) from the CIGAR string
  numbers <- as.numeric(unlist(regmatches(cigar, gregexpr("\\d+", cigar))))
  operations <- unlist(regmatches(cigar, gregexpr("[A-Z]", cigar)))
  
  # Total alignment span is the sum of all numbers
  total_span <- sum(numbers)
  
  # Total number of matches is the sum of numbers before an 'M'
  match_numbers <- numbers[operations == "M"]
  total_matches <- sum(match_numbers)
  
  # Return both totals
  list(total_span = total_span, total_matches = total_matches)
}
# Apply the function to column V13 and calculate similarity
sd_data <- sd_data %>%
  rowwise() %>%
  mutate(
    totals = list(calculate_totals(V13)),
    total_span = totals$total_span,
    total_matches = totals$total_matches,
    similarity = (total_matches / total_span) * 100
  ) %>%
  ungroup()


### Add a column indicating if SDs are on the same chromosome
sd_data <- sd_data %>%
  mutate(same_chromosome = V1 == V4)


###############################################################################
############################## FILTER DATA ####################################
###############################################################################

# Filter out rows where either V1 or V4 contain scaffold identifiers
sd_data_filtered <- sd_data %>%
  filter(!grepl("^sNNf", V1) & !grepl("^sNNf", V4))

# Filter rows based on similarity threshold of 90%
sd_data_filtered <- filter(sd_data_filtered, similarity >= 90)

#Filter copies that are less than 1000 bp
sd_data_filtered <- sd_data_filtered %>%
  filter((V3 - V2) > 999, (V6 - V5) > 999)

sd_data_filtered <- sd_data_filtered %>%
  mutate(
    V1 = ifelse(grepl("^Nhr", V1), sub("^Nhr", "chr", V1), V1),
    V4 = ifelse(grepl("^Nhr", V4), sub("^Nhr", "chr", V4), V4)
  )





###############################################################################
###########################  HISTOGRAMS OF LENGTHS  ###########################
###############################################################################


#Histogram of SD length of raw data
hist(sd_data$mean_length, breaks = 10000, xlim = c(0, 3000), main = "Raw data SD length")
#Histogram of length of filtered data
hist(sd_data_filtered$mean_length, breaks = 3000, xlim = c(0, 3000), main = "Filtered SD length")


# Create a histogram with two colors based on the 'same_chromosome' column
hist(sd_data_filtered$mean_length[sd_data_filtered$same_chromosome == TRUE], 
     breaks = 5000, 
     xlim = c(0, 3000),
     ylim = c(0, 80000),
     col = "blue", 
     main = "Filtered SD length - Same Chromosome",
     xlab = "bp")


# Add histogram for 'same_chromosome == FALSE' with a different color
hist(sd_data_filtered$mean_length[sd_data_filtered$same_chromosome == FALSE], 
     breaks = 3000, 
     xlim = c(0, 3000),
     ylim = c(0, 75000),
     col = rgb(1, 0.5, 0, 0.5), 
     add = TRUE)

# Add legend
legend("topright", legend = c("Same Chromosome", "Different Chromosome"), fill = c("blue", "orange"))




###############################################################################
########################  HISTOGRAMS OF SIMILARITY  ###########################
###############################################################################


# Histogram of Percent Identity of raw data 
ggplot(sd_data, aes(x = similarity)) + 
  geom_histogram(color = "black", fill = "blue", bins = 100) +
  labs(title = "Distribution of Percent Identity of raw data",
       x = "Percent Identity",
       y = "Frequency") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) 

# Histogram of Percent Identity without SDs that are less than 1000 pb
ggplot(sd_data_filtered, aes(x = similarity)) + 
  geom_histogram(color = "black", fill = "blue", bins = 100) +
  labs(title = "Distribution of Percent Identity without SDs that are less than 1000 pb",
       x = "Percent Identity",
       y = "Frequency") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) 




###############################################################################
#######################  SAME CHROMOSOME HISTOGRAM  ###########################
###############################################################################

calculate_sd_counts <- function(sd_data_filtered) {
  # Calculate counts for the first copy
  counts_v1 <- sd_data_filtered %>%
    group_by(V1) %>%
    summarize(
      intra_chromosomal = sum(same_chromosome == TRUE),
      inter_chromosomal = sum(same_chromosome == FALSE)
    ) %>%
    rename(Chr = V1)
  
  # Calculate counts for the second copy
  counts_v4 <- sd_data_filtered %>%
    group_by(V4) %>%
    summarize(
      intra_chromosomal = sum(same_chromosome == TRUE), # This will be 0 for V4 if same_chromosome is TRUE
      inter_chromosomal = sum(same_chromosome == FALSE)
    ) %>%
    rename(Chr = V4)
  
  # Combine counts from V1 and V4, summing inter-chromosomal counts
  combined_counts <- bind_rows(counts_v1, counts_v4) %>%
    group_by(Chr) %>%
    summarize(
      intra_chromosomal = sum(intra_chromosomal),
      inter_chromosomal = sum(inter_chromosomal),
      .groups = 'drop' # Drop grouping
    ) %>%
    # Calculate total SDs per chromosome
    mutate(total_SDs = intra_chromosomal + inter_chromosomal,
           # Calculate percentage of intra-chromosomal SDs
           percent_intra = (intra_chromosomal / total_SDs) * 100,
           # Calculate percentage of inter-chromosomal SDs
           percent_inter = (inter_chromosomal / total_SDs) * 100)
  
  return(combined_counts)
}

# Use the function with your data frame
sd_counts_per_chr <- calculate_sd_counts(sd_data_filtered)
print(sd_counts_per_chr)


# Extract numeric parts from the 'V1' column
numeric_part <- as.numeric(sub("chr", "", sd_counts_per_chr$Chr))

# Order the dataframe by this numeric sequence
ordered_indices <- order(numeric_part)

# Reorder the dataframe
sd_counts_per_chr <- sd_counts_per_chr[ordered_indices, ]

# Convert V1 to a factor with levels in the current order of V1
sd_counts_per_chr$Chr <- factor(sd_counts_per_chr$Chr, levels = unique(sd_counts_per_chr$Chr))


# Plot
ggplot(sd_counts_per_chr, aes(x = Chr, y = percent_intra)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Percentage of SDs in the Same Chromosome",
       x = "Chromosome",
       y = "Percentage of SDs in the Same Chromosome") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for readability




# Transform the data to long format for ggplot
sd_counts_long <- tidyr::pivot_longer(sd_counts_per_chr,
                                      cols = c("percent_intra", "percent_inter"),
                                      names_to = "SD_Type",
                                      values_to = "Percentage")

# Adjust SD_Type for clear labeling
sd_counts_long$SD_Type <- recode(sd_counts_long$SD_Type,
                                 percent_intra = "Intra-chromosomal",
                                 percent_inter = "Inter-chromosomal")

# Plot with stacked bars
ggplot(sd_counts_long, aes(x = Chr, y = Percentage, fill = SD_Type)) +
  geom_bar(stat = "identity") + # Default position is stack
  scale_fill_manual(values = c("Intra-chromosomal" = "blue", "Inter-chromosomal" = "orange")) +
  labs(title = "Proportion of Intra- vs. Inter-chromosomal SDs per Chromosome",
       x = "Chromosome",
       y = "Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



############################################################################################
############ NOW WITH CORRECTED DATA ?????????????? ########################################
############################################################################################



calculate_sd_counts <- function(sd_data_filtered) {
  # Initial counts for each chromosome based on the first copy (V1)
  counts_v1 <- sd_data_filtered %>%
    group_by(V1) %>%
    summarize(
      intra_chromosomal = sum(same_chromosome == TRUE),
      inter_chromosomal = sum(same_chromosome == FALSE),
      .groups = 'drop'
    ) %>%
    rename(Chr = V1)
  
  # Initial counts for each chromosome based on the second copy (V4), noting that
  # intra_chromosomal should be 0 for these since it's inter-chromosomal if V1 != V4
  counts_v4 <- sd_data_filtered %>%
    group_by(V4) %>%
    summarize(
      inter_chromosomal = sum(same_chromosome == FALSE),
      .groups = 'drop'
    ) %>%
    mutate(intra_chromosomal = 0) %>% # Ensure intra_chromosomal is 0 for V4
    rename(Chr = V4)
  
  # Combine counts from V1 and V4, summing inter-chromosomal counts correctly
  combined_counts <- bind_rows(counts_v1, counts_v4) %>%
    group_by(Chr) %>%
    summarize(
      intra_chromosomal = sum(intra_chromosomal),
      inter_chromosomal = sum(inter_chromosomal),
      .groups = 'drop'
    ) %>%
    # Calculate total SDs per chromosome
    mutate(total_SDs = intra_chromosomal + inter_chromosomal,
           # Calculate percentage of intra-chromosomal SDs
           percent_intra = (intra_chromosomal / total_SDs) * 100,
           # Calculate percentage of inter-chromosomal SDs
           percent_inter = (inter_chromosomal / total_SDs) * 100)
  
  return(combined_counts)
}

# Use the function with your data frame
sd_counts_per_chr2 <- calculate_sd_counts(sd_data_filtered)
print(sd_counts_per_chr2)

# Extract numeric parts from the 'V1' column
numeric_part <- as.numeric(sub("chr", "", sd_counts_per_chr2$Chr))

# Order the dataframe by this numeric sequence
ordered_indices <- order(numeric_part)

# Reorder the dataframe
sd_counts_per_chr2 <- sd_counts_per_chr2[ordered_indices, ]

# Convert V1 to a factor with levels in the current order of V1
sd_counts_per_chr2$Chr <- factor(sd_counts_per_chr2$Chr, levels = unique(sd_counts_per_chr2$Chr))



# Plot
ggplot(sd_counts_per_chr2, aes(x = Chr, y = percent_intra)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Percentage of SDs in the Same Chromosome",
       x = "Chromosome",
       y = "Percentage of SDs in the Same Chromosome") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for readability


#########################################################################################################


# the length of the totaol SDs of chromosome 1 by the total length of chromosome 1 


sd_positions_merged$length_copy <- sd_positions_merged$end - sd_positions_merged$st


# Assuming sd_positions_merged is your dataframe
# Create the length_copy column
sd_positions_merged$length_copy <- sd_positions_merged$end - sd_positions_merged$st

# Sum the length_copy for each chromosome
total_length_copy_per_chr <- aggregate(length_copy ~ Chr, data = sd_positions_merged, FUN = sum)

# Now you can merge this back into your length_chr dataframe if you have one
# Assuming length_chr is your dataframe with chromosome lengths
length_chr <- merge(length_chr, total_length_copy_per_chr, by = 'Chr', all.x = TRUE)

# Replace NAs with 0 if there are chromosomes without any SDs
length_chr$length_copy[is.na(length_chr$length_copy)] <- 0



# Calculate the percentage
length_chr$percentage_length_copy <- (length_chr$length_copy / length_chr$Length) * 100

# Replace any potential NaN or Inf values with 0 (which may occur if Length is 0)
length_chr$percentage_length_copy[is.na(length_chr$percentage_length_copy) | is.infinite(length_chr$percentage_length_copy)] <- 0

#Plot

ggplot(length_chr, aes(x = Chr, y = percentage_length_copy)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Percentage of SDs length in the the chromosome ",
       x = "Chromosome",
       y = "Percentage of SDs length") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for readability


##############################################################################################

#Histogram of SD length of raw data
hist(same_chr_sd$mean_length, breaks = 5000, xlim = c(0, 3000), main = "Same chr Length")
mean(same_chr_sd$mean_length)


#Histogram of length of filtered data
hist(different_chr_sd$mean_length, breaks = 3000, xlim = c(0, 3000), main = "Diff chr Length")
mean(different_chr_sd$mean_length)



###############################################################################
############ PLOT GRAPH PER CHROMOSMES ########################################
###############################################################################

## CREATE DATA FRAME FOR THE LENGTH OF CHROMOSOMES ################

length_chr <- read.table('Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt', header = FALSE, sep = "\t")
head(length_chr)

# Ensure lenDF is prepared correctly
length_chr <- length_chr %>%
  rename(Chr = V1, Length = V2)

length_chr <- length_chr %>%
  filter(!grepl("^>scaf", Chr))

length_chr$Chr <- sub(">", "", length_chr$Chr)  # Remove ">" from chromosome names


############### CREATE DATA FRAME OF EACH COPY #######################

# Create a dataframe for the first copy
first_copy_df <- sd_data_filtered %>%
  select(Chr = V1, Start = V2, End = V3)

# Create a dataframe for the second copy
second_copy_df <- sd_data_filtered %>%
  select(Chr = V4, Start = V5, End = V6)

# Combine the two dataframes
sd_positions <- bind_rows(first_copy_df, second_copy_df)
head(sd_positions)

# Adjusting sd_positions column names
sd_positions <- sd_positions %>%
  rename(st = Start, end = End)

# If your chromosome names do not start with "chr", adjust the sub function accordingly.
numeric_part <- as.numeric(sub("chr", "", sd_positions$Chr))

# Create a new dataframe with an additional column for the numeric part
sd_positions <- sd_positions %>%
  mutate(Chr_numeric = numeric_part)

# Sort by Chr_numeric, then by Start
sd_positions_sorted <- sd_positions %>%
  arrange(Chr_numeric, st) %>%
  select(-Chr_numeric)  # Remove the Chr_numeric column after sorting

# Export the sorted dataframe to a BED file
write.table(sd_positions_sorted, "sd_positions_sorted.bed",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


##################### MERGE WITH BEDTOOLS ############################################

sd_positions_merged <- read.table('sd_positions_merged.bed', header = FALSE, sep = "\t")

############################################################


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

plot_all_chr_rows_len(length_chr, sd_positions_merged, "blue")






###############################################################################
##########################  MERGED BAR PLOT  ##################################
#############################################################################
# FILTER DATA BY INTRA AND INTER

# Filter rows where same_chromosome is TRUE
sd_data_intra <- filter(sd_data_filtered, same_chromosome == TRUE)

# Create a dataframe for the first copy
first_copy_intra <- sd_data_intra %>%
  select(Chr = V1, Start = V2, End = V3)

# Create a dataframe for the second copy
second_copy_intra <- sd_data_intra %>%
  select(Chr = V4, Start = V5, End = V6)

# Combine the two dataframes
sd_positions_intra <- bind_rows(first_copy_intra, second_copy_intra)
head(sd_positions_intra)

# Adjusting sd_positions column names
sd_positions_intra <- sd_positions_intra %>%
  rename(st = Start, end = End)

# If your chromosome names do not start with "chr", adjust the sub function accordingly.
numeric_part <- as.numeric(sub("chr", "", sd_positions_intra$Chr))

# Create a new dataframe with an additional column for the numeric part
sd_positions_intra <- sd_positions_intra %>%
  mutate(Chr_numeric = numeric_part)

# Sort by Chr_numeric, then by Start
sd_positions_intra <- sd_positions_intra %>%
  arrange(Chr_numeric, st) %>%
  select(-Chr_numeric)  # Remove the Chr_numeric column after sorting

# Export the sorted dataframe to a BED file
write.table(sd_positions_intra, "sd_positions_intra.bed",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)



#####################################################
# NOW INTER


# Filter rows where same_chromosome is FALSE
sd_data_inter <- filter(sd_data_filtered, same_chromosome == FALSE)

# Create a dataframe for the first copy
first_copy_inter <- sd_data_inter %>%
  select(Chr = V1, Start = V2, End = V3)

# Create a dataframe for the second copy
second_copy_inter <- sd_data_inter %>%
  select(Chr = V4, Start = V5, End = V6)

# Combine the two dataframes
sd_positions_inter <- bind_rows(first_copy_inter, second_copy_inter)
head(sd_positions_inter)

# Adjusting sd_positions column names
sd_positions_inter <- sd_positions_inter %>%
  rename(st = Start, end = End)

# If your chromosome names do not start with "chr", adjust the sub function accordingly.
numeric_part <- as.numeric(sub("chr", "", sd_positions_inter$Chr))

# Create a new dataframe with an additional column for the numeric part
sd_positions_inter <- sd_positions_inter %>%
  mutate(Chr_numeric = numeric_part)

# Sort by Chr_numeric, then by Start
sd_positions_inter <- sd_positions_inter %>%
  arrange(Chr_numeric, st) %>%
  select(-Chr_numeric)  # Remove the Chr_numeric column after sorting

# Export the sorted dataframe to a BED file
write.table(sd_positions_inter, "sd_positions_inter.bed",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


################ MERGE WITH BEDTOOLS ##################
#Install bedtools in your conda env 

#merge sds that overlap 
#bedtools merge -i sd_positions_intra.bed > sd_positions_intra_merged.bed
#bedtools merge -i sd_positions_inter.bed > sd_positions_inter_merged.bed

#find SDs that are present in both your intrachromosomal and interchromosomal merged datasets. These are your mixed cases
#bedtools intersect -a sd_positions_intra_merged.bed -b sd_positions_inter_merged.bed > mixed_sd_cases.bed

#find SDs in the intrachromosomal dataset that do not overlap with the interchromosomal dataset.
#bedtools intersect -v -a sd_positions_intra_merged.bed -b sd_positions_inter_merged.bed > pure_intra_cases.bed

#find SDs in the interchromosomal dataset that do not overlap with the intrachromosomal dataset.
#bedtools intersect -v -a sd_positions_inter_merged.bed -b sd_positions_intra_merged.bed > pure_inter_cases.bed

#########################################################################

#upload the intra data in R
sd_intra <- read.table('pure_intra_cases.bed', header = FALSE, sep = "\t")

sd_intra <- sd_intra %>%
  mutate(V4 = V3 - V2)

#upload the inter data in R
sd_inter <- read.table('pure_inter_cases.bed', header = FALSE, sep = "\t")

sd_inter <- sd_inter %>%
  mutate(V4 = V3 - V2)

#upload mixed data in R
sd_mixed <- read.table('mixed_sd_cases.bed', header = FALSE, sep = "\t")
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



# Plot
ggplot(total_length, aes(x = V1, y = Total_Length, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Intra" = "blue", "Intersection" = "chartreuse3", "Inter" = "orange")) +
  labs(title = "Total Length of SDs by Chromosome",
       x = "Chromosome",
       y = "Total Length",
       fill = "SD Category") +
  theme_minimal()









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
plot_all_chr_rows_len <- function(lenDF, regDF){
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

plot_all_chr_rows_len(length_chr, combined_sd_df)
# Add a legend to the plot

legend("bottomright",            # Position of the legend
       legend = c("Intra", "Inter", "Intersection"),  # Labels for the legend
       fill = c("blue", "orange", "chartreuse3"), # Colors corresponding to the labels
       title = "SD Types",    # Title for the legend
       cex = 0.8,             # Character expansion factor for the legend text
       bg = 'white')          # Background color of the legend box




######################################################################################


######################################################################################
# COMPLETE THE GRAPH WITH THE REPEAT REGIONS 
######################################################################
#ON BASH : 

# EXTRACT POSITION FILES FOR TRF MASKED REGIONS

#grep -v "^$" Branchiostoma_lanceolatum.BraLan3_genome.fa.masked.2.7.7.80.10.50.15.dat | awk '/Sequence:/ {chromosome=$2} {print chromosome, $1, $2}' > repeats_info.txt
#tail -n +6 repeats_info.txt | grep -v "Sequence" | grep -v "Parameters" |grep -v "scaf" > repeats_trf.txt

repeats_trf <- read.table('/home/nanobio/AmphioxusCNVs/AmphioxusCNVs/repeats_trf.txt', header = FALSE, sep = " ")

#EXTRACT POSITIONS FOR THE REPEATMASKER MASKED REGIONS

#awk '{if (NR > 1) print $5, $6, $7}' Branchiostoma_lanceolatum.BraLan3_genome.fa.masked.2.7.7.80.10.50.15.mask.out | tail -n +3 | grep -v "scaf" > repeats_repeatmasker.txt

#when uploading data that doesn't have collumns just put sep " " 
repeats_repeatmasker <- read.table('/home/nanobio/AmphioxusCNVs/AmphioxusCNVs/results/mask_genome_with_repeatmasker2/repeats_repeatmasker.txt', header = FALSE, sep = " ")




#################################################################################
############## CORRECTED CHROMOSOME PLOT 2 ########################################
#################################################################################


# Combine into one dataframe and add a 'Type' column
combined_sd_df2 <- bind_rows(
  sd_intra %>% mutate(Type = "Intra"),
  sd_inter %>% mutate(Type = "Inter"),
  sd_mixed %>% mutate(Type = "Mixed"),
  repeats_repeatmasker %>% mutate(Type = "repeatmasker"),
  repeats_trf %>% mutate(Type = "trf")
)

# Now you can use your combined_sd_df with the function

combined_sd_df2 <- combined_sd_df2 %>%
  rename(Chr = V1, st = V2, end = V3)


#Create function and plot 
plot_all_chr_rows_len <- function(lenDF, regDF){
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
  
}

draw_chr_row <- function(h, p, chromosomeLength, regDF){
  for(r in 1:nrow(regDF)){
    sts <- regDF$st[r] / chromosomeLength * chromosomeLength
    ends <- regDF$end[r] / chromosomeLength * chromosomeLength
    type <- regDF$Type[r]
    
    # Assign color based on the type
    color <- ifelse(type == "Intra", "blue", 
                    ifelse(type == "Inter", "orange", 
                           ifelse(type == "Mixed", "chartreuse3",
                                  ifelse(type == "trf", "red", "red"))))
    
    polygon(c(sts, ends, ends, sts), c(p-h/2, p-h/2, p+h/2, p+h/2), col=color, border=NA)
  }
}


plot_all_chr_rows_len(length_chr, combined_sd_df2)
# Add a legend to the plot

legend("bottomright",            # Position of the legend
       legend = c("Intra", "Inter", "Intersection", "Repeats"),  # Labels for the legend
       fill = c("blue", "orange", "chartreuse3", "red"), # Colors corresponding to the labels
       title = "SD Types",    # Title for the legend
       cex = 0.8,             # Character expansion factor for the legend text
       bg = 'white')          # Background color of the legend box








