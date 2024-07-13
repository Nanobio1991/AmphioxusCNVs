
#Libraries 
library(dplyr)
library(ggplot2)
library(tidyr)


# Load the data
sd_data <- read.table('results/finding_SDs/Branchiostoma_lanceolatum.BraLan3_SDs.bedpe', header = FALSE, sep = "\t")


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

sd_data_filtered <- sd_data_filtered %>%
  select(-V7, -V8, -V9, -V10, -V11, -V12, -V13, -V14)



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

# Adjusting sd_positions column names
sd_positions_intra <- sd_positions_intra %>%
  rename(st = Start, end = End)

numeric_part <- as.numeric(sub("chr", "", sd_positions_intra$Chr))

# Create a new dataframe with an additional column for the numeric part
sd_positions_intra <- sd_positions_intra %>%
  mutate(Chr_numeric = numeric_part)

# Sort by Chr_numeric, then by Start
sd_positions_intra <- sd_positions_intra %>%
  arrange(Chr_numeric, st) %>%
  select(-Chr_numeric)  # Remove the Chr_numeric column after sorting

# Export the sorted dataframe to a BED file
write.table(sd_positions_intra, "results/filtered_sd_data/sd_positions_intra.bed",
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

# Adjusting sd_positions column names
sd_positions_inter <- sd_positions_inter %>%
  rename(st = Start, end = End)

numeric_part <- as.numeric(sub("chr", "", sd_positions_inter$Chr))

# Create a new dataframe with an additional column for the numeric part
sd_positions_inter <- sd_positions_inter %>%
  mutate(Chr_numeric = numeric_part)

# Sort by Chr_numeric, then by Start
sd_positions_inter <- sd_positions_inter %>%
  arrange(Chr_numeric, st) %>%
  select(-Chr_numeric)  # Remove the Chr_numeric column after sorting

# Export the sorted dataframe to a BED file
write.table(sd_positions_inter, "results/filtered_sd_data/sd_positions_inter.bed",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


################ MERGE WITH BEDTOOLS ##################
#Install bedtools in conda env 

#merge sds that overlap 
#bedtools merge -i sd_positions_intra.bed > sd_positions_intra_merged.bed
#bedtools merge -i sd_positions_inter.bed > sd_positions_inter_merged.bed

#find SDs that are present in both intrachromosomal and interchromosomal merged datasets. These are mixed cases
#bedtools intersect -a sd_positions_intra_merged.bed -b sd_positions_inter_merged.bed > mixed_sd_cases.bed

#find SDs in the intrachromosomal dataset that do not overlap with the interchromosomal dataset.
#bedtools intersect -v -a sd_positions_intra_merged.bed -b sd_positions_inter_merged.bed > pure_intra_cases.bed

#find SDs in the interchromosomal dataset that do not overlap with the intrachromosomal dataset.
#bedtools intersect -v -a sd_positions_inter_merged.bed -b sd_positions_intra_merged.bed > pure_inter_cases.bed

#########################################################################




###############################################################################
############ PLOT GRAPH PER CHROMOSMES ########################################
###############################################################################


############### CREATE DATA FRAME OF EACH COPY #######################

# Create a dataframe for the first copy
first_copy_df <- sd_data_filtered %>%
  select(Chr = V1, Start = V2, End = V3)

# Create a dataframe for the second copy
second_copy_df <- sd_data_filtered %>%
  select(Chr = V4, Start = V5, End = V6)

# Combine the two dataframes
sd_positions <- bind_rows(first_copy_df, second_copy_df)

# Adjusting sd_positions column names
sd_positions <- sd_positions %>%
  rename(st = Start, end = End)

numeric_part <- as.numeric(sub("chr", "", sd_positions$Chr))

# Create a new dataframe with an additional column for the numeric part
sd_positions <- sd_positions %>%
  mutate(Chr_numeric = numeric_part)

# Sort by Chr_numeric, then by Start
sd_positions_sorted <- sd_positions %>%
  arrange(Chr_numeric, st) %>%
  select(-Chr_numeric)  # Remove the Chr_numeric column after sorting

# Export the sorted dataframe to a BED file
write.table(sd_positions_sorted, "results/filtered_sd_data/sd_positions_sorted.bed",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

##################### MERGE WITH BEDTOOLS ############################################
#bedtools merge -i sd_positions_sorted.bed > sd_positions_merged.bed


