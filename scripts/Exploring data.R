#Load environment 
load("Amphioxus_SDs.RData")

#Libraries 
library(dplyr)
library(ggplot2)

# Load the data
sd_data <- read.table('Branchiostoma_lanceolatum.BraLan3_SDs.bedpe', header = FALSE, sep = "\t")


#Adding a collumn with the mean length of each SDs
sd_data <- sd_data %>%
  mutate(mean_length = ((V3 - V2 + 1) + (V6 - V5 + 1)) / 2)


#Adding collumns for similarity 

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


# Add a column indicating if SDs are on the same chromosome
sd_data <- sd_data %>%
  mutate(same_chromosome = V1 == V4)

###############################################################################
############################## FILTER DATA ####################################
###############################################################################


# Filter rows based on similarity threshold of 90%
filtered_sd_data <- filter(sd_data, similarity >= 90)

#Filter copies that are less than 1000 bp
final_filtered_sd_data <- filtered_sd_data %>%
  filter((V3 - V2) > 999, (V6 - V5) > 999)

#Filter copies that are less than 1000 bp on raw data 
filtered_sd_data2 <- sd_data %>%
  filter((V3 - V2) > 999, (V6 - V5) > 999)


###############################################################################
###########################  HISTOGRAMS OF LENGTHS  ###########################
###############################################################################


#Histogram of SD length of raw data
hist(sd_data$mean_length, breaks = 10000, xlim = c(0, 3000), main = "Raw data SD length")
#Histogram of length of filtered data
hist(filtered_sd_data$mean_length, breaks = 3000, xlim = c(0, 3000), main = "Filtered SD length")


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
ggplot(filtered_sd_data2, aes(x = similarity)) + 
  geom_histogram(color = "black", fill = "blue", bins = 100) +
  labs(title = "Distribution of Percent Identity without SDs that are less than 1000 pb",
       x = "Percent Identity",
       y = "Frequency") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 100, by = 10)) 



# Histogram of Percent Identity with more detailed breaks
ggplot(filtered_sd_data, aes(x = similarity)) + 
  geom_histogram(color = "black", fill = "blue", bins = 100) + # Adjusted for finer detail
  labs(title = "Distribution of Percent Identity without SDs that are less than 90% similar",
       x = "Percent Identity",
       y = "Frequency") +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 100, by = 10))

# Histogram of Percent Identity
ggplot(final_filtered_sd_data, aes(x = similarity)) + 
  geom_histogram(fill = "blue", color = "black", bins = 100) +
  labs(title = "Distribution of Percent Identity without SDs that are less than 90% similar and less than 1000 bp",
       x = "Percent Identity",
       y = "Frequency") +
  theme_minimal()




###############################################################################
#######################  SAME CHROMOSOME HISTOGRAM  ###########################
###############################################################################

#Add a column to see if the duplication is in the same chromosome 
percentage_same_chromosome <- sd_data %>%
  group_by(V1) %>%
  summarise(
    total = n(),
    same_chr_total = sum(same_chromosome, na.rm = TRUE),
    percentage = (same_chr_total / total) * 100
  )



# Keep only the chromosomes !!! Not the scaffolds 
chromosome_similarity <- percentage_same_chromosome %>%
  filter(!grepl("^sNNf", V1))

# Extract numeric parts from the 'V1' column
numeric_part <- as.numeric(sub("Nhr", "", chromosome_similarity$V1))

# Order the dataframe by this numeric sequence
ordered_indices <- order(numeric_part)

# Reorder the dataframe
chromosome_similarity <- chromosome_similarity[ordered_indices, ]

# Convert V1 to a factor with levels in the current order of V1
chromosome_similarity$V1 <- factor(chromosome_similarity$V1, levels = unique(chromosome_similarity$V1))

# Plot
ggplot(chromosome_similarity, aes(x = V1, y = percentage)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Percentage of SDs in the Same Chromosome by Chromosome",
       x = "Chromosome",
       y = "Percentage of SDs in the Same Chromosome") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for readability











