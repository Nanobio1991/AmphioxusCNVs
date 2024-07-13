#set working directory 
setwd("//wsl.localhost/Ubuntu/home/nanobio/AmphioxusCNVs/AmphioxusCNVs")

options(scipen = 999)

#livraries
library(readxl)
library(writexl)
library(tidyr)
library(dplyr)
library(ggplot2)

#upload data
cnv_merged <- read_excel("results/CNVpytor/merged_cnvs.xlsx")

#Separate REGION collumn into 3 colllumns for Chr, St, End positions
cnv_merged <- cnv_merged %>%
  separate(col = REGION, into = c("Chr", "st", "end"), sep = "[:-]")

# Convert 'st' and "end" columns to numeric 
cnv_merged$st <- as.numeric(as.character(cnv_merged$st))
cnv_merged$end <- as.numeric(as.character(cnv_merged$end))


# Round up numbers in the dataframe
cnv_merged <- cnv_merged %>%
  mutate(across(where(is.numeric), round))

# Create dataframe for the first population
cnv_pop1 <- cnv_merged[, -c(23:41)]

# Create dataframe for the second population
cnv_pop2 <- cnv_merged[, -c(6:22)]


#Delete rows that contains the same number of copies because they are not CNVs
cnv_merged <- cnv_merged %>%
  filter(apply(.[, -c(1:5)], 1, function(x) length(unique(x)) > 1))

cnv_pop1 <- cnv_pop1 %>%
  filter(apply(.[, -c(1:5)], 1, function(x) length(unique(x)) > 1))

cnv_pop2 <- cnv_pop2 %>%
  filter(apply(.[, -c(1:5)], 1, function(x) length(unique(x)) > 1))

# Add a collumn for gain loss and gain/loss
cnv_merged <- cnv_merged %>%
  mutate(Classification = apply(.[, -c(1:5)], 1, function(x) {
    if(all(x >= 2)) {
      "gain"
    } else if(all(x <= 2)) {
      "loss"
    } else {
      "gain/loss"
    }
  }))

table(cnv_merged$Classification)

# Add a collumn for gain loss and gain/loss
cnv_pop1 <- cnv_pop1 %>%
  mutate(Classification = apply(.[, -c(1:5)], 1, function(x) {
    if(all(x >= 2)) {
      "gain"
    } else if(all(x <= 2)) {
      "loss"
    } else {
      "gain/loss"
    }
  }))

table(cnv_pop1$Classification)

# Add a collumn for gain loss and gain/loss
cnv_pop2 <- cnv_pop2 %>%
  mutate(Classification = apply(.[, -c(1:5)], 1, function(x) {
    if(all(x >= 2)) {
      "gain"
    } else if(all(x <= 2)) {
      "loss"
    } else {
      "gain/loss"
    }
  }))

table(cnv_pop2$Classification)



# CREATE DATA FRAME FOR THE LENGTH OF CHROMOSOMES ################
length_chr <- read.table('data/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt', header = FALSE, sep = "\t")
length_chr <- length_chr %>%
  rename(Chr = V1, Length = V2)
length_chr <- length_chr %>%
  filter(!grepl("^>scaf", Chr))
length_chr$Chr <- sub(">", "", length_chr$Chr)  # Remove ">" from chromosome names



################################################################################################
################################### CHROMOSOME PLOT ############################################
################################################################################################

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
    type <- regDF$Classification[r]
    
    # Assign color based on the type
    color <- ifelse(type == "gain", "blue", ifelse(type == "loss", "orange", "chartreuse3"))
    
    polygon(c(sts, ends, ends, sts), c(p-h/2, p-h/2, p+h/2, p+h/2), col=color, border=NA)
  }
}


pdf("results/plots/Types_od_CNVs_across_chromosomes.pdf",width = 15, height = 8)


plot_all_chr_rows_len(length_chr, cnv_merged, "Types of Copy Number Variations Across chromosomes")
# Add a legend to the plot

legend("bottomright",            # Position of the legend
       legend = c("gain", "loss", "gain/loss"),  # Labels for the legend
       fill = c("blue", "orange", "chartreuse3"), # Colors corresponding to the labels
       title = "CNV Types",    # Title for the legend
       cex = 0.8,             # Character expansion factor for the legend text
       bg = 'white')          # Background color of the legend box


dev.off()



###################################################
#######@@#### PRIVATE CNVS ########################
###################################################


# Load the data
cnv_merged <- read_excel("results/CNVpytor/merged_cnvs.xlsx")

# Round up numbers in the dataframe
cnv_merged <- cnv_merged %>%
  mutate(across(where(is.numeric), round))

# Remove rows where the copy number does not vary across individuals within each row
cnv_merged <- cnv_merged %>%
  filter(apply(.[, -c(1:3)], 1, function(x) length(unique(x)) > 1))




# Split data into two populations
pop1 <- cnv_merged[, c(1:3, 4:20)]
pop2 <- cnv_merged[, c(1:3, 21:39)]

# Find regions in Pop1 where all individuals have the same CNV number
uniform_pop1 <- pop1 %>%
  filter(apply(.[, 4:20], 1, function(x) length(unique(x)) == 1))

# Use the REGION column from uniform_pop1 to filter Pop2
pop2_matched <- pop2 %>%
  filter(REGION %in% uniform_pop1$REGION)


# Find regions in Pop1 where all individuals have the same CNV number
uniform_pop2 <- pop2 %>%
  filter(apply(.[, 4:22], 1, function(x) length(unique(x)) == 1))

# Use the REGION column from uniform_pop1 to filter Pop2
pop1_matched <- pop1 %>%
  filter(REGION %in% uniform_pop2$REGION)



#Separate REGION collumn into 3 colllumns for Chr, St, End positions
pop1_matched <- pop1_matched %>%
  separate(col = REGION, into = c("Chr", "st", "end"), sep = "[:-]")

# Convert 'st' and "end" columns to numeric 
pop1_matched$st <- as.numeric(as.character(pop1_matched$st))
pop1_matched$end <- as.numeric(as.character(pop1_matched$end))


#Separate REGION collumn into 3 colllumns for Chr, St, End positions
pop2_matched <- pop2_matched %>%
  separate(col = REGION, into = c("Chr", "st", "end"), sep = "[:-]")

# Convert 'st' and "end" columns to numeric 
pop2_matched$st <- as.numeric(as.character(pop2_matched$st))
pop2_matched$end <- as.numeric(as.character(pop2_matched$end))

# Add a collumn for gain loss and gain/loss
cnv_merged <- cnv_merged %>%
  mutate(Classification = apply(.[, -c(1:3)], 1, function(x) {
    if(all(x >= 2)) {
      "gain"
    } else if(all(x <= 2)) {
      "loss"
    } else {
      "gain/loss"
    }
  }))


# Add a collumn for gain loss and gain/loss
pop1_matched <- pop1_matched %>%
  mutate(Classification = apply(.[, -c(1:5)], 1, function(x) {
    if(all(x >= 2)) {
      "gain"
    } else if(all(x <= 2)) {
      "loss"
    } else {
      "gain/loss"
    }
  }))

table(pop1_matched$Classification)

# Add a collumn for gain loss and gain/loss
pop2_matched <- pop2_matched %>%
  mutate(Classification = apply(.[, -c(1:5)], 1, function(x) {
    if(all(x >= 2)) {
      "gain"
    } else if(all(x <= 2)) {
      "loss"
    } else {
      "gain/loss"
    }
  }))

table(pop2_matched$Classification)


pdf("results/plots/private_cnvs.pdf",width = 15, height = 8)

par(mfrow=c(1, 2)) 

#Plot
plot_all_chr_rows_len(length_chr, pop1_matched, "Private CNVs Banyuls-sur-Mer")


# Add a legend to the plot
legend("bottomright",            # Position of the legend
       legend = c("gain", "loss", "gain/loss"),  # Labels for the legend
       fill = c("blue", "orange", "chartreuse3"), # Colors corresponding to the labels
       title = "CNV Cassification",    # Title for the legend
       cex = 0.8,             # Character expansion factor for the legend text
       bg = 'white')          # Background color of the legend box



#Plot 
plot_all_chr_rows_len(length_chr, pop2_matched, "Private CNVs Roscoff")


# Add a legend to the plot
legend("bottomright",            # Position of the legend
       legend = c("gain", "loss", "gain/loss"),  # Labels for the legend
       fill = c("blue", "orange", "chartreuse3"), # Colors corresponding to the labels
       title = "CNV Cassification",    # Title for the legend
       cex = 0.8,             # Character expansion factor for the legend text
       bg = 'white')          # Background color of the legend box


par(mfrow=c(1, 1)) 


dev.off()



################################################################################################
############################ LENGHT DISTRIBUTION PLOT ##########################################
################################################################################################

pdf("results/plots/Lenght_distribution_CNVs.pdf", width = 10, height = 8)


# Create the density plot with lines for each classification
ggplot(cnv_merged, aes(x = SIZE, color = Classification)) +
  geom_density() +
  scale_x_continuous(breaks = seq(0, 200000, by = 25000), limits = c(0, 200000)) +
  scale_color_manual(values = c("gain" = "blue", "loss" = "orange", "gain/loss" = "chartreuse3")) +
  labs(title = "CNV Size Distribution", x = "Region length (bp)", y = "Density") +
  theme_minimal() +
  theme(legend.title = element_blank())


dev.off()




