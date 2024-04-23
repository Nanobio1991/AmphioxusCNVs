#set working directory 
setwd("//wsl.localhost/Ubuntu/home/nanobio/AmphioxusCNVs/AmphioxusCNVs")

#livraries
library(readxl)
library(tidyr)
library(dplyr)

#upload data
cnv_merged <- read_excel("results/CNVpytor/merged_cnvs.xlsx")

#Separate REGION collumn into 3 colllumns for Chr, St, End positions
cnv_merged <- cnv_merged %>%
  separate(col = REGION, into = c("Chr", "st", "end"), sep = "[:-]")


# Round up numbers in the dataframe
cnv_merged <- cnv_merged %>%
  mutate(across(where(is.numeric), round))


# Assuming cnv_merged is your dataframe
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


# Convert 'st' and "end" columns to numeric 
cnv_merged$st <- as.numeric(as.character(cnv_merged$st))
cnv_merged$end <- as.numeric(as.character(cnv_merged$end))





## CREATE DATA FRAME FOR THE LENGTH OF CHROMOSOMES ################
length_chr <- read.table('data/Branchiostoma_lanceolatum.BraLan3_chr_lengths.txt', header = FALSE, sep = "\t")
head(length_chr)
length_chr <- length_chr %>%
  rename(Chr = V1, Length = V2)
length_chr <- length_chr %>%
  filter(!grepl("^>scaf", Chr))
length_chr$Chr <- sub(">", "", length_chr$Chr)  # Remove ">" from chromosome names



#### CHROMOSOME GRAPH


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
  title(main = plot_title, col.main = "black", font.main = 4)
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


pdf("results/plots/cnvs_cassification.pdf",width = 15, height = 8)


plot_all_chr_rows_len(length_chr, cnv_merged, "Copy Number Variations Across chromosomes")
# Add a legend to the plot

legend("bottomright",            # Position of the legend
       legend = c("gain", "loss", "gain/loss"),  # Labels for the legend
       fill = c("blue", "orange", "chartreuse3"), # Colors corresponding to the labels
       title = "CNV Cassification",    # Title for the legend
       cex = 0.8,             # Character expansion factor for the legend text
       bg = 'white')          # Background color of the legend box


dev.off()



