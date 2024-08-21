# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read the allMatchedMotifs.tsv file
cerercarie_allStemMotifs <- read.csv("C:/Users/Suryansh Raghuvanshi/Desktop/cercarie/allMatchedMotifs.tsv", sep="\t", header=TRUE)



# Manually input the MEME-JASPAR ID pairs
stem2_mejsp <- data.frame(
  motif_alt_id = c("STREME-1", "STREME-2", "STREME-3", "STREME-4", "STREME-5"),
  jaspar_id = c("MA0264.1", "MA0537.1", "MA1450.1", NA, "MA1444.1"),
  jaspar_name = c("ceh-22", "blmp-1", "lin-54", "", "cebp-1")
)

# Read the gene marker order file
stem_markerOrder <- read.table("C:/Users/Suryansh Raghuvanshi/Desktop/Cerercarie_Stem_Marker.txt", sep="", header=TRUE)

# Check the first few rows of the data
head(stem_markerOrder)

# Check column names of stem_markerOrder
colnames(stem_markerOrder)

# Ensure start and stop columns are numeric
allStemMotifs_new$start <- as.numeric(allStemMotifs_new$start)
allStemMotifs_new$stop <- as.numeric(allStemMotifs_new$stop)

# Check for any non-numeric values and handle them
allStemMotifs_new <- allStemMotifs_new[!is.na(allStemMotifs_new$start) & !is.na(allStemMotifs_new$stop), ]

# Merge the data based on MEME ids
allStemMotifs_new <- merge(allStemMotifs_new, stem2_mejsp, by="motif_alt_id", all.x=TRUE)
allStemMotifs_new <- merge(allStemMotifs_new, stem_markerOrder, by.x="sequence_name", by.y="Gene", all.x=TRUE)

# Check for successful merge and correct column names
head(allStemMotifs_new)
colnames(allStemMotifs_new)

# Add a column for JASPAR information
allStemMotifs_new$jaspar <- paste0(allStemMotifs_new$motif_alt_id, " (", allStemMotifs_new$jaspar_name, ")")

# Plotting with base lines for all genes
pstrul <- ggplot(data = allStemMotifs_new, aes(x = reorder(sequence_name, -Order))) +
  
  # Base layer for horizontal lines (gene range)
  geom_segment(aes(x = sequence_name, xend = sequence_name, y = -1000, yend = 0), 
               color = "grey", size = 0.5) +
  
  # Adding the motif segments
  geom_segment(aes(x = sequence_name, xend = sequence_name, y = start*(-1), yend = stop*(-1), color = jaspar), 
               size = 5) +
  
  coord_flip() +
  labs(x = "Genes", y = "Position") +
  ylim(-1000, 0) +
  theme_classic() +
  labs(color = "Motifs") +
  theme(axis.ticks.y = element_blank(), axis.line.y = element_blank())

# Print the updated plot
print(pstrul)
