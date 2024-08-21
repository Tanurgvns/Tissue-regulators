# Load necessary libraries
library(ggplot2)

# Read the combined motifs file
miracidia_allMatchedMotifs <- read.table("C:/Users/Suryansh Raghuvanshi/Downloads/Miracidia/allMatchedMotifs.tsv", sep="\t", header=TRUE)


# Manually input the MEME-JASPAR ID pairs 
stem2_mejsp <- data.frame(
  motif_alt_id = c("STREME-1", "STREME-2", "STREME-3", "STREME-4", "STREME-5"),
  jaspar_id = c( NA, "MA1699.1", "MA0264.1", "MA1704.1"),
  jaspar_name = c("", "ceh-38", "ceh-22", "", "zip-8")
)

# Read the gene marker order file
stem_markerOrder <- read.table("C:/Users/Suryansh Raghuvanshi/Downloads/Miraidia_ marker.txt", sep="", header=TRUE)

# Merge the data based on MEME ids
allStemMotifs_new <- merge(allMatchedMotifs, stem2_mejsp, by="motif_alt_id", all.x=TRUE)
allStemMotifs_new <- merge(allStemMotifs_new, stem_markerOrder, by.x="sequence_name", by.y="Gene", all.x=TRUE)

# Ensure start and stop columns are numeric
allStemMotifs_new$start <- as.numeric(allStemMotifs_new$start)
allStemMotifs_new$stop <- as.numeric(allStemMotifs_new$stop)

# Extract only the gene ID from sequence_name
allStemMotifs_new$gene_id <- sapply(strsplit(allStemMotifs_new$sequence_name, "\\|"), `[`, 1)

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
