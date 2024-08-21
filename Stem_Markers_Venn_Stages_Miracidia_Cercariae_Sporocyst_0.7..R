# Load necessary libraries
library(Seurat)      # For single-cell RNA-seq data analysis
library(dplyr)       # For data manipulation
library(pheatmap)    # For heatmap generation
library(VennDiagram) # For creating Venn diagrams
library(ggplot2)     # For data visualization

# Load the integrated dataset and specific stage datasets
run1.combined.v10 <- readRDS("C:/Users/Suryansh Raghuvanshi/Downloads/run1.combined.v10_integrated_55.rds")
cercariae_26PC <- readRDS("C:/Users/Suryansh Raghuvanshi/AppData/Local/Temp/052d59d4-3721-418d-8660-80f8a7ded786_Tanushka_obj.zip.786/Tanushka_obj/cercariae_26PC.rds")
sporocysts3 <- readRDS("C:/Users/Suryansh Raghuvanshi/Downloads/sporocysts3.rds")

# Set the default assay to "integrated" for the integrated dataset
DefaultAssay(run1.combined.v10) <- "integrated"

# Rename cluster identities in the integrated dataset to specific labels, including stem cell types
run1.combined.v10 <- RenameIdents(
  object = run1.combined.v10,
  "0" = "Muscle 1",
  "1" = "Stem C",
  "2" = "Stem D",
  "3" = "Stem F",
  "4" = "Stem A",
  "5" = "Stem B",
  "6" = "Neuron 1",
  "7" = "Muscle 2",
  "8" = "Neuron 4",
  "9" = "Parenchyma 1",
  "10" = "Stem E",
  "11" = "Parenchyma 2",
  "12" = "Protonephridia",
  "13" = "Tegument",
  "14" = "Neuron 2",
  "15" = "Neuron 5",
  "16" = "Stem G",
  "17" = "Ciliary plate",
  "18" = "Neuron 3"
)

# Store the current identities in metadata for future reference
run1.combined.v10[["merge.ident"]] <- Idents(object = run1.combined.v10)
run1.combined.v10[["may.ident"]] <- Idents(object = run1.combined.v10)

# Plot UMAP with cluster labels
DimPlot(run1.combined.v10, reduction = "umap", label = TRUE, repel = TRUE)

# Merge all Stem cell identities into a single "Stem" label
run1.combined.v10 <- RenameIdents(
  object = run1.combined.v10,
  "Stem A" = "Stem",
  "Stem B" = "Stem",
  "Stem C" = "Stem",
  "Stem D" = "Stem",
  "Stem E" = "Stem",
  "Stem F" = "Stem",
  "Stem G" = "Stem"
)

# Identify all markers for each cluster using the ROC test, focusing on markers with min.pct >= 0.5
Allcell_marker <- FindAllMarkers(
  run1.combined.v10,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.5,
  test.use = "roc",
  return.thresh = 0
)

# Filter markers for the Stem cluster in the integrated dataset with AUC >= 0.7
stem_clusters <- filter(Allcell_marker, cluster == "Stem" & myAUC >= 0.7)
miracidia_stem <- stem_clusters %>%
  select(gene) %>%
  distinct()

# Replace dashes in gene names with underscores
miracidia_stem$gene <- gsub('\\-', '_', miracidia_stem$gene)

# Plot UMAP for the Cercariae dataset with cluster labels
DimPlot(cercariae_26PC, reduction = "umap", label = TRUE)

# Identify all markers for each cluster in the Cercariae dataset using the ROC test
Allcell_marker_cercariae <- FindAllMarkers(
  cercariae_26PC,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.5,
  test.use = "roc",
  return.thresh = 0
)

# Filter markers for the Stem cluster in the Cercariae dataset with AUC >= 0.7
cercariae_stemgerm <- filter(Allcell_marker_cercariae, cluster == "Stem" & myAUC >= 0.7)
cercariae_stem <- cercariae_stemgerm %>%
  select(gene) %>%
  distinct()

# Replace dashes in gene names with underscores
cercariae_stem$gene <- gsub('\\-', '_', cercariae_stem$gene)

# Identify all markers for each cluster in the Sporocysts dataset using the ROC test
Allcell_marker_sporocysts <- FindAllMarkers(
  sporocysts3,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.5,
  test.use = "roc",
  return.thresh = 0
)

# Filter markers for the Stem/germinal cluster in the Sporocysts dataset with AUC >= 0.7
sporocysts_stemgerm <- filter(Allcell_marker_sporocysts, cluster == "Stem/germinal" & myAUC >= 0.7)
sporocysts_stem <- sporocysts_stemgerm %>%
  select(gene) %>%
  distinct()

# Replace dashes in gene names with underscores
sporocysts_stem$gene <- gsub('\\-', '_', sporocysts_stem$gene)

# Extract unique gene lists for each developmental stage
genes_miracidia <- unique(miracidia_stem$gene)
genes_cercariae <- unique(cercariae_stem$gene)
genes_sporocysts <- unique(sporocysts_stem$gene)

# Identify common and unique genes across the three stages
common_genes_all <- Reduce(intersect, list(genes_miracidia, genes_cercariae, genes_sporocysts))
common_genes_mir_cerca <- intersect(genes_miracidia, genes_cercariae)
common_genes_mir_sporo <- intersect(genes_miracidia, genes_sporocysts)
common_genes_cerca_sporo <- intersect(genes_cercariae, genes_sporocysts)

unique_genes_miracidia <- setdiff(genes_miracidia, union(genes_cercariae, genes_sporocysts))
unique_genes_cercariae <- setdiff(genes_cercariae, union(genes_miracidia, genes_sporocysts))
unique_genes_sporocysts <- setdiff(genes_sporocysts, union(genes_miracidia, genes_cercariae))

# Display the number of unique genes for each stage and the common genes across all stages
cat("Unique genes in Miracidia:", length(unique_genes_miracidia), "\n")
cat("Unique genes in Cercariae:", length(unique_genes_cercariae), "\n")
cat("Unique genes in Sporocysts:", length(unique_genes_sporocysts), "\n")
print(common_genes_all)

# Create a Venn Diagram to visualize the overlap of gene lists across the stages
venn.plot <- venn.diagram(
  x = list(
    Miracidia = genes_miracidia,
    Cercariae = genes_cercariae,
    Sporocysts = genes_sporocysts
  ),
  category.names = c("Miracidia", "Cercariae", "Sporocysts"),
  filename = NULL,
  fill = c("yellow", "green", "pink"),  # Set colors for each category
  alpha = 0.5,                           # Set transparency level
  cat.cex = 1.5,                         # Set size of category labels
  main.cex = 1.5                         # Set size of the main title
)

# Draw the Venn diagram
grid.draw(venn.plot)
