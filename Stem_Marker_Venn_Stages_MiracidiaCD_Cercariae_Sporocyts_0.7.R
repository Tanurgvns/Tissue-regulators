# Load required libraries
library(Seurat)
library(dplyr)
library(pheatmap)
library(VennDiagram)
library(ggplot2)

# Load data
run1.combined.v10 <- readRDS("C:/Users/Suryansh Raghuvanshi/Downloads/run1.combined.v10_integrated_55.rds")
cercariae_26PC <- readRDS("C:/Users/Suryansh Raghuvanshi/AppData/Local/Temp/052d59d4-3721-418d-8660-80f8a7ded786_Tanushka_obj.zip.786/Tanushka_obj/cercariae_26PC.rds")
sporocysts3 <- readRDS("C:/Users/Suryansh Raghuvanshi/Downloads/sporocysts3.rds")

DefaultAssay(run1.combined.v10) <- "integrated"

# Rename cluster 
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


run1.combined.v10[["merge.ident"]] <- Idents(object = run1.combined.v10)
run1.combined.v10[["may.ident"]] <- Idents(object = run1.combined.v10)
DimPlot(run1.combined.v10, reduction = "umap", label = TRUE, repel = TRUE)

# Combine Stem C and Stem D into a new cluster Stem CD
run1.combined.v10 <- RenameIdents(
  object = run1.combined.v10,
  "Stem C" = "Stem CD",
  "Stem D" = "Stem CD"
)

# Find all markers using ROC test with min.pct=0.5
Allcell_marker <- FindAllMarkers(
  run1.combined.v10,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.5,
  test.use = "roc",
  return.thresh = 0
)

# Filter markers for Stem CD cluster only with myAUC >= 0.7
stem_CD_clusters <- filter(Allcell_marker, cluster == "StemCD" & myAUC >= 0.7)
miracidia_stem <- stem_CD_clusters %>%
  select(gene) %>%
  distinct()

# Remove dashes in gene list
miracidia_stem$gene <- gsub('\\-', '_', miracidia_stem$gene)

# Plot UMAP for Cercariae
DimPlot(cercariae_26PC, reduction = "umap", label = TRUE)

# Find all markers for Cercariae using ROC test
Allcell_marker_cercariae <- FindAllMarkers(
  cercariae_26PC,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.5,
  test.use = "roc",
  return.thresh = 0
)

# Filter markers for Stem cluster only with myAUC >= 0.7 in Cercariae
cercariae_stemgerm <- filter(Allcell_marker_cercariae, cluster == "Stem" & myAUC >= 0.7)
cercariae_stem <- cercariae_stemgerm %>%
  select(gene) %>%
  distinct()
cercariae_stem$gene <- gsub('\\-', '_', cercariae_stem$gene)

# Find all markers for Sporocysts using ROC test
Allcell_marker_sporocysts <- FindAllMarkers(
  sporocysts3,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.5,
  test.use = "roc",
  return.thresh = 0
)

# Filter markers for Stem/germinal cluster only with myAUC >= 0.7 in Sporocysts
sporocysts_stemgerm <- filter(Allcell_marker_sporocysts, cluster == "Stem/germinal" & myAUC >= 0.7)
sporocysts_stem <- sporocysts_stemgerm %>%
  select(gene) %>%
  distinct()
sporocysts_stem$gene <- gsub('\\-', '_', sporocysts_stem$gene)

# Extract unique gene lists for each stage
genes_miracidia <- unique(miracidia_stem$gene)
genes_cercariae <- unique(cercariae_stem$gene)
genes_sporocysts <- unique(sporocysts_stem$gene)

# Find common and unique genes
common_genes_all <- Reduce(intersect, list(genes_miracidia, genes_cercariae, genes_sporocysts))
common_genes_mir_cerca <- intersect(genes_miracidia, genes_cercariae)
common_genes_mir_sporo <- intersect(genes_miracidia, genes_sporocysts)
common_genes_cerca_sporo <- intersect(genes_cercariae, genes_sporocysts)

unique_genes_miracidia <- setdiff(genes_miracidia, union(genes_cercariae, genes_sporocysts))
unique_genes_cercariae <- setdiff(genes_cercariae, union(genes_miracidia, genes_sporocysts))
unique_genes_sporocysts <- setdiff(genes_sporocysts, union(genes_miracidia, genes_cercariae))

cat("Unique genes in Miracidia:", length(unique_genes_miracidia), "\n")
cat("Unique genes in Cercariae:", length(unique_genes_cercariae), "\n")
cat("Unique genes in Sporocysts:", length(unique_genes_sporocysts), "\n")
print(common_genes_all)

# Create a Venn Diagram with colors
venn.plot <- venn.diagram(
  x = list(
    Miracidia = genes_miracidia,
    Cercariae = genes_cercariae,
    Sporocysts = genes_sporocysts
  ),
  category.names = c("Miracidia Stem CD", "Cercariae", "Sporocysts"),
  filename = NULL,
  fill = c("yellow", "green", "pink"), # Set colors for each category
  alpha = 0.5, # Transparency level
  cat.cex = 1.5, # Size of category labels
  main.cex = 1.5 # Size of the main title
)

# Draw the Venn diagram
grid.draw(venn.plot)
