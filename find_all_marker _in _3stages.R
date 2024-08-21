---
title:"findmarker_gene_across_stages"
author:"Tanushka Raghuvanshi"

# Load necessary libraries
library(DESeq2)
library(Seurat)  # Seurat v4.4.9
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)

#Miracidia####

# Load data
miracidia <- readRDS("C:/Users/Suryansh Raghuvanshi/Downloads/run1.combined.v10_integrated_55.rds")

# Set default assay to integrated
DefaultAssay(miracidia) <- "integrated"

# Rename cluster identities
miracidia <- RenameIdents(object = miracidia,
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
                          "18" = "Neuron 3")

# Save the identities to metadata columns
miracidia[["merge.ident"]] <- Idents(object = miracidia)
miracidia[["may.ident"]] <- Idents(object = miracidia)
table(Idents(miracidia))

# Visualize UMAP with cluster labels
p1 <- DimPlot(miracidia, reduction = "umap", label = TRUE)
p2 <- DimPlot(miracidia, reduction = "umap", label = TRUE, repel = TRUE)

# Merge all Stem cell identities into a single "Stem" label
miracidia <- RenameIdents(object = miracidia,
                          "Stem A" = "Stem",
                          "Stem B" = "Stem",
                          "Stem C" = "Stem",
                          "Stem D" = "Stem",
                          "Stem E" = "Stem",
                          "Stem F" = "Stem",
                          "Stem G" = "Stem")

# Identify all markers for each cluster using ROC test with min.pct >= 0.5
Allcell_marker <- FindAllMarkers(miracidia,
                                 only.pos = TRUE,
                                 min.pct = 0.5,
                                 logfc.threshold = 0.5,
                                 test.use = "roc",
                                 return.thresh = 0)

# Filter markers for the Stem cluster in the integrated dataset with different AUC thresholds
# AUC >= 0.8(5 marker stem)
stem_clusters_0.8 <- filter(Allcell_marker, cluster == "Stem" & myAUC >= 0.8) %>%
  select(gene) %>%
  distinct()

# AUC >= 0.7(32 marker stem)
stem_clusters_0.7 <- filter(Allcell_marker, cluster == "Stem" & myAUC >= 0.7) %>%
  select(gene) %>%
  distinct()
top_30_miracidia_stem_0.7 <- stem_cluster_0.7 %>% #top30 
  top_n(n = 30, wt = myAUC)

# AUC >= 0.6(69 marker stem)
stem_clusters_0.6 <- filter(Allcell_marker, cluster == "Stem" & myAUC >= 0.6) %>%
  select(gene) %>%
  distinct()

# Remove dashes in gene lists
stem_clusters_0.8$gene <- gsub('\\-', '_', stem_clusters_0.8$gene)
stem_clusters_0.7$gene <- gsub('\\-', '_', stem_clusters_0.7$gene)
stem_clusters_0.6$gene <- gsub('\\-', '_', stem_clusters_0.6$gene)

# Combine Stem A and Stem B into a new cluster Stem AB
miracidia <- RenameIdents(object = miracidia,
                          "Stem A" = "Stem AB",
                          "Stem B" = "Stem AB")

# Find all markers for the newly combined cluster
Allcell_marker <- FindAllMarkers(miracidia,
                                 only.pos = TRUE,
                                 min.pct = 0.5,
                                 logfc.threshold = 0.5,
                                 test.use = "roc",
                                 return.thresh = 0)

# Filter markers for Stem AB cluster with AUC >= 0.7#23
stem_AB_clusters <- filter(Allcell_marker, cluster == "Stem AB" & myAUC >= 0.7) %>%
  select(gene) %>%
  distinct()

# Remove dashes in gene list for Stem AB
stem_AB_clusters$gene <- gsub('\\-', '_', stem_AB_clusters$gene)

# Combine Stem C and Stem D into a new cluster Stem CD
miracidia <- RenameIdents(object = miracidia,
                          "Stem C" = "Stem CD",
                          "Stem D" = "Stem CD")

# Find all markers for the newly combined cluster
Allcell_marker <- FindAllMarkers(miracidia,
                                 only.pos = TRUE,
                                 min.pct = 0.5,
                                 logfc.threshold = 0.5,
                                 test.use = "roc",
                                 return.thresh = 0)

# Filter markers for Stem CD cluster with AUC >= 0.7#35
stem_CD_clusters <- filter(Allcell_marker, cluster == "Stem CD" & myAUC >= 0.7) %>%
  select(gene) %>%
  distinct()

# Remove dashes in gene list for Stem CD
stem_CD_clusters$gene <- gsub('\\-', '_', stem_CD_clusters$gene)

#Cercariae Stage####

# Load the Cercariae data
cercariae_26PC <- readRDS("C:/Users/Suryansh Raghuvanshi/AppData/Local/Temp/052d59d4-3721-418d-8660-80f8a7ded786_Tanushka_obj.zip.786/Tanushka_obj/cercariae_26PC.rds")

# Check the cluster identities
cluster_identities <- Idents(cercariae_26PC)
print(cluster_identities)
print(table(cluster_identities))  # Print table of cluster identities

# Visualize UMAP with cluster labels
p1 <- DimPlot(cercariae_26PC, reduction = "umap", label = TRUE)
p2 <- DimPlot(cercariae_26PC, reduction = "umap", label = FALSE, repel = TRUE)

# Find all markers for Cercariae using the ROC test
Allcell_marker_cercariae <- FindAllMarkers(
  cercariae_26PC,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.5,
  test.use = "roc",
  return.thresh = 0
)

# Filter markers for the Stem cluster with myAUC >= 0.8 in Cercariae (82 markers)
cercariae_stem_0.8 <- Allcell_marker_cercariae %>%
  filter(cluster == "Stem" & myAUC >= 0.8) %>%
  select(gene) %>%
  distinct()


# Filter markers for the Stem cluster with myAUC >= 0.7 in Cercariae (206 markers)
cercariae_stem_0.7 <- Allcell_marker_cercariae %>%
  filter(cluster == "Stem" & myAUC >= 0.7) %>%
  select(gene) %>%
  distinct()

top_30_cercariae_stem_0.7 <- cercariae_stem_0.7 %>% #top30 
  top_n(n = 30, wt = myAUC)

# Filter markers for the Stem cluster with myAUC >= 0.6 in Cercariae (350 markers)
cercariae_stem_0.6 <- Allcell_marker_cercariae %>%
  filter(cluster == "Stem" & myAUC >= 0.6) %>%
  select(gene) %>%
  distinct()

# Remove dashes in gene names for all the filtered stem gene lists
top_30_cercariae_stem_0.7 $gene <- gsub('\\-', '_', top_30_cercariae_stem_0.7 $gene)
cercariae_stem_0.8$gene <- gsub('\\-', '_', cercariae_stem_0.8$gene)
cercariae_stem_0.6$gene <- gsub('\\-', '_', cercariae_stem_0.6$gene)

#Sporocyts Stage####
  
# Load the Sporocysts data
sporocysts3 <- readRDS("C:/Users/Suryansh Raghuvanshi/Downloads/sporocysts3.rds")

# Find all markers for Sporocysts using the ROC test
Allcell_marker_sporocysts <- FindAllMarkers(
  sporocysts3,
  only.pos = TRUE,
  min.pct = 0.5,
  logfc.threshold = 0.5,
  test.use = "roc",
  return.thresh = 0
)


# Filter markers for the Stem/germinal cluster with AUC >= 0.8 (28 markers)
sporocysts_stem_0.8 <- Allcell_marker_sporocysts %>%
  filter(cluster == "Stem/germinal" & myAUC >= 0.8) %>%
  select(gene) %>%
  distinct()

# Filter markers for the Stem/germinal cluster with AUC >= 0.7 (189 markers)
sporocysts_stem_0.7 <- Allcell_marker_sporocysts %>%
  filter(cluster == "Stem/germinal" & myAUC >= 0.7) %>%
  select(gene) %>%
  distinct()

# Get top 30 markers for Stem/germinal cluster at AUC >= 0.7
top_30_sporocysts_stem_0.7 <- sporocysts_stem_0.7 %>%
  top_n(n = 30, wt = myAUC)

# Replace dashes in gene names with underscores for the top 30 markers
top_30_sporocysts_stem_0.7$gene <- gsub('\\-', '_', top_30_sporocysts_stem_0.7$gene)


# Filter markers for the Stem/germinal cluster with AUC >= 0.6 (417 markers)
sporocysts_stem_0.6 <- Allcell_marker_sporocysts %>%
  filter(cluster == "Stem/germinal" & myAUC >= 0.6) %>%
  select(gene) %>%
  distinct()

# Replace dashes in gene names with underscores for all filtered lists
sporocysts_stem_0.7$gene <- gsub('\\-', '_', sporocysts_stem_0.7$gene)
sporocysts_stem_0.8$gene <- gsub('\\-', '_', sporocysts_stem_0.8$gene)
sporocysts_stem_0.6$gene <- gsub('\\-', '_', sporocysts_stem_0.6$gene)
