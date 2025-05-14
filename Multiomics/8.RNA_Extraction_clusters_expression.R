##Imputation 

#Install extra packages
required.packages <- c("dplyr", "tidyr", "Seurat", "reticulate")

if (!all(required.packages %in% installed.packages()[,"Package"])){
  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]], force = T)
}

library(dplyr)
library(tidyr)
library(reticulate)
library(Seurat)
library(paletteer)
set.seed(76)

setwd("~/Desktop/sn_clones/")

#Define plot variables for condition
group.cols <- paletteer_d("rcartocolor::Prism")
names(group.cols) <- c("Clone2", "Clone3", "Clone4", "Clone5", "Clone7",
                       "Clone9", "Clone11", "Clone13")

#Non-batch corrected data
#Load in Data files
sample <- readRDS("~/OneDrive/PhD/Thesis/Thesis Plots/Results Chapter 3/Sub Chapter 1/datafiles/nb_seurat_AMT.rds")

#Extract matrix of normalised expression data
exp_mat <- as.matrix(GetAssayData(sample, assay = "SCT", slot = "data"))

#Retrieve PCA cell embeddings
pca_values <- sample@reductions[["pca"]]@cell.embeddings

## Cluster-based Mean Expression Table ----
DefaultAssay(sample) <- "SCT"

#We have already defined clusters and neighbours etc so lets skip to extracting this information
cluster_order <- paste0("C", as.character(sample$seurat_clusters))
sample$seurat_clusters <- paste0("C", as.character(sample$seurat_clusters))
cells_in_cluster <- sample$seurat_clusters; names(cells_in_cluster) <- rownames(sample@meta.data)
saveRDS(cells_in_cluster, "datafiles/cells_in_clusters.RDS")

DimPlot(sample, group.by = "seurat_clusters", label = TRUE) +
  NoLegend()

## Create Expression Table ----
#Create matrix with genes as rownames and seurat cluster names as colnames
exp_mat_clusters <- as.data.frame(matrix(NA, nrow = nrow(exp_mat)))
rownames(exp_mat_clusters) <- rownames(exp_mat)
exp_mat_clusters[, unique(cluster_order)] = NA
exp_mat_clusters = subset(exp_mat_clusters, select = -c(V1) ) 

#For each cell cluster compute the average expression for each gene
for(cluster in colnames(exp_mat_clusters)){ 
  cells_in_cluster <- rownames(sample@meta.data[which(sample$seurat_clusters == cluster),])
  exp_mat_clusters[,cluster] <- apply(exp_mat[,cells_in_cluster], 1, mean)
}

saveRDS(exp_mat_clusters, "datafiles/exp_mat_clusters.RDS")

## Add Cluster Annotation ---- 
cluster_annot <- data.frame(cluster = unique(cluster_order), sample = NA, mean_PC1 = NA, mean_PC2 = NA, stringsAsFactors = FALSE)
rownames(cluster_annot) <- cluster_annot$cluster

for(cluster in cluster_annot$cluster){ 
  # Add sample of origin
  cluster_annot[cluster,"sample"] <- paste(unique(sample@meta.data[which(sample$seurat_clusters==cluster),"orig.ident"]), collapse=", ")
  # Add mean PC
  cluster_annot[cluster, "mean_PC1"] <- mean(pca_values[rownames(sample@meta.data[which(sample$seurat_clusters == cluster),]), "PC_1"])
  cluster_annot[cluster, "mean_PC2"] <- mean(pca_values[rownames(sample@meta.data[which(sample$seurat_clusters == cluster),]), "PC_2"])
}

## Project Cluster Mean Expression on Individual Cells ----
exp_mat_clusters_per_cell <- exp_mat
for(cluster in colnames(exp_mat_clusters)){
  cells_in_cluster <- rownames(sample@meta.data[which(sample$seurat_clusters == cluster),])
  exp_mat_clusters_per_cell[,cells_in_cluster] <- exp_mat_clusters[,cluster]
}

saveRDS(exp_mat_clusters_per_cell, "~/OneDrive/PhD/Thesis/Thesis Plots/Results Chapter 3/Sub Chapter 1/datafiles/exp_mat_clusters_per_cell.RDS")

