#1.snATAC-seq Preprocessing

#Initial R Env with correct installation of packages - ArchR_Renv
renv::activate()

#Set Working Directory and Arch R parameters
setwd("~/Desktop/sn_clones/")

addArchRThreads(threads = 4) 
addArchRGenome("hg38")
set.seed(76)

#Load in ArchR Project file
project <- readRDS("~/Desktop/sn_clones/project_clones/Save-ArchR-Project.rds")

## Dimensionality Reduction Using LSI ----
project <- addIterativeLSI(ArchRProj = project,
                           useMatrix = "TileMatrix",
                           name = "LSI_ATAC",
                           iterations = 2,
                           clusterParams = list(#See Seurat::FindClusters
                             resolution = c(0.2),
                             n.start = 10),
                           varFeatures = 25000,
                           dimsToUse = 1:30,
                           force = TRUE)

## Harmony Batch Correction ----
project <- addHarmony(ArchRProj = project,
                              reducedDims = "LSI_ATAC",
                              name = "Harmony_ATAC",
                              groupBy = "Sample",
                              force = TRUE)

## Cluster Identification ----
project <- addClusters(input = project,
                           reducedDims = "Harmony_ATAC",
                           method = "Seurat",
                           name = "Clusters_Harmony_ATAC",
                           resolution = 0.6,
                           force = TRUE)

project <- addClusters(input = project,
                               reducedDims = "LSI_ATAC",
                               method = "Seurat",
                               name = "Clusters_LSI_ATAC",
                               resolution = 0.6,
                               force = TRUE)

## UMAP Visualisation ----
project <- addUMAP(ArchRProj = project, 
                       reducedDims = "Harmony_ATAC", 
                       name = "UMAP_Harmony", 
                       nNeighbors = 30, 
                       minDist = 0.5, 
                       metric = "cosine",
                       force = TRUE)

project <- addUMAP(ArchRProj = project, 
                           reducedDims = "LSI_ATAC", 
                           name = "UMAP_LSI", 
                           nNeighbors = 30, 
                           minDist = 0.5, 
                           metric = "cosine",
                           force = TRUE)

UMAP_Harmony_df <- getEmbedding(ArchRProj = project, embedding = "UMAP_Harmony", returnDF = TRUE)
UMAP_LSI_df <- getEmbedding(ArchRProj = project, embedding = "UMAP_LSI", returnDF = TRUE)
saveRDS(UMAP_Harmony_df, "datafiles/UMAP_Harmony_coordinates.RDS")
saveRDS(UMAP_LSI_df, "datafiles/UMAP_LSI_coordinates.RDS")

pdf("plots/ArchR_processing_results.pdf")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "UMAP_Harmony")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters_Harmony_ATAC", embedding = "UMAP_Harmony")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP_Harmony")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "UMAP_LSI")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters_LSI_ATAC", embedding = "UMAP_LSI")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAP_LSI")
dev.off()

## Get marker genes for clusters of LSI and harmony data ----
#NB. if there are issues with GeneMatrix, re-add with the code lines below
#Add GeneScoreMatrix for ArchR project
#getAvailableMatrices(ArchRProj = project)
#project <- addGeneScoreMatrix(project, force = TRUE)

#Identify marker genes for clusters for LSI
markersGS_LSI <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters_LSI_ATAC",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

#Isolate significant marker genes
#Matching previously used FDR and logsFC cut offs i.e. loose thresholds to visualise as much data as possible
markerList_LSI <- getMarkers(markersGS_LSI, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

#Identify values to look at marker genes for clusters
heatmapGS_LSI <- plotMarkerHeatmap(
  seMarker = markersGS_LSI, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  transpose = TRUE)

#Plot heatmap for each clustering method
install.packages("magick")
library(magick)

plotPDF(heatmapGS_LSI, name = "GeneScores-Marker-Heatmap-LSI", width = 14, height = 6, ArchRProj = project, addDOC = FALSE)

#Convert Summarised heatmap of marker genes to dataframe for saving
markerList_LSI <- as.data.frame(markerList_LSI)
write.csv(markerList_LSI, "datafiles/LSI_cluster_marker_genes.csv")


## Save ArchR Project ----
saveArchRProject(ArchRProj = project, outputDirectory = "~/Desktop/sn_clones/project_clones_working", load = TRUE)
meta.df <- getCellColData(project)
saveRDS(meta.df, "datafiles/metadata_clustering.RDS")

