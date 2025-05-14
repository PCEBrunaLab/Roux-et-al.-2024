#1.snATAC-seq Preprocessing

#Initial R Env with correct installation of packages - ArchR_Renv
renv::activate()

#Set Working Directory and Arch R parameters - i.e. where your lock file is located
setwd("~/Desktop/sn_clones/")

#Often run into issues if you try and parallalise processing so do not change these settings or try and introduce library(parallel)
addArchRThreads(threads = 4) 
addArchRGenome("hg38")

#Set hdf5 file locking to false to prevent downstream issues
set.seed(76)

## Create Arrow Files for ArchR ----
#Define Samples for Arrow file creation
# Create Arrow files from 10x cellranger output atac_fragment.tsv.gz files
Inputs <- list("~/Desktop/sn_clones/datafiles/Clone 2/atac_fragments.tsv.gz",
               "~/Desktop/sn_clones/datafiles/Clone 3/atac_fragments.tsv.gz",
               "~/Desktop/sn_clones/datafiles/Clone 4/atac_fragments.tsv.gz",
               "~/Desktop/sn_clones/datafiles/Clone 5/atac_fragments.tsv.gz",
               "~/Desktop/sn_clones/datafiles/Clone 7/atac_fragments.tsv.gz",
               "~/Desktop/sn_clones/datafiles/Clone 9/atac_fragments.tsv.gz",
               "~/Desktop/sn_clones/datafiles/Clone 11/atac_fragments.tsv.gz",
               "~/Desktop/sn_clones/datafiles/Clone 13/atac_fragments.tsv.gz")
inputFiles <- paste0(Inputs)

#Generate ArrowFiles
#If Seurat is broken then you will generate error messages saying HDF5 file is inaccessible but often due to broken Seurat
#For some reason everytime I now run this it breaks, copy in previously generated arrowfiles and see if we can go from there
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,
                               sampleNames = as.character(list("Clone_2", "Clone_3", "Clone_4", "Clone_5", "Clone_7", "Clone_9", "Clone_11", "Clone_13")),
                               minTSS = 4, #Dont set this too high because you can always increase later
                               minFrags = 1000,   
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
                               force = TRUE)

## Infer Doublets ----
#NB. For this to work correctly install ggplot older version from source
#Older and newer versions of Matrix will bug out the Doublet Scoring, remove Matrix and update to newest version then restart R session to use most updated version
doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                               knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbour search with doublet projection.
                               LSIMethod = 1) 

#Clone2 UMAP Projection R^2 - 0.77796; Correlation of UMAP Projection is below 0.9 (normally this is ~0.99); little heterogeneity in sample; setting as -1
#Clone3 UMAP Projection R^2 - 0.99535
#Clone4 UMAP Projection R^2 - 0.98431 ; Filtering 1 dims correlated > 0.75 to log10(depth + 1)
#Clone5 UMAP Projection R^2 - 0.9968
#Clone7 UMAP Projection R^2 - 0.99675
#Clone9 UMAP Projection R^2 - 0.97636
#Clone11 UMAP Projection R^2 - 0.97700
#Clone13 UMAP Projection R^2 - 0.63854; Correlation of UMAP Projection is below 0.9 (normally this is ~0.99); little heterogeneity in sample; setting as -1

## Create ArchR Project ----
project <- ArchRProject(ArrowFiles = ArrowFiles,
                        outputDirectory = "~/Desktop/sn_clones/project_clones/",
                        copyArrows = TRUE) #This is recommened so that if you modify the Arrow files you have an original copy for later usage.

## Display QC Metrics ----
meta.df <- getCellColData(project) #Extract metadata of cells
saveRDS(meta.df, file = "datafiles/metadata.RDS")

#Visualise Number of Fragments Vs TSS Enrichment
df <- getCellColData(project, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(x = df[,1], y = df[,2], 
             colorDensity = TRUE,
             continuousSet = "sambaNight",
             xlabel = "Log10 Unique Fragments",
             ylabel = "TSS Enrichment",
             xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
             ylim = c(0, quantile(df[,2], probs = 0.99))) + 
  geom_hline(yintercept = 4, lty = "dashed") + 
  geom_vline(xintercept = 3, lty = "dashed")

ggsave(plot = p, "plots/QCmetrics.pdf", width = 12, height = 7)

#Visualise Fragment Size Distributions
pdf("Fragment_size_distribution.pdf")
plotFragmentSizes(ArchRProj = project)
dev.off()

#Plot All Other QC Metrics
p1 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p2 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p3 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "ReadsInTSS",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p4 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "ReadsInPromoter",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p5 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "ReadsInBlacklist",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p6 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "PromoterRatio",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p7 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "NucleosomeRatio",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p8 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nMultiFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p9 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nMonoFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p10 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p11 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nDiFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p12 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "DoubletScore",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p13 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "DoubletEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

p14 <- plotGroups(
  ArchRProj = project, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "BlacklistRatio",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE)

plotPDF(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, name = "All_QC_metrics.pdf",
        ArchRProj = project, addDOC = FALSE, width = 4, height = 4)

## Filer out doublets ----
project <- filterDoublets(project, filterRatio = 1.2) 

#Filtering 5576 cells from ArchRProject!
#Clone_5 : 1119 of 9659 (11.6%)
#Clone_2 : 0 of 11060 (0%)
#Clone_3 : 498 of 6448 (7.7%)
#Clone_7 : 275 of 4790 (5.7%)
#Clone_9 : 1449 of 10989 (13.2%)
#Clone_4 : 447 of 6104 (7.3%)
#Clone_11 : 1788 of 12207 (14.6%)
#Clone_13 : 0 of 1404 (0%)

#Save ArchR Project
saveArchRProject(ArchRProj = project, outputDirectory = "~/Desktop/sn_clones/project_clones/", load = TRUE)

