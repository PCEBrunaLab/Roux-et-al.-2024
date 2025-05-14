#1.snATAC-seq Preprocessing

#Work outside of REnv in R version 4.4 to use up to date packages for Signac processing

#Set Working Directory and Arch R parameters
setwd("~/OneDrive/PhD/Thesis/Thesis Plots/Results Chapter 3/Sub Chapter 1/")
#NB. for some reason setwd is not changing so just add ../ to file paths 

#Install extra packages
required.packages <- c("dplyr", "tidyr", "Signac", "Seurat", "GenomeInfoDb", "EnsDb.Hsapiens.v86",
                       "BSgenome.Hsapiens.UCSC.hg38", "patchwork", "MASS", "viridis", "tidyverse",
                       "DropletUtils")

if (!all(required.packages %in% installed.packages()[,"Package"])){
  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]])
}

library(dplyr)
library(tidyr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(MASS)
library(viridis)
library(tidyverse)
library(DropletUtils)
set.seed(76)

#Load in data for each sample
counts.Clone2 <- Read10X(data.dir = "datafiles/cellranger_outputs/Clone2/filtered_feature_bc_matrix/")
counts.Clone2 <- counts.Clone2[[2]]

#Read in all necessary cellranger output files per clone
read_all_clone_data <- function(clone_names = c("Clone 2", "Clone 3", "Clone 4", "Clone 5",
                                                "Clone 7", "Clone 9", "Clone 11", "Clone 13"),
                                base_dir = "datafiles/cellranger outputs") {
  # Use lapply to read each dataset
  data_list <- lapply(clone_names, function(clone) {
    # Construct the full path to the filtered_feature_bc_matrix folder
    file_path <- file.path(base_dir, clone, "filtered_feature_bc_matrix")
    Read10X(data.dir = file_path)
  })
  # Name each list element by its clone ID for easier reference
  names(data_list) <- clone_names
  return(data_list)
}

all_clone_data <- read_all_clone_data()

#Remove RNA assay from each object
#Add sample variable to all objects
for(clone in names(all_clone_data)) {
  sce_obj <- all_clone_data[[clone]]
  
  #Remove RNA assay from each object
  sce_obj <- sce_obj[[2]]
  
  # Save the modified object back into the list
  all_clone_data[[clone]] <- sce_obj
}


#Gene annotation for TSS enrichment
BiocManager::install("biovizBase", force = T)
library(biovizBase)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

#Load cells to keep after Arch QC 
archr_cells <- readRDS("datafiles/metadata_clustering.RDS")
archr_cells <- gsub("#", "_", rownames(archr_cells))

## Generation Signac Objects ----
#Create Chromatin assay in Signac for each sample
chrom_assay.Clone2 <- CreateChromatinAssay(counts = counts.Clone2,
                                    sep = c(":", "-"),
                                    genome = NULL, 
                                    fragments = "cellranger_outputs/Clone2/atac_fragments.tsv.gz",
                                    min.cells = 0, # default is 10
                                    min.features = 0, # default is 200
                                    annotation=annotation)

chrom_assay.Clone3 <- CreateChromatinAssay(counts = counts.Clone3,
                                           sep = c(":", "-"),
                                           genome = NULL, 
                                           fragments = "cellranger_outputs/Clone3/atac_fragments.tsv.gz",
                                           min.cells = 0, # default is 10
                                           min.features = 0, # default is 200
                                           annotation=annotation)

chrom_assay.Clone4 <- CreateChromatinAssay(counts = counts.Clone4,
                                           sep = c(":", "-"),
                                           genome = NULL, 
                                           fragments = "cellranger_outputs/Clone4/atac_fragments.tsv.gz",
                                           min.cells = 0, # default is 10
                                           min.features = 0, # default is 200
                                           annotation=annotation)

chrom_assay.Clone5 <- CreateChromatinAssay(counts = counts.Clone5,
                                           sep = c(":", "-"),
                                           genome = NULL, 
                                           fragments = "cellranger_outputs/Clone5/atac_fragments.tsv.gz",
                                           min.cells = 0, # default is 10
                                           min.features = 0, # default is 200
                                           annotation=annotation)

chrom_assay.Clone7 <- CreateChromatinAssay(counts = counts.Clone7,
                                           sep = c(":", "-"),
                                           genome = NULL, 
                                           fragments = "cellranger_outputs/Clone7/atac_fragments.tsv.gz",
                                           min.cells = 0, # default is 10
                                           min.features = 0, # default is 200
                                           annotation=annotation)

chrom_assay.Clone9 <- CreateChromatinAssay(counts = counts.Clone9,
                                           sep = c(":", "-"),
                                           genome = NULL, 
                                           fragments = "cellranger_outputs/Clone9/atac_fragments.tsv.gz",
                                           min.cells = 0, # default is 10
                                           min.features = 0, # default is 200
                                           annotation=annotation)

chrom_assay.Clone11 <- CreateChromatinAssay(counts = counts.Clone11,
                                           sep = c(":", "-"),
                                           genome = NULL, 
                                           fragments = "cellranger_outputs/Clone11/atac_fragments.tsv.gz",
                                           min.cells = 0, # default is 10
                                           min.features = 0, # default is 200
                                           annotation=annotation)

chrom_assay.Clone13 <- CreateChromatinAssay(counts = counts.Clone13,
                                           sep = c(":", "-"),
                                           genome = NULL, 
                                           fragments = "cellranger_outputs/Clone13/atac_fragments.tsv.gz",
                                           min.cells = 0, # default is 10
                                           min.features = 0, # default is 200
                                           annotation=annotation)



#Create Seurat Objects from chromatin assays
sample_atac.Clone2 <- CreateSeuratObject(counts = chrom_assay.Clone2,
                                  assay = "peaks",
                                  meta.data = NULL)

sample_atac.Clone3 <- CreateSeuratObject(counts = chrom_assay.Clone3,
                                         assay = "peaks",
                                         meta.data = NULL)

sample_atac.Clone4 <- CreateSeuratObject(counts = chrom_assay.Clone4,
                                         assay = "peaks",
                                         meta.data = NULL)

sample_atac.Clone5 <- CreateSeuratObject(counts = chrom_assay.Clone5,
                                         assay = "peaks",
                                         meta.data = NULL)

sample_atac.Clone7 <- CreateSeuratObject(counts = chrom_assay.Clone7,
                                         assay = "peaks",
                                         meta.data = NULL)

sample_atac.Clone9 <- CreateSeuratObject(counts = chrom_assay.Clone9,
                                         assay = "peaks",
                                         meta.data = NULL)

sample_atac.Clone11 <- CreateSeuratObject(counts = chrom_assay.Clone11,
                                         assay = "peaks",
                                         meta.data = NULL)

sample_atac.Clone13 <- CreateSeuratObject(counts = chrom_assay.Clone13,
                                         assay = "peaks",
                                         meta.data = NULL)


rm(chrom_assay.Clone2, chrom_assay.Clone3, chrom_assay.Clone4 ,chrom_assay.Clone5,
   chrom_assay.Clone7, chrom_assay.Clone9, chrom_assay.Clone11, chrom_assay.Clone13,
   counts.Clone2, counts.Clone3, counts.Clone4 ,counts.Clone5,
   counts.Clone7, counts.Clone9, counts.Clone11, counts.Clone13)

#Add CellID and sample to meta data
sample_atac.Clone2$Sample <- "Clone_2"
sample_atac.Clone2$CellID <- paste0(sample_atac.Clone2$Sample, "_" ,colnames(sample_atac.Clone2))

sample_atac.Clone3$Sample <- "Clone_3"
sample_atac.Clone3$CellID <- paste0(sample_atac.Clone3$Sample, "_" ,colnames(sample_atac.Clone3))

sample_atac.Clone4$Sample <- "Clone_4"
sample_atac.Clone4$CellID <- paste0(sample_atac.Clone4$Sample, "_" ,colnames(sample_atac.Clone4))

sample_atac.Clone5$Sample <- "Clone_5"
sample_atac.Clone5$CellID <- paste0(sample_atac.Clone5$Sample, "_" ,colnames(sample_atac.Clone5))

sample_atac.Clone7$Sample <- "Clone_7"
sample_atac.Clone7$CellID <- paste0(sample_atac.Clone7$Sample, "_" ,colnames(sample_atac.Clone7))

sample_atac.Clone9$Sample <- "Clone_9"
sample_atac.Clone9$CellID <- paste0(sample_atac.Clone9$Sample, "_" ,colnames(sample_atac.Clone9))

sample_atac.Clone11$Sample <- "Clone_11"
sample_atac.Clone11$CellID <- paste0(sample_atac.Clone11$Sample, "_" ,colnames(sample_atac.Clone11))

sample_atac.Clone13$Sample <- "Clone_13"
sample_atac.Clone13$CellID <- paste0(sample_atac.Clone13$Sample, "_" ,colnames(sample_atac.Clone13))


#Save ATAC datasets
saveRDS(sample_atac.Clone2, "datafiles/Clone2_signac_atac.RDS")
saveRDS(sample_atac.Clone3, "datafiles/Clone3_signac_atac.RDS")
saveRDS(sample_atac.Clone4, "datafiles/Clone4_signac_atac.RDS")
saveRDS(sample_atac.Clone5, "datafiles/Clone5_signac_atac.RDS")
saveRDS(sample_atac.Clone7, "datafiles/Clone7_signac_atac.RDS")
saveRDS(sample_atac.Clone9, "datafiles/Clone9_signac_atac.RDS")
saveRDS(sample_atac.Clone11, "datafiles/Clone11_signac_atac.RDS")
saveRDS(sample_atac.Clone13, "datafiles/Clone13_signac_atac.RDS")

gc()

## QC Metrics Prior to QC ----
#Load in datafiles
sample_atac.Clone2 <- readRDS("datafiles/Clone2_signac_atac.RDS")
sample_atac.Clone3 <- readRDS("datafiles/Clone3_signac_atac.RDS")
sample_atac.Clone4 <- readRDS("datafiles/Clone4_signac_atac.RDS")
sample_atac.Clone5 <- readRDS("datafiles/Clone5_signac_atac.RDS")
sample_atac.Clone7 <- readRDS("datafiles/Clone7_signac_atac.RDS")
sample_atac.Clone9 <- readRDS("datafiles/Clone9_signac_atac.RDS")
sample_atac.Clone11 <- readRDS("datafiles/Clone11_signac_atac.RDS")
sample_atac.Clone13 <- readRDS("datafiles/Clone13_signac_atac.RDS")

#Calculate Total Fragment Counts per Cell Barcode
fragment_counts.Clone2 <- CountFragments(fragments = "cellranger_outputs/Clone2/atac_fragments.tsv.gz")
fragment_counts.Clone3 <- CountFragments(fragments = "cellranger_outputs/Clone3/atac_fragments.tsv.gz")
fragment_counts.Clone4 <- CountFragments(fragments = "cellranger_outputs/Clone4/atac_fragments.tsv.gz")
fragment_counts.Clone5 <- CountFragments(fragments = "cellranger_outputs/Clone5/atac_fragments.tsv.gz")
fragment_counts.Clone7 <- CountFragments(fragments = "cellranger_outputs/Clone7/atac_fragments.tsv.gz")
fragment_counts.Clone9 <- CountFragments(fragments = "cellranger_outputs/Clone9/atac_fragments.tsv.gz")
fragment_counts.Clone11 <- CountFragments(fragments = "cellranger_outputs/Clone11/atac_fragments.tsv.gz")
fragment_counts.Clone13 <- CountFragments(fragments = "cellranger_outputs/Clone13/atac_fragments.tsv.gz")

fragment_counts.Clone2$CellID <- paste0("Clone_2", "_", fragment_counts.Clone2$CB)
fragment_counts.Clone3$CellID <- paste0("Clone_3", "_", fragment_counts.Clone3$CB)
fragment_counts.Clone4$CellID <- paste0("Clone_4", "_", fragment_counts.Clone4$CB)
fragment_counts.Clone5$CellID <- paste0("Clone_5", "_", fragment_counts.Clone5$CB)
fragment_counts.Clone7$CellID <- paste0("Clone_7", "_", fragment_counts.Clone7$CB)
fragment_counts.Clone9$CellID <- paste0("Clone_9", "_", fragment_counts.Clone9$CB)
fragment_counts.Clone11$CellID <- paste0("Clone_11", "_", fragment_counts.Clone11$CB)
fragment_counts.Clone13$CellID <- paste0("Clone_13", "_", fragment_counts.Clone13$CB)

fragment_counts.Clone2$Sample <- "Clone2"
fragment_counts.Clone3$Sample <- "Clone3"
fragment_counts.Clone4$Sample <- "Clone4"
fragment_counts.Clone5$Sample <- "Clone5"
fragment_counts.Clone7$Sample <- "Clone7"
fragment_counts.Clone9$Sample <- "Clone9"
fragment_counts.Clone11$Sample <- "Clone11"
fragment_counts.Clone13$Sample <- "Clone13"

rownames(fragment_counts.Clone2) <- fragment_counts.Clone2$CellID
rownames(fragment_counts.Clone3) <- fragment_counts.Clone3$CellID
rownames(fragment_counts.Clone4) <- fragment_counts.Clone4$CellID
rownames(fragment_counts.Clone5) <- fragment_counts.Clone5$CellID
rownames(fragment_counts.Clone7) <- fragment_counts.Clone7$CellID
rownames(fragment_counts.Clone9) <- fragment_counts.Clone9$CellID
rownames(fragment_counts.Clone11) <- fragment_counts.Clone11$CellID
rownames(fragment_counts.Clone13) <- fragment_counts.Clone13$CellID

fragment_counts.Clone2 <- fragment_counts.Clone2[sample_atac.Clone2$CellID,]
fragment_counts.Clone3 <- fragment_counts.Clone3[sample_atac.Clone3$CellID,]
fragment_counts.Clone4 <- fragment_counts.Clone4[sample_atac.Clone4$CellID,]
fragment_counts.Clone5 <- fragment_counts.Clone5[sample_atac.Clone5$CellID,]
fragment_counts.Clone7 <- fragment_counts.Clone7[sample_atac.Clone7$CellID,]
fragment_counts.Clone9 <- fragment_counts.Clone9[sample_atac.Clone9$CellID,]
fragment_counts.Clone11 <- fragment_counts.Clone11[sample_atac.Clone11$CellID,]
fragment_counts.Clone13 <- fragment_counts.Clone13[sample_atac.Clone13$CellID,]

sample_atac.Clone2$nb_fragments <- fragment_counts.Clone2$frequency_count # add to signac object
sample_atac.Clone3$nb_fragments <- fragment_counts.Clone3$frequency_count # add to signac object
sample_atac.Clone4$nb_fragments <- fragment_counts.Clone4$frequency_count # add to signac object
sample_atac.Clone5$nb_fragments <- fragment_counts.Clone5$frequency_count # add to signac object
sample_atac.Clone7$nb_fragments <- fragment_counts.Clone7$frequency_count # add to signac object
sample_atac.Clone9$nb_fragments <- fragment_counts.Clone9$frequency_count # add to signac object
sample_atac.Clone11$nb_fragments <- fragment_counts.Clone11$frequency_count # add to signac object
sample_atac.Clone13$nb_fragments <- fragment_counts.Clone13$frequency_count # add to signac object

#compute Nucleosome Signal Score Per Cell
sample_atac.Clone2 <- NucleosomeSignal(object = sample_atac.Clone2)
sample_atac.Clone3 <- NucleosomeSignal(object = sample_atac.Clone3)
sample_atac.Clone4 <- NucleosomeSignal(object = sample_atac.Clone4)
sample_atac.Clone5 <- NucleosomeSignal(object = sample_atac.Clone5)
sample_atac.Clone7 <- NucleosomeSignal(object = sample_atac.Clone7)
sample_atac.Clone9 <- NucleosomeSignal(object = sample_atac.Clone9)
sample_atac.Clone11 <- NucleosomeSignal(object = sample_atac.Clone11)
sample_atac.Clone13 <- NucleosomeSignal(object = sample_atac.Clone13)

#Compute TSS Enrichment Score Per Cell
sample_atac.Clone2 <- TSSEnrichment(object = sample_atac.Clone2, fast = FALSE)
sample_atac.Clone3 <- TSSEnrichment(object = sample_atac.Clone3, fast = FALSE)
sample_atac.Clone4 <- TSSEnrichment(object = sample_atac.Clone4, fast = FALSE)
sample_atac.Clone5 <- TSSEnrichment(object = sample_atac.Clone5, fast = FALSE)
sample_atac.Clone7 <- TSSEnrichment(object = sample_atac.Clone7, fast = FALSE)
sample_atac.Clone9 <- TSSEnrichment(object = sample_atac.Clone9, fast = FALSE)
sample_atac.Clone11 <- TSSEnrichment(object = sample_atac.Clone11, fast = FALSE)
sample_atac.Clone13 <- TSSEnrichment(object = sample_atac.Clone13, fast = FALSE)

## Visualisation of Quality Control Metrics ----
df_metrics.Clone2 <- sample_atac.Clone2@meta.data[,c("nb_fragments", "TSS.enrichment")]
df_metrics.Clone3 <- sample_atac.Clone3@meta.data[,c("nb_fragments", "TSS.enrichment")]
df_metrics.Clone4 <- sample_atac.Clone4@meta.data[,c("nb_fragments", "TSS.enrichment")]
df_metrics.Clone5 <- sample_atac.Clone5@meta.data[,c("nb_fragments", "TSS.enrichment")]
df_metrics.Clone7 <- sample_atac.Clone7@meta.data[,c("nb_fragments", "TSS.enrichment")]
df_metrics.Clone9 <- sample_atac.Clone9@meta.data[,c("nb_fragments", "TSS.enrichment")]
df_metrics.Clone11 <- sample_atac.Clone11@meta.data[,c("nb_fragments", "TSS.enrichment")]
df_metrics.Clone13 <- sample_atac.Clone13@meta.data[,c("nb_fragments", "TSS.enrichment")]

df_metrics.Clone2$log10_nb_fragments <- log10(df_metrics.Clone2$nb_fragments)
df_metrics.Clone3$log10_nb_fragments <- log10(df_metrics.Clone3$nb_fragments)
df_metrics.Clone4$log10_nb_fragments <- log10(df_metrics.Clone4$nb_fragments)
df_metrics.Clone5$log10_nb_fragments <- log10(df_metrics.Clone5$nb_fragments)
df_metrics.Clone7$log10_nb_fragments <- log10(df_metrics.Clone7$nb_fragments)
df_metrics.Clone9$log10_nb_fragments <- log10(df_metrics.Clone9$nb_fragments)
df_metrics.Clone11$log10_nb_fragments <- log10(df_metrics.Clone11$nb_fragments)
df_metrics.Clone13$log10_nb_fragments <- log10(df_metrics.Clone13$nb_fragments)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

df_metrics.Clone2$density <- get_density(df_metrics.Clone2$log10_nb_fragments, df_metrics.Clone2$TSS.enrichment, n = 5000) # n: create a square n by n grid to compute density
df_metrics.Clone3$density <- get_density(df_metrics.Clone3$log10_nb_fragments, df_metrics.Clone3$TSS.enrichment, n = 5000)
df_metrics.Clone4$density <- get_density(df_metrics.Clone4$log10_nb_fragments, df_metrics.Clone4$TSS.enrichment, n = 5000)
df_metrics.Clone5$density <- get_density(df_metrics.Clone5$log10_nb_fragments, df_metrics.Clone5$TSS.enrichment, n = 5000)
df_metrics.Clone7$density <- get_density(df_metrics.Clone7$log10_nb_fragments, df_metrics.Clone7$TSS.enrichment, n = 5000)
df_metrics.Clone9$density <- get_density(df_metrics.Clone9$log10_nb_fragments, df_metrics.Clone9$TSS.enrichment, n = 5000)
df_metrics.Clone11$density <- get_density(df_metrics.Clone11$log10_nb_fragments, df_metrics.Clone11$TSS.enrichment, n = 5000)
df_metrics.Clone13$density <- get_density(df_metrics.Clone13$log10_nb_fragments, df_metrics.Clone13$TSS.enrichment, n = 5000)

#Plot TSS Enrichment vs log10(number of fragments)
p1 <- ggplot(df_metrics.Clone2) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  ggtitle("Clone 2 Fragments")+
  scale_color_viridis() +
  theme_bw()

p2 <- ggplot(df_metrics.Clone3) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  ggtitle("Clone 3 Fragments")+
  scale_color_viridis() +
  theme_bw()

p3 <- ggplot(df_metrics.Clone4) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  ggtitle("Clone 4 Fragments")+
  scale_color_viridis() +
  theme_bw()

p4 <- ggplot(df_metrics.Clone5) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  ggtitle("Clone 5 Fragments")+
  scale_color_viridis() +
  theme_bw()

p5 <- ggplot(df_metrics.Clone7) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  ggtitle("Clone 7 Fragments")+
  scale_color_viridis() +
  theme_bw()

p6 <- ggplot(df_metrics.Clone9) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  ggtitle("Clone 9 Fragments")+
  scale_color_viridis() +
  theme_bw()

p7 <- ggplot(df_metrics.Clone11) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  ggtitle("Clone 11 Fragments")+
  scale_color_viridis() +
  theme_bw()

p8 <- ggplot(df_metrics.Clone13) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  ggtitle("Clone 13 Fragments")+
  scale_color_viridis() +
  theme_bw()

p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8
ggsave("plots/Signac_TSSVsFrags.png", width = 12, height = 7)

#Calculate Nucleosome Signal Score
sample_atac.Clone2$nucleosome_group <- ifelse(sample_atac.Clone2$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_atac.Clone2, group.by = 'nucleosome_group')

sample_atac.Clone3$nucleosome_group <- ifelse(sample_atac.Clone3$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_atac.Clone3, group.by = 'nucleosome_group')

sample_atac.Clone4$nucleosome_group <- ifelse(sample_atac.Clone4$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_atac.Clone4, group.by = 'nucleosome_group')

sample_atac.Clone5$nucleosome_group <- ifelse(sample_atac.Clone5$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_atac.Clone5, group.by = 'nucleosome_group')

sample_atac.Clone7$nucleosome_group <- ifelse(sample_atac.Clone7$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_atac.Clone7, group.by = 'nucleosome_group')

sample_atac.Clone9$nucleosome_group <- ifelse(sample_atac.Clone9$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_atac.Clone9, group.by = 'nucleosome_group')

sample_atac.Clone11$nucleosome_group <- ifelse(sample_atac.Clone11$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_atac.Clone11, group.by = 'nucleosome_group')

sample_atac.Clone13$nucleosome_group <- ifelse(sample_atac.Clone13$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = sample_atac.Clone13, group.by = 'nucleosome_group')


p1 <- VlnPlot(object = sample_atac.Clone2,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)

p2 <- VlnPlot(object = sample_atac.Clone3,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)

p3 <- VlnPlot(object = sample_atac.Clone4,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)

p4 <- VlnPlot(object = sample_atac.Clone5,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)

p5 <- VlnPlot(object = sample_atac.Clone7,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)

p6 <- VlnPlot(object = sample_atac.Clone9,
              features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
              pt.size = 0.1, alpha = 0.5,
              ncol = 3)

p7 <- VlnPlot(object = sample_atac.Clone11,
              features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
              pt.size = 0.1, alpha = 0.5,
              ncol = 3)

p8 <- VlnPlot(object = sample_atac.Clone13,
              features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
              pt.size = 0.1, alpha = 0.5,
              ncol = 3)


p1 / p2 / p3 / p4 / p5 / p6 / p7 / p8
ggsave("plots/vln_QC_metrics.png", width = 7, height = 48)

## Quality Control of Cells from ArchR Project ----
#POT
#Calculate intersection of cellranger cells in peak matrix and ArchR cells (post QC)
#Currently CellID/colnames dont match so lets change cellranger ones to match ArchR
colnames(sample_atac.Clone2) == sample_atac.Clone2$CellID #FALSE
colnames(sample_atac.Clone2) <- sample_atac.Clone2$CellID
counts <- sample_atac.Clone2@assays$peaks$counts
intersect_cellranger_ArchR <- intersect(archr_cells, colnames(counts)) #ArchR 55618 cells 9265 intersecting cells
counts <- counts[ ,intersect_cellranger_ArchR] 

#Extract metadata for Arch QC'd cells
qc_meta <- sample_atac.Clone2@meta.data #12011 cells
qc_meta <- qc_meta %>% dplyr::filter(CellID %in% intersect_cellranger_ArchR) #9265 cells

#Create Seurat object with ArchR QC'd cells
sample_atac_qc.Clone2 <- CreateSeuratObject(counts = counts, meta.data = qc_meta, assay = "peaks")

rm(counts, intersect_cellranger_ArchR, qc_meta)

#Clone 3
colnames(sample_atac.Clone3) == sample_atac.Clone3$CellID #FALSE
colnames(sample_atac.Clone3) <- sample_atac.Clone3$CellID
counts <- sample_atac.Clone3@assays$peaks$counts
intersect_cellranger_ArchR <- intersect(archr_cells, colnames(counts)) #ArchR 5078 cells 55618 intersecting cells
counts <- counts[ ,intersect_cellranger_ArchR] 

#Extract metadata for Arch QC'd cells
qc_meta <- sample_atac.Clone3@meta.data #6257 cells
qc_meta <- qc_meta %>% dplyr::filter(CellID %in% intersect_cellranger_ArchR) #5078 cells

#Create Seurat object with ArchR QC'd cells
sample_atac_qc.Clone3 <- CreateSeuratObject(counts = counts, meta.data = qc_meta, assay = "peaks")

rm(counts, intersect_cellranger_ArchR, qc_meta)

#Clone 4
colnames(sample_atac.Clone4) == sample_atac.Clone4$CellID #FALSE
colnames(sample_atac.Clone4) <- sample_atac.Clone4$CellID
counts <- sample_atac.Clone4@assays$peaks$counts
intersect_cellranger_ArchR <- intersect(archr_cells, colnames(counts)) #ArchR 55618 cells 9265 intersecting cells
counts <- counts[ ,intersect_cellranger_ArchR] 

#Extract metadata for Arch QC'd cells
qc_meta <- sample_atac.Clone4@meta.data #12011 cells
qc_meta <- qc_meta %>% dplyr::filter(CellID %in% intersect_cellranger_ArchR) #9265 cells

#Create Seurat object with ArchR QC'd cells
sample_atac_qc.Clone4 <- CreateSeuratObject(counts = counts, meta.data = qc_meta, assay = "peaks")

rm(counts, intersect_cellranger_ArchR, qc_meta)

#Clone 5
colnames(sample_atac.Clone5) == sample_atac.Clone5$CellID #FALSE
colnames(sample_atac.Clone5) <- sample_atac.Clone5$CellID
counts <- sample_atac.Clone5@assays$peaks$counts
intersect_cellranger_ArchR <- intersect(archr_cells, colnames(counts)) #ArchR 55618 cells 9265 intersecting cells
counts <- counts[ ,intersect_cellranger_ArchR] 

#Extract metadata for Arch QC'd cells
qc_meta <- sample_atac.Clone5@meta.data #12011 cells
qc_meta <- qc_meta %>% dplyr::filter(CellID %in% intersect_cellranger_ArchR) #9265 cells

#Create Seurat object with ArchR QC'd cells
sample_atac_qc.Clone5 <- CreateSeuratObject(counts = counts, meta.data = qc_meta, assay = "peaks")

rm(counts, intersect_cellranger_ArchR, qc_meta)

#Clone 7
colnames(sample_atac.Clone7) == sample_atac.Clone7$CellID #FALSE
colnames(sample_atac.Clone7) <- sample_atac.Clone7$CellID
counts <- sample_atac.Clone7@assays$peaks$counts
intersect_cellranger_ArchR <- intersect(archr_cells, colnames(counts)) #ArchR 55618 cells 9265 intersecting cells
counts <- counts[ ,intersect_cellranger_ArchR] 

#Extract metadata for Arch QC'd cells
qc_meta <- sample_atac.Clone7@meta.data #12011 cells
qc_meta <- qc_meta %>% dplyr::filter(CellID %in% intersect_cellranger_ArchR) #9265 cells

#Create Seurat object with ArchR QC'd cells
sample_atac_qc.Clone7 <- CreateSeuratObject(counts = counts, meta.data = qc_meta, assay = "peaks")

rm(counts, intersect_cellranger_ArchR, qc_meta)

#Clone 9
colnames(sample_atac.Clone9) == sample_atac.Clone9$CellID #FALSE
colnames(sample_atac.Clone9) <- sample_atac.Clone9$CellID
counts <- sample_atac.Clone9@assays$peaks$counts
intersect_cellranger_ArchR <- intersect(archr_cells, colnames(counts)) #ArchR 55618 cells 9265 intersecting cells
counts <- counts[ ,intersect_cellranger_ArchR] 

#Extract metadata for Arch QC'd cells
qc_meta <- sample_atac.Clone9@meta.data #12011 cells
qc_meta <- qc_meta %>% dplyr::filter(CellID %in% intersect_cellranger_ArchR) #9265 cells

#Create Seurat object with ArchR QC'd cells
sample_atac_qc.Clone9 <- CreateSeuratObject(counts = counts, meta.data = qc_meta, assay = "peaks")

rm(counts, intersect_cellranger_ArchR, qc_meta)

#Clone 11
colnames(sample_atac.Clone11) == sample_atac.Clone11$CellID #FALSE
colnames(sample_atac.Clone11) <- sample_atac.Clone11$CellID
counts <- sample_atac.Clone11@assays$peaks$counts
intersect_cellranger_ArchR <- intersect(archr_cells, colnames(counts)) #ArchR 55618 cells 9265 intersecting cells
counts <- counts[ ,intersect_cellranger_ArchR] 

#Extract metadata for Arch QC'd cells
qc_meta <- sample_atac.Clone11@meta.data #12011 cells
qc_meta <- qc_meta %>% dplyr::filter(CellID %in% intersect_cellranger_ArchR) #9265 cells

#Create Seurat object with ArchR QC'd cells
sample_atac_qc.Clone11 <- CreateSeuratObject(counts = counts, meta.data = qc_meta, assay = "peaks")

rm(counts, intersect_cellranger_ArchR, qc_meta)

#Clone 13
colnames(sample_atac.Clone13) == sample_atac.Clone13$CellID #FALSE
colnames(sample_atac.Clone13) <- sample_atac.Clone13$CellID
counts <- sample_atac.Clone13@assays$peaks$counts
intersect_cellranger_ArchR <- intersect(archr_cells, colnames(counts)) #ArchR 55618 cells 9265 intersecting cells
counts <- counts[ ,intersect_cellranger_ArchR] 

#Extract metadata for Arch QC'd cells
qc_meta <- sample_atac.Clone13@meta.data #12011 cells
qc_meta <- qc_meta %>% dplyr::filter(CellID %in% intersect_cellranger_ArchR) #9265 cells

#Create Seurat object with ArchR QC'd cells
sample_atac_qc.Clone13 <- CreateSeuratObject(counts = counts, meta.data = qc_meta, assay = "peaks")

rm(counts, intersect_cellranger_ArchR, qc_meta)


#Visualise these metrics for each sample - POT
df_metrics <- sample_atac_qc.Clone2@meta.data[,c("nb_fragments", "TSS.enrichment")] 
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) 

p1 <- ggplot(df_metrics) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  scale_color_viridis() +
  theme_bw()

p2 <- VlnPlot(object = sample_atac_qc.Clone2,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)
p1 + p2
ggsave("plots/Clone2_signac_qc.png", width = 18, height = 7)

saveRDS(sample_atac.Clone2, "datafiles/Clone2_signac_preQC.RDS")
saveRDS(sample_atac_qc.Clone2, "datafiles/Clone2_signac_postQC.RDS")

#Clone 3
df_metrics <- sample_atac_qc.Clone3@meta.data[,c("nb_fragments", "TSS.enrichment")] 
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) 

p1 <- ggplot(df_metrics) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  scale_color_viridis() +
  theme_bw()

p2 <- VlnPlot(object = sample_atac_qc.Clone3,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)
p1 + p2
ggsave("plots/Clone3_signac_qc.png", width = 18, height = 7)

saveRDS(sample_atac.Clone3, "datafiles/Clone3_signac_preQC.RDS")
saveRDS(sample_atac_qc.Clone3, "datafiles/Clone3_signac_postQC.RDS")

#Clone 4
df_metrics <- sample_atac_qc.Clone4@meta.data[,c("nb_fragments", "TSS.enrichment")] 
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) 

p1 <- ggplot(df_metrics) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  scale_color_viridis() +
  theme_bw()

p2 <- VlnPlot(object = sample_atac_qc.Clone4,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)
p1 + p2
ggsave("plots/Clone4_signac_qc.png", width = 18, height = 7)

saveRDS(sample_atac_qc.Clone4, "datafiles/Clone4_signac_preQC.RDS")
saveRDS(sample_atac_qc.Clone4, "datafiles/Clone4_signac_postQC.RDS")

#Clone 5 ######
df_metrics <- sample_atac_qc.Clone5@meta.data[,c("nb_fragments", "TSS.enrichment")] 
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) 

p1 <- ggplot(df_metrics) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  scale_color_viridis() +
  theme_bw()

p2 <- VlnPlot(object = sample_atac_qc.Clone5,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)
p1 + p2
ggsave("plots/Clone5_signac_qc.png", width = 18, height = 7)

saveRDS(sample_atac_qc.Clone5, "datafiles/Clone5_signac_preQC.RDS")
saveRDS(sample_atac_qc.Clone5, "datafiles/Clone5_signac_postQC.RDS")

#Clone 7
df_metrics <- sample_atac_qc.Clone7@meta.data[,c("nb_fragments", "TSS.enrichment")] 
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) 

p1 <- ggplot(df_metrics) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  scale_color_viridis() +
  theme_bw()

p2 <- VlnPlot(object = sample_atac_qc.Clone7,
        features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
        pt.size = 0.1, alpha = 0.5,
        ncol = 3)
p1 + p2
ggsave("plots/Clone7_signac_qc.png", width = 18, height = 7)

saveRDS(sample_atac_qc.Clone7, "datafiles/Clone7_signac_preQC.RDS")
saveRDS(sample_atac_qc.Clone7, "datafiles/Clone7_signac_postQC.RDS")

#Clone 9
df_metrics <- sample_atac_qc.Clone9@meta.data[,c("nb_fragments", "TSS.enrichment")] 
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) 

p1 <- ggplot(df_metrics) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  scale_color_viridis() +
  theme_bw()

p2 <- VlnPlot(object = sample_atac_qc.Clone9,
              features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
              pt.size = 0.1, alpha = 0.5,
              ncol = 3)
p1 + p2
ggsave("plots/Clone9_signac_qc.png", width = 18, height = 7)

saveRDS(sample_atac_qc.Clone9, "datafiles/Clone9_signac_preQC.RDS")
saveRDS(sample_atac_qc.Clone9, "datafiles/Clone9_signac_postQC.RDS")

#Clone 11
df_metrics <- sample_atac_qc.Clone11@meta.data[,c("nb_fragments", "TSS.enrichment")] 
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) 

p1 <- ggplot(df_metrics) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  scale_color_viridis() +
  theme_bw()

p2 <- VlnPlot(object = sample_atac_qc.Clone11,
              features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
              pt.size = 0.1, alpha = 0.5,
              ncol = 3)
p1 + p2
ggsave("plots/Clone11_signac_qc.png", width = 18, height = 7)

saveRDS(sample_atac_qc.Clone11, "datafiles/Clone11_signac_preQC.RDS")
saveRDS(sample_atac_qc.Clone11, "datafiles/Clone11_signac_postQC.RDS")

#Clone 13
df_metrics <- sample_atac_qc.Clone13@meta.data[,c("nb_fragments", "TSS.enrichment")] 
df_metrics$log10_nb_fragments <- log10(df_metrics$nb_fragments)
df_metrics$density <- get_density(df_metrics$log10_nb_fragments, df_metrics$TSS.enrichment, n = 5000) 

p1 <- ggplot(df_metrics) + 
  geom_point(aes(log10_nb_fragments, TSS.enrichment, color = density)) + 
  scale_color_viridis() +
  theme_bw()

p2 <- VlnPlot(object = sample_atac_qc.Clone13,
              features = c("nb_fragments", "TSS.enrichment", "nucleosome_signal"),
              pt.size = 0.1, alpha = 0.5,
              ncol = 3)
p1 + p2
ggsave("plots/Clone13_signac_qc.png", width = 18, height = 7)

saveRDS(sample_atac_qc.Clone13, "datafiles/Clone13_signac_preQC.RDS")
saveRDS(sample_atac_qc.Clone13, "datafiles/Clone13_signac_postQC.RDS")
