#1.snATAC-seq Preprocessing

#Work outside of REnv in R version 4.4 to use up to date packages for Signac processing

#Set Working Directory and Arch R parameters
setwd("~/OneDrive/PhD/Thesis/Thesis Plots/Results Chapter 3/Sub Chapter 1/")

#Install extra packages
required.packages <- c("dplyr", "tidyr", "Signac", "Seurat", "GenomeInfoDb", "EnsDb.Hsapiens.v86",
                       "BSgenome.Hsapiens.UCSC.hg38", "patchwork", "MASS", "viridis", "hdf5r")

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
library(hdf5r)
set.seed(76)

## Peak Calling with MACS3 POT ----
#Generate combined peak set across all samples
#Retrieve gene annotations for TSS enrichment
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevels(annotation) <- paste0("chr", seqlevels(annotation)) 

#Call Peaks using MACS2 per sample
macs2_path = "/Users/sianhamer/MACS3env/bin/macs3"

#Read in all necessary cellranger output files per clone
read_all_clone_data_rds <- function(clone_names = c("Clone 2", "Clone 3", "Clone 4", "Clone 5",
                                                    "Clone 7", "Clone 9", "Clone 11", "Clone 13"),
                                    base_dir = "datafiles") {
  data_list <- lapply(clone_names, function(clone) {
    # Remove spaces from the clone name to match the file naming convention
    file_name <- paste0(gsub(" ", " ", clone), "_signac_postQC.RDS")
    file_path <- file.path(base_dir, file_name)
    readRDS(file_path)
  })
  names(data_list) <- clone_names
  return(data_list)
}

all_clone_data <- read_all_clone_data_rds()


#Loop through each clone file and conduct individual peak calling
# Loop over each Seurat object in the list for peak calling
for (clones in names(all_clone_data)) {
  
  #Extract each object
  seurat_obj <- all_clone_data[[clones]]
  
  colnames(seurat_obj) <- paste0(clones, "_", colnames(seurat_obj))
  colnames(seurat_obj) <- gsub(" ", "_", colnames(seurat_obj))
  
  # Read 10X filtered feature matrix and metadata
  counts <- Read10X_h5(paste0("datafiles/cellranger outputs/", clones, "/filtered_feature_bc_matrix.h5"))[[2]]
  colnames(counts) <- paste0(clones, "_", colnames(counts))
  colnames(counts) <- gsub(" ", "_", colnames(counts))
  metadata <- read.csv(file = paste0("datafiles/cellranger outputs/", clones, "/per_barcode_metrics.csv"),
                       header = TRUE, row.names = 1)
  rownames(metadata) <- paste0(clones, "_", rownames(metadata))
  rownames(metadata) <- gsub(" ", "_", rownames(metadata))
  fragpath <- paste0("datafiles/cellranger outputs/", clones, "/atac_fragments.tsv.gz")
  
  # Create Seurat object
  tmp <- CreateSeuratObject(counts = counts, assay = "ATAC", meta.data = metadata)
  
  # Extract cells for filtering
  cells_to_filter <- colnames(tmp)
  
  # Find intersecting cells
  intersect_cells <- intersect(colnames(seurat_obj), colnames(tmp))
  counts <- tmp@assays$ATAC$counts[, intersect_cells]
  
  # Remove sample prefix in colnames for compatibility with fragment file
  clones <- gsub(" ", "_", clones)
  colnames(seurat_obj) <- sub("(Clone)_([0-9]*)_([ACGT]*-1)", "\\3", colnames(seurat_obj))
  colnames(counts) <- sub("(Clone)_([0-9]*)_([ACGT]*-1)", "\\3", colnames(counts))
  
  # Create Chromatin Assay
  seurat_obj[["ATAC"]] <- CreateChromatinAssay(counts = counts, sep = c(":", "-"),
                                                fragments = fragpath, annotation = annotation)
  DefaultAssay(seurat_obj) <- "ATAC"
  
  # Call peaks using MACS2
  peaks <- CallPeaks(object = seurat_obj, macs2.path = macs2_path)
  
  # Filter peaks to retain only standard chromosomes and remove blacklist regions
  peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  chroms_to_keep <- which(seqnames(peaks) %in% paste0("chr", 1:22))
  peaks <- peaks[chroms_to_keep, ]
  
  # Create peak-cell counts matrix
  macs2_counts <- FeatureMatrix(fragments = Fragments(seurat_obj), features = peaks, cells = colnames(seurat_obj))
  
  # Create chromatin assay
  sample_atac2 <- CreateChromatinAssay(macs2_counts, fragments = fragpath, annotation=annotation)
  
  # Create new Seurat object with updated chromatin assay
  sample_atac_macs2 <- CreateSeuratObject(sample_atac2, assay = "ATAC", meta.data = seurat_obj@meta.data)
  Idents(sample_atac_macs2) <- sample_atac_macs2$Sample
  
  # Save results
  saveRDS(peaks, paste0("datafiles/", clones, "_peaks_granges_per_subtype.RDS"))
  saveRDS(macs2_counts, paste0("datafiles/", clones, "_peak_cell_counts_matrix.RDS"))
  saveRDS(sample_atac_macs2, paste0("datafiles/", clones, "_signac_MACS2.RDS"))
}

