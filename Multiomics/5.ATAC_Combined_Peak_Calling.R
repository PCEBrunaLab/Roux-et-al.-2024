#Combined Peak Calling

#Work outside of REnv in R version 4.4 to use up to date packages for Signac processing

#Set Working Directory and Arch R parameters
setwd("~/OneDrive/PhD/Thesis/Thesis Plots/Results Chapter 3/Sub Chapter 1/")

## Load in Packages ----
required.packages <- c("dplyr", "tidyr", "Signac", "Seurat", "future", "GenomeInfoDb", "EnsDb.Hsapiens.v86",
                       "BSgenome.Hsapiens.UCSC.hg38", 'GenomicRanges', "patchwork", "MASS", "viridis")

if (!all(required.packages %in% installed.packages()[,"Package"])){
  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]], force = TRUE)
}

library(dplyr)
library(tidyr)
library(Signac)
library(Seurat)
library(future)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)

#Load in datafiles
#Read in all necessary cellranger output files per clone
clone_names <- c("Clone_2", "Clone_3", "Clone_4", "Clone_5",
                 "Clone_7", "Clone_9", "Clone_11", "Clone_13")
base_dir <- "datafiles"

# Loop through each clone and load the corresponding RDS file
for (clone in clone_names) {
  # Construct the file path
  file_name <- paste0(clone, "_peaks_granges_per_subtype.RDS")
  file_path <- file.path(base_dir, file_name)
  
  # Check if file exists before loading
  gr_obj <- readRDS(file_path)
  
  # Assign the object to the global environment with the clone name
  assign(clone, gr_obj, envir = .GlobalEnv)
}

#Retrieve Gene Annotation 
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) 
seqlevels(annotation) <- paste0("chr", seqlevels(annotation)) 

## Generate common peak set across samples ----
#Create reduce combined peaks
combined.peaks.reduced <- GenomicRanges::reduce(x = c(Clone_2, Clone_3, Clone_4, Clone_5, Clone_7,
                                       Clone_9, Clone_11, Clone_13)) 

#Create disjoined combined peaks
combined.peaks.disjoin <- GenomicRanges::disjoin(x = c(Clone_2, Clone_3, Clone_4, Clone_5, Clone_7,
                                        Clone_9, Clone_11, Clone_13)) 

#Fitler out bad peaks based on length for reduced
peakwidths <- width(combined.peaks.reduced) 
min(peakwidths) #200 minimun
max(peakwidths) #4482 maximum

#NB. no need for filtering as our datarang is already between 20 and 10,000
#combined.peaks.reduced <- combined.peaks.reduced[peakwidths  < 10000 & peakwidths > 20]
combined.peaks.reduced # 252084 peaks in total

saveRDS(combined.peaks.reduced, "datafiles/combined_peaks_reduced.RDS")

#Repeat for combined peaks with disjoin
#Filter out bad peaks based on length for reduced
peakwidths <- width(combined.peaks.disjoin) 
min(peakwidths) #1 minimun
max(peakwidths) #2735 maximum

#NB. this time we need to filter as we have peak less than 20
combined.peaks.disjoin <- combined.peaks.disjoin[peakwidths  < 10000 & peakwidths > 20]
combined.peaks.disjoin # 877979 peaks in total

saveRDS(combined.peaks.disjoin, "datafiles/combined_peaks_disjoin.RDS")

#NB. Signac merging samples tutorial uses reduce so lets proceed with this one for now for downstream analysis

#Convert Granges to BED
combined.peaks_df <- data.frame(chr=seqnames(combined.peaks.reduced),
                                start=start(combined.peaks.reduced),
                                end=end(combined.peaks.reduced))
combined.peaks_df$name <- paste(combined.peaks_df$chr, combined.peaks_df$start, combined.peaks_df$end, sep="_")
combined.peaks_df$score=1

write.table(combined.peaks_df, file="datafiles/combined_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


