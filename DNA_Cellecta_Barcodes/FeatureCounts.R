#FeatureCounts.R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

#install packages needed for downstream analysis
required.packages <- c("dplyr", "glue", "Rsubread", "ggplot2", "devtools", "optparse", "Polychrome")
if (!all(required.packages %in% installed.packages()[,"Package"])){
  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]])
}

#Load packages into library
library(dplyr)
library(glue)
library(Rsubread)
library(ggplot2)
library(devtools)
library(optparse)
library(reshape2)
library(Polychrome)

#Set Working Directory
setwd("/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J63_DNA_sequencing/")

#Load in .tsv with all Cellecta barcode sequences
saf <- data.table::fread("bwa_index/saf_barcodes.tsv", data.table = F)

#List bam files for processing
cell_files <- list.files("bam_files", pattern = "*.bam", full.names = T)

#Feature counts
counts_sc <- featureCounts(cell_files,isPairedEnd = TRUE ,annot.ext = saf, nthreads= 2, verbose = FALSE, 
                           maxMOp = 16, minMQS = 30,                          
                           countChimericFragments = FALSE, allowMultiOverlap = FALSE, nonSplitOnly = TRUE,
                           requireBothEndsMapped = TRUE, countMultiMappingReads = FALSE, primaryOnly = TRUE)

#Remove unneccesary information
df_sc <- counts_sc$counts %>% reshape2::melt() 
colnames(df_sc) <- c("real_bc44", "sample_id", "N")
df_sc %>% head()

#Filter for barcodes which are present at least once
df_sc <- df_sc %>% filter(N > 0)

#Add percentage representation of barcode within each sample
df_sc <- df_sc %>% group_by(sample_id) %>% mutate(N_tot = sum(N) ) %>% ungroup %>% mutate(perc = N/N_tot)

#Save summarised information for cellecta barcodes
df_sc %>% data.table::fwrite("all_barcodes.tsv")