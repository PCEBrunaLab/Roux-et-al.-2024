#install BiocManager to install further packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

#install packages needed for downstream analysis
required.packages <- c("dplyr", "glue", "Rsubread", "ggplot2", "devtools", "optparse")

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

#Install easypar package from GitHub
devtools::install_github("caravagnalab/easypar", force = TRUE)

#Write function align_barcodes which defines two reads for each sample, an output folder to place bamfiles and a reference genome to compare to.
align_barcodes <- function(sample) {
  library(glue)
  R1 <- paste0(sample, "_R1_001.fastq.gz")
  R2 <- paste0(sample, "_R2_001.fastq.gz")
  output <- paste0("bam_files/",gsub("^.*/", "", sample), ".bam")
  genome <- "bwa_index/barcodes_full.fasta"
  command <- glue("bwa mem -T 30 -t 12 {genome} {R1} {R2} | samtools sort -@12 -o {output} - ")
  system(command)
} 

#List files within working directory to define number of samples you are working with, output this number
samples <- list.files("fastq/NB/", recursive = T, full.names = T)
print(samples)
samples <- gsub("J63_[0-9]*/", "", samples) %>% unique()
print(samples)
sample_length <- samples %>% length()
print(sample_length)

#Use EasyPar package to run align_barcodes function, using samples defined above
easypar::run_SLURM(FUN = align_barcodes, PARAMS = data.frame(samples = samples), N_simultaneous_jobs = 50, per_task = 5)

#Define the files in which you can find the names for the sequencing samples
samples_sc_rna_seq <- list.files("fastq/NB/",pattern = "*.fastq.gz", recursive = T, full.names = T)

#Substitute the names generated above without the file extras
samples_sc_rna_seq <- gsub("_R[12]_001.fastq.gz", "", samples_sc_rna_seq) %>% unique()

#Apply the naming files to the bam files generated from align_barcodes
ll <- lapply(samples_sc_rna_seq, align_barcodes)


