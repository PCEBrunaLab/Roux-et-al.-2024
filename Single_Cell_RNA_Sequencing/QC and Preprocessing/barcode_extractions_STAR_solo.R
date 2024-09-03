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

print("start")
start <- Sys.time()

parser <- OptionParser()

parser <- add_option(parser, c("-i", "--sample"), type="character",
                     help="Sample input directory")


opt <- parse_args(parser)

#Install easypar package from GitHub
devtools::install_github("caravagnalab/easypar", force = TRUE)

#Write function align_barcodes which defines two reads for each sample, an output folder to place bamfiles and a reference genome to compare to.
align_barcodes_scRNA <- function(sample) {
  library(glue)
  R1 <- paste0(sample, "_R1_001.fastq.gz")
  R2 <- paste0(sample, "_R2_001.fastq.gz")
  output <- "outputs/Cellecta17/"
  genome <- "GenomeDir"
  command <- glue("STAR --runThreadN 48 --genomeDir {genome} --readFilesIn {R2} {R1} --readFilesCommand zcat --outFileNamePrefix {output} --soloType Droplet --soloFeatures GeneFull --soloCBwhitelist None --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 60 --outFilterScoreMin 30 --outFilterMismatchNmax  6 --outSAMmode Full --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM
")
  system(command)
}


#Define the files in which you can find the names for the sequencing samples
samples_sc_rna_seq <- list.files("fastq",pattern = "*.fastq.gz", recursive = T, full.names = T)
print(samples_sc_rna_seq)

#Substitute the names generated above without the file extras
samples_sc_rna_seq <- gsub("_R[12]_001.fastq.gz", "", samples_sc_rna_seq) %>% unique()
print(samples_sc_rna_seq)

#Show how many samples
sample_length <- samples_sc_rna_seq %>% length()
print(sample_length)

#Use EasyPar package to run align_barcodes function, using samples defined above
easypar::run_SLURM(FUN = align_barcodes_scRNA, PARAMS = data.frame(samples = samples_sc_rna_seq), N_simultaneous_jobs = 50, per_task = 5)

print("checkpoint")

#Apply the naming files to the bam files generated from align_barcodes
ll <- lapply(samples_sc_rna_seq, align_barcodes_scRNA)
