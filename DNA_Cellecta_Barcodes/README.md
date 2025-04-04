# Cellecta Barcode Extraction and Processing from DNA

## Overview

This pipeline is designed for the extraction and processing of Cellecta cellular barcodes from DNA sequencing. It consists of several bash scripts which will generate custom reference for fastq alignment, algn fastqs to generate bam files for barcode extraction, and an R script which will summarise Cellecta barcode presence within samples. Following this, there is also an R analysis script which takes summary of Cellecta barcode information and analyses inline with "Dynamic Plasticity Systems Direct Early Adaptation to Treatment in Neuroblastoma" manuscript.

## Pipeline Structure

The pipeline is structured as a series of SLURM jobs, allowing for efficient resource management and parallel processing in a cluster environment. Followed by R analysis scripts to analysis downstream of processing:

1. **Custom Reference Generation from Barcode Whitelist (bwa_index_generation.R):** 
   - Generates genome reference for fastq alignment of DNA sequencing for Cellecta barcodes.

2. **Fastq Alignment and Bam Generation (bwa_aligment.sh & bwa_alignment.R):**
   - Aligns fastq to previously generated custom reference using BWA-mem and generates bam files using SAMtools.

3. **Summary of Cellecta barcode counts (FeatureCounts.R):**
   - Summarising of cellecta barcode information per sample using FeatureCounts.

4. **Analysis of Cellecta Barcodes (ClonalDynamics.R):**
   - Final script conducting analysis as per manuscript (Roux, Hamer & Shea, 2024) to identify clonal dynamics.

## Running the Pipeline

To run the pipeline, run R script (locally or on cluster) to conduct full analysis pipeline.

bwa_index_generation.R

Next, submit each SLURM script to the cluster using the `sbatch` command. Ensure that all scripts are correctly configured with the necessary file paths and resource requirements before submission.

sbatch bwa_aligment.sh 

Following SLURM scripts, run R scripts (locally or on cluster) to conduct full analysis pipeline.

FeatureCounts.R \
ClonalDynamics.R

