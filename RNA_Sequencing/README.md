
---
title: "Bulk RNA-seq Data Analysis Pipeline"
output: 
  pdf_document: default
  html_document: default
---

# Bulk RNA-seq Data Analysis Pipeline

## Overview

This document details the scripts for processing and analyzing Bulk RNA-seq data in a High-Performance Computing (HPC) environment. The pipeline includes steps for quality control, alignment, quantification, and differential expression analysis.

## Requirements

- FastQC
- STAR
- Samtools
- HTSeq
- R 

## Raw Data Characteristics

### Data Format

- **Format:** Processes FASTQ format data, typical for Bulk RNA-seq outputs.

## Pipeline Structure

These scripts are designed for a High-Performance Computing (HPC) environment using the SLURM Workload Manager. The pipeline is organized as a series of SLURM jobs and R scripts for efficient processing. Each script corresponds to a stage in the RNA-seq data analysis workflow:


1. **Quality Control (QC.sh)**
   - Uses FastQC for quality checks on raw sequence data.
   
2. **Genome Indexing with STAR (GenomeGenerate.sh)**
   - Creates an index of the reference genome using STAR.

3. **Read Alignment (STARMapping.sh)**
   - Aligns reads to a reference genome using STAR.

4. **Sorting and Indexing (BamIndex.sh)**
   - Utilizes Samtools for sorting and indexing aligned reads.

5. **Quantification with HTSeq (HTSeq.sh)**
   - Quantifies gene expression levels using HTSeq.

6. **Merging counts matrix (Merged_HTSeq.R)**
  - Merges multiple count files into one merged file with counts for all samples.

7. **Differential Expression Analysis (DEG_SK_N_SH&Organoid.R)**
   - Performs differential expression analysis using DESeq2 in R and also visualize  results.


## Running the Pipeline

Each script within the pipeline is designed to be submitted as an individual job to the SLURM scheduler. It's crucial to run each script in the correct order, as the output from one step serves as the input for the subsequent step. To ensure this sequence, you can set job dependencies when submitting your jobs to SLURM. This means that a job will only start once the previous job has successfully completed.
Make sure that all dependencies, such as software modules or libraries, are loaded at the start of each script. Additionally, double-check that all file paths are correctly set, directing each script to the appropriate input files and specifying where output files should be saved. By carefully managing the input and output of each stage, the pipeline will process the data through each step in the pipeline smoothly and efficiently.



