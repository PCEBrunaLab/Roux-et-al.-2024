---
title: "shallow Whole Genome Sequencing Analysis Pipeline"
output: 
  pdf_document: default
  html_document: default
---
# Whole Genome Sequencing Analysis

## Overview

This pipeline is designed for the processing and analysis of shallow Whole Genome Sequencing (WGS) data. It consists of several Bash scripts, each handling a specific part of the analysis process, including quality control, sequence alignment, sorting and indexing, and duplicate removal. These scripts are tailored to run on a High-Performance Computing (HPC) environment utilizing the SLURM Workload Manager.The R script performs quantitative DNA sequencing analysis for chromosomal aberrations on multiple WGS samples using the QDNAseq package. 

## Requirements 
FastQC,
MultiQC,
BWA,
SAMtools,
Picard Tools,
QDNAseq package in R.

## Raw Data Characteristics

### Data Format

- **Format:** The pipeline processes data in FASTQ format, which is the standard output from most next-generation sequencing (NGS) platforms.

## Pipeline Structure

The pipeline is structured as a series of individual SLURM jobs and R script, allowing for efficient resource management and parallel processing in a cluster environment. Each script corresponds to a specific stage in the sWGS data analysis workflow:

1. **Quality Control (QC.sh):** 
   - Performs quality checks on raw sequence data using FastQC and aggregates the reports using MultiQC.

2. **BWA Alignment (BWA_Alignment.sh):**
   - Aligns sequencing reads to a reference genome using BWA and processes the alignments using SAMtools.

3. **Samtools Sorting and Indexing (Samtools_Sort_Index.sh):**
   - Sorts and indexes the aligned reads using SAMtools.

4. **Remove Duplicates (Remove_Duplicates_Picard.sh):**
   - Removes duplicate reads from the aligned data using Picard Tools.

5. **QDNAseq Analysis (QDNAseq_Analysis.R):**
   - Analyzes processed data to identify chromosomal aberrations using the QDNAseq package.

## Running the Pipeline

Each script within the pipeline is designed to be submitted as an individual job to the SLURM scheduler. It's crucial to run each script in the correct order, as the output from one step serves as the input for the subsequent step. To ensure this sequence, you can set job dependencies when submitting your jobs to SLURM. This means that a job will only start once the previous job has successfully completed.
Make sure that all dependencies, such as software modules or libraries, are loaded at the start of each script. Additionally, double-check that all file paths are correctly set, directing each script to the appropriate input files and specifying where output files should be saved. By carefully managing the input and output of each stage, the pipeline will process the data through each step in the pipeline smoothly and efficiently.

