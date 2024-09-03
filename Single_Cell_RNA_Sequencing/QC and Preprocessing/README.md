# Single Cell RNA Sequencing Processing

## Overview

This pipeline is designed for the extraction and processing of Cellecta cellular barcodes from single-cell RNA sequencing. The first set fo scripts detail the extractions, error correction and processing of cellecta barcodes from fastq. This consists of a single R script which works through the analysis of these clones as per manuscript "Dynamic Plasticity Systems Direct Early Adaptation to Treatment in Neuroblastoma" manuscript.

Summary excel spreadsheet has been provided to indicate input and output files from each stage of analysis for each of use (workflow_summary.csv)

## Pipeline Structure

The pipeline is structured as a series of SLURM jobs, allowing for efficient resource management and parallel processing in a cluster environment. Followed by R analysis scripts to analysis downstream of processing:

1. **Count matrices generation with CellRanger multi (cellranger_multi.sh & config.csv):**
   - Analyse cell multiplexed samples from 10x, performs alignment, filteringm barcode counting and UMI counting.

Two different methods were used for Cellecta barcode extraction; first method 2A through 2A - 3 were used for cellecta barcode extraction for dataset comprising of SK-N-SH treatment with onyl Cisplatin and BRD4i for 48hrs, second method 2B was used for extended dataset in SK-N-SH cells comprising of EZH2i and BRD4i samples and extended treatment lengths up to 7days. One method was used for cellecta barcode extractions and not both.
2A. **Cellecta Barcode Processing from FASTQ (barcode_process.sh & process_sc_barcodes.py):** 
   - Looping over FASTQ files from Illumina sequencing and validate Cellecta barcodes.

2A-1. **Error Correction of Cellecta Barcodes (error_correct_bcs.sh & barcode_error_correct.R):**
   - Perform error correction on barcodes using hashing.

2A-2. **Extract Celelcta barcodes from bam files (extract_barcodes.sh & extract_barcodes.py):**
   - Extract the cellranger corrected barcodes from the alignment bam files - map to the observed/sequenced cellranger barcodes.

2A-3. **Compile list of uniquly observed barcodes (compress_barcodes.sh & compress_barcodes.py):**
   - Run over the cellbarcodes and compress them down to the unique observed barcodes.


2B. **Cellecta Barcode Processing from FASTQ using STARsolo (submit_barcode_extractions_STARsolo.sh & barcode_extractions_STARsolo.R):** 
   - Extracting Cellecta barcodes using whitelist of possible barcode sequences from fastq using STARsolo.


3. **Identify and drop empty droplets with EmptyDrops (emptydrops.sh & EmptyDrops.R):**
   - Takes in all the 10x experiment cell calls and normalises them together. EmptyDrops removes poor quality cells.

4. **Summarise Cell and CMO information (cmo_summary.sh & cmo_stats.R):**
   - Compute summarises for each sample of CMOs and cellIDs.

5. **Demultiplex samples and assign CMOs to pools (demux.sh & demux.R):**
   - Uses the joint distribution of counts over all CMOs to assign cells from pools to correct CMOs.

6. **Sample Normalisation and generation of single SCE object (normalise.sh & norm.R):**
   - Takes in RDS objects from each pool and normalises across samples to generation a single RDS object.

7. **Associate Cellecta Barcodes with SCE object (single_cell_cellecta_barcodes.R):**
   - Identify cell quality in terms of cellecta barcodes, then collate metadata for both sample varibles and cellecta barcodes for each cell.

8. **scRNA-seq analysis workflow using Seurat (analysis_seurat.R):**

    - Data preparation, QC and filtering
    - Normalisation, scaling and variable features using `SCTransform()`
    - Dimensionality reduction
    - Clustering and marker finding
    - AMT phenotyping



## Running the Pipeline

To run the pipeline, submit each SLURM script to the cluster using the `sbatch` command. Ensure that all scripts are correctly configured with the necessary file paths and resource requirements before submission. The first few scripts require a large amount of processing power and must be run on the cluster.

sbatch cellranger_multi.sh \
sbatch barcode_process.sh \
sbatch error_correct_bcs.sh \
sbatch extract_barcodes.sh \
sbatch compress_barcodes.sh \
sbatch emptydrops.sh \
sbatch cmo_summary.sh \
sbatch demux.sh \
sbatch normalise.sh 

Following SLURM scripts, run R scripts (locally or on cluster) to conduct full analysis pipeline.

single_cell_cellecta_barcodes.R \
analysis_seurat.R
