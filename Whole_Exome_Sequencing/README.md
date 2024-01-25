# Exome Sequencing Analysis Pipeline

## Overview

This document outlines a series of bash scripts used for processing and analyzing Exome Sequencing data. These scripts are designed for a High-Performance Computing (HPC) environment using the SLURM Workload Manager.

## Requirements
  - FastQC
  - MultiQC
  - BWA
  - SAMtools
  - Picard
  - GATK
  - MuTect2
  - BCFtools
  - ANNOVAR 
  - maftools
  

## Raw Data Characteristics

### Data Format
- **Format:** The pipeline processes data in FASTQ format, which is the standard output from most next-generation sequencing (NGS) platforms.

## Pipeline Structure

The pipeline is structured as a series of individual SLURM jobs, allowing for efficient resource management and parallel processing in a cluster environment. Each script corresponds to a specific stage in the WES data analysis workflow:

1. **Quality Control (QC.sh)** 
   - Performs quality checks on raw sequence data using FastQC and aggregates the reports using MultiQC.

2. **GATK Resource Bundle Download (GATKBundle_2.sh)**
  - Downloads the Genome Analysis Toolkit (GATK) resource bundle for human genome reference hg38.

3. **Genome Indexing with BWA (GenomeIndex_3.sh)**
  - Aligns sequencing reads to a reference genome using BWA and processes the alignments using SAMtools.

4. **Aligning Sequences with BWA (BWA_MAPPING_4.sh)**
 - Aligns sequencing reads to the reference genome using BWA and processes the output with SAMtools.

5. **Sorting and Indexing BAM Files with Samtools (Alignment Post-Processing_ Sorting_5.sh)**
  - Sorts and indexes the aligned BAM files using SAMtools for efficient data access and downstream processing.

6. **Samtools Coverage and Mapping Statistics (Alignment Post-Processing_ Coverage_6.sh)**
  - Computes coverage and mapping statistics for the aligned reads, providing insights into the quality of the alignment using SAMtools.

7. **Remove Duplicate Reads with Picard (Picard_7.sh)**
  - Removes duplicate reads from the aligned BAM files, a crucial step in ensuring the quality of variant calls, using Picard Tools.

8. **Add Read Groups with Picard (Add_ReadGroups_Picard.sh)**
  - Adds read group information to the BAM files, which is necessary for many downstream analyses, using Picard's AddOrReplaceReadGroups tool.

9. **Samtools Indexing of Recalibrated BAM Files (Samtools_index_RG_9.sh)**
  - Indexes recalibrated BAM files using SAMtools, necessary for efficient data retrieval in subsequent analyses.

10. **Apply Base Quality Score Recalibration with GATK (BQSR_10.sh)**
  - Applies base quality score recalibration using GATK, improving the accuracy of variant calls.

11. **Variant Calling with MuTect2 (Mutect2_11.sh)**
  - Performs variant calling on the processed BAM files using GATK's MuTect2, identifying potential somatic mutations.

12. **Filtering MuTect2 Calls (FilterMutectCalls_12.sh)**
  - Filters the output of MuTect2 using GATK's FilterMutectCalls, refining the list of variant calls.

13. **Germline Variant Filtering for Strelka (GermlineFilter_13.sh)**
  - Filters out known germline variants from variant calls using BCFtools, focusing on somatic mutations.

14. **Selecting Variants from Filtered VCF Files (Select_Variants_GATK.sh)**
  - Extracts variants from filtered VCF files 

15. **Applying Additional Filtering Criteria (Mutect_FilterCriteria_15.sh)**
  - Applies additional filtering criteria to the VCF files to refine variant calls

16. **Targeted Variant Selection with SureSelect (TTARGET_BED_SureSelect_16.sh)**
  - Selects variants based on a targeted BED file, focusing analysis on specific genomic regions 

17. **ANNOVAR Database Download (ANNOVAR_dbDownload_17.sh)**
  - Downloads essential databases for variant annotation using ANNOVAR

18. **Preparing ANNOVAR Input (ANNOVAR_inputformat_18.sh)**
  - Converts VCF files to ANNOVAR input format, preparing them for comprehensive annotation.

19. **Variant Annotation with ANNOVAR (ANNOVAR_annotation_19.sh)**
  -  Performs functional annotation of variants using ANNOVAR

20. **Conversion of ANNOVAR Output to MAF Format (ANNOVAR_to_MAF_20.sh)**
  - Converts the annotated variant data from ANNOVAR's format to the Mutation Annotation Format (MAF), facilitating        further analysis and visualization. 
  
21. **calculate variant allele frequencies (VAFs) (VAF_calculation_21.sh)**
  - The VAF results for each input file to use for downstream analysis. 


## Running the Pipeline

Each script within the pipeline is designed to be submitted as an individual job to the SLURM scheduler. It's crucial to run each script in the correct order, as the output from one step serves as the input for the subsequent step. To ensure this sequence, you can set job dependencies when submitting your jobs to SLURM. This means that a job will only start once the previous job has successfully completed.
Make sure that all dependencies, such as software modules or libraries, are loaded at the start of each script. Additionally, double-check that all file paths are correctly set, directing each script to the appropriate input files and specifying where output files should be saved. By carefully managing the input and output of each stage, the pipeline will process the data through each step in the pipeline smoothly and efficiently.

## Using Mutation Annotation Format (MAF) for Analysis

At the end of the Exome Sequencing Analysis Pipeline, the generated MAF files can be utilized for comprehensive mutation analysis. The Mutation Annotation Format is widely used for summarizing, analyzing, annotating, and visualizing genomic mutations.

## Preparation of Input Files for Heatmap Visualization
Utilize Excel to prepare and format the input files for heatmap visualization using the VAF_HeatmapPlot.R script. 

