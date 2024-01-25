# Single Cell Clonal Dynamics Analysis

## Overview

This analysis is designed to explore and understand the clonal dynamics which occur with our Single Cell RNA-sequencing using Cellecta barcodes.
This consists of a single R script which works through the analysis of these clones as per manuscript "Dynamic Plasticity Systems Direct Early Adaptation to Treatment in Neuroblastoma" manuscript.

Summary excel spreadsheet has been provided to indicate input and output files from each stage of analysis for each of use.

## Pipeline Structure

The pipeline is structured as a series of SLURM jobs, allowing for efficient resource management and parallel processing in a cluster environment. Followed by R analysis scripts to analysis downstream of processing:

1. **Analysis of Cellecta Clonal Dynamics at Single Cell Resolution (ClonalDynamics_sc.R):**
   - Analysis of cellecta clones under different treatment conditions at single cell resolution to understand plasticity systems.


## Running the Pipeline

To run the pipeline, run R script (locally or on cluster) to conduct full analysis pipeline.

ClonalDynamics_sc.R
