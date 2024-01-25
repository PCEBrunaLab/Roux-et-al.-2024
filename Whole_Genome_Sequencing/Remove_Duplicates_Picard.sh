#!/bin/bash
#SBATCH --job-name=Remove_Duplicate
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Paths
SCRATCH_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WGS/BWA
Picard=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WGS/Picard

# Change to the job directory
cd $SCRATCH_DIR

# Load Samtools module
module load picard-tools/2.23.8

# Iterate over all bam files in the SCRATCH_DIR directory
for BAM_FILE in $SCRATCH_DIR/*_sorted.bam; do
  # Get the base name of the bam file (without the .bam extension)
  BASENAME=$(basename "$BAM_FILE" .bam)

  # Remove duplicates with Picard and save output to Picard directory
  picard MarkDuplicates INPUT=$BAM_FILE OUTPUT=${Picard}/${BASENAME}_dedup_reads.bam METRICS_FILE=${Picard}/${BASENAME}_output.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=16000

done


