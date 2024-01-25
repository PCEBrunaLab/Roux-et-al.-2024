#!/bin/bash
#SBATCH --job-name=Samtools
#SBATCH --ntasks=1
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Paths
SCRATCH_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata
STAR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/STAR

# Change to the job directory
cd $SCRATCH_DIR

# Running STAR
module load SAMtools/1.11

# Iterate over all bam files in the STAR directory
for BAM_FILE in $STAR/*.bam; do
  # Get the base name of the bam file
  BASENAME=$(basename "$BAM_FILE" .bam)

  # Index the bam file with samtools and output to the Samtools directory
  samtools index ${BAM_FILE} ${BASENAME}.bam.bai
done