#!/bin/bash
#SBATCH --job-name=Readgroup
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Paths
Picard_Dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/Picard_output
ReadGroup_Dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/ReadGroup

# Change to the job directory
cd $Picard_Dir

# Load Picard module
module load picard-tools/2.23.8

# Iterate over all bam files in the Picard directory
for BAM_FILE in $Picard_Dir/*.bam; do
  # Get the base name of the bam file 
  BASENAME=$(basename "$BAM_FILE" .bam)

  # Create output BAM filename
  output_file="$ReadGroup_Dir/${BASENAME}_RG.bam"

  # Run Picard's AddOrReplaceReadGroups
  picard AddOrReplaceReadGroups \
  I="$BAM_FILE" \
  O="$output_file" \
  RGID=group1 \
  RGLB=lib1 \
  RGPL=illumina \
  RGPU=unit1 \
  RGSM=sample1

done
