#!/bin/bash
#SBATCH --job-name=Samtools_index_RG
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Paths
ReadGroup_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/ReadGroup

# Change to the job directory
cd $ReadGroup_dir

# Load Samtools module
module load SAMtools/1.11

# Iterate over all bam files in the BWA directory
for BAM_FILE in $ReadGroup_dir/*recalibrated.bam; do
  # Get the base name of the bam file 
  BASENAME=$(basename "$BAM_FILE" .bam)

  # Check the number of reads and coverage
  samtools index ${BAM_FILE} ${ReadGroup_dir}/${BASENAME}.bai


done
