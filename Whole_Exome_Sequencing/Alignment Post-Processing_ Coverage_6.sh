#!/bin/bash
#SBATCH --job-name=Samtools_Coverage
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Paths
SCRATCH_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/MERGED_WES
BWA=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/BWA

# Change to the job directory
cd $SCRATCH_DIR

# Load Samtools module
module load SAMtools/1.11

# Iterate over all bam files in the BWA directory
for BAM_FILE in $BWA/*.bam; do
  # Get the base name of the bam file 
  BASENAME=$(basename "$BAM_FILE" .bam)

  # Check the number of reads and coverage
  samtools flagstat ${BAM_FILE} > ${BWA}/${BASENAME}_coverage_mapping

done
