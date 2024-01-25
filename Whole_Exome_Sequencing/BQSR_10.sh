#!/bin/bash
#SBATCH --job-name=Apply_BQSR
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init

# Paths
ReadGroup_Dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/ReadGroup
REF_GENOME=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/GenomeReference/Homo_sapiens_assembly38.fasta

# Load GATK module
conda activate gatk4.1.9.0

# Iterate over all bam files in the SCRATCH_DIR directory
for BAM_FILE in $ReadGroup_Dir/*.bam; do
  # Get the base name of the bam file (without the .bam extension)
  BASENAME=$(basename "$BAM_FILE" .bam)

  # Apply the recalibration with GATK
  gatk ApplyBQSR -R $REF_GENOME -I $BAM_FILE --bqsr-recal-file ${ReadGroup_Dir}/${BASENAME}_recal.table -O ${ReadGroup_Dir}/${BASENAME}_recalibrated.bam
done
