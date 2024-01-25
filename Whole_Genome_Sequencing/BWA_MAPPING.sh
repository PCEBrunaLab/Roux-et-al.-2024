#!/bin/bash
#SBATCH --job-name=bwa
#SBATCH --ntasks=1
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=120:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Parent directory where all sample folders are located
PARENT_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WGS/clean

# Other paths
REF_GENOME=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/GenomeReference/Homo_sapiens_assembly38.fasta

# Load BWA and Samtools
module load BWA/0.7.17
module load SAMtools/1.11

# Iterate over all directories (samples) in the parent directory
for DIR in $PARENT_DIR/*; do
  if [ -d "$DIR" ]; then
    # Set sample-specific paths
    READ1=$(find $DIR -type f -name '*_1.fq.gz')  # Searches for any file ending with '*_1.fq.gz'
    READ2=$(find $DIR -type f -name '*_2.fq.gz')  # Searches for any file ending with '*_2.fq.gz'
    BASENAME=$(basename "$READ1" _1.fq.gz)

    PREFIX=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WGS/BWA/${BASENAME}  

    # Check if input files exist
    if [ ! -f "$READ1" ] || [ ! -f "$READ2" ]; then
      echo "Input files do not exist for $DIR"
      continue
    fi

    # Run BWA
    bwa mem -t $SLURM_CPUS_PER_TASK $REF_GENOME $READ1 $READ2 | samtools view -Sb - > ${PREFIX}.bam

  fi
done
