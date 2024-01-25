#!/bin/bash
#SBATCH --job-name=star
#SBATCH --ntasks=1
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16000
#SBATCH --cpus-per-task=12
#SBATCH --time=120:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Parent directory where all sample folders are located
PARENT_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/clean/

# Other paths
GENOMEDIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/GenomeReference
ANNOTATION=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/GenomeReference/annotation.gtf
STAR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping

# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init
SCRATCH_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata

# Change to the job directory
cd $SCRATCH_DIR

# Running STAR
conda activate star2.7.6a


# Iterate over all directories (samples) in the parent directory
for DIR in $PARENT_DIR/*; do
  if [ -d "$DIR" ]; then
    # Set sample-specific paths
    READ1=$(find $DIR -type f -name '*_1.fq.gz')  # Searches for any file ending with '*_1.fq.gz'
    READ2=$(find $DIR -type f -name '*_2.fq.gz')  # Searches for any file ending with '*_2.fq.gz'
    BASENAME=$(basename "$READ1" _1.fq.gz)

    PREFIX=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/${BASENAME}_out  # STAR outputs will be written here

    # Check if input files exist
    if [ ! -f "$READ1" ] || [ ! -f "$READ2" ]; then
      echo "Input files do not exist for $DIR"
      continue
    fi

    # Run STAR
    STAR --runThreadN $SLURM_CPUS_PER_TASK \
    --runMode alignReads \
    --genomeDir $GENOMEDIR \
    --readFilesIn $READ1 $READ2 \
    --outFileNamePrefix ${PREFIX} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate 

    # Rename the BAM file to have the base name of the FASTQ file
    mv "${BASENAME}_outAligned.sortedByCoord.out.bam" "${STAR}/${BASENAME}.bam"
    fi
done