#!/bin/bash
#SBATCH --job-name=HTSeq
#SBATCH --ntasks=1
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16000
#SBATCH --cpus-per-task=48
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init

# Set the directory where your bam files are located
STAR_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/STAR

# Set the directory where your output files will be located
HTSEQ_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/HTSeq
mkdir -p $HTSEQ_DIR

# Path to your GTF file
GTF_PATH=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/GenomeReference/Homo_sapiens.GRCh38.109.gtf

# Activate the conda environment with HTSeq
conda activate htseq0.12.4

# Iterate over all bam files in the STAR directory
for BAM_FILE in $STAR_DIR/*.bam; do
  # Get the base name of the bam file 
  BASENAME=$(basename "$BAM_FILE" .bam)

# Run HTSeq count and output to the HTSeq directory
  htseq-count \
  --format bam \
  --stranded yes \
  --order pos \
  --additional-attr gene_name \
  ${BAM_FILE} ${GTF_PATH} > ${HTSEQ_DIR}/${BASENAME}_counts.txt
done
