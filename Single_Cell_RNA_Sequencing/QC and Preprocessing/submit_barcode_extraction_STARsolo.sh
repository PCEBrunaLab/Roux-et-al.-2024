#! /usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=2:30:00
#SBATCH --partition=compute

# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate ...

# Set working directory
BASE_DIR=...
cd $BASE_DIR

OUTPUT_DIR=$(echo $BASE_DIR/outputs/)
FASTQ=$(echo $BASE_DIR/fastq/)

srun Rscript $BASE_DIR"/barcode_extractions_STARsolo.R" --sample=$FASTQ
