#! /usr/bin/bash
#SBATCH --job=Barcode_Extraction_Single_Cell8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=2:30:00
#SBATCH --partition=compute
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/M62_singlecell/scRNA-seq/cellecta/logs/%x.%j-.out
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/M62_singlecell/scRNA-seq/cellecta/logs/%x.%j-.err

#SBATCH --mail-user=sian.hamer@icr.ac.uk
#SBATCH --mail-type=ALL

# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate hotspot

# Set working directory
BASE_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/M62_singlecell/scRNA-seq/cellecta
cd $BASE_DIR

OUTPUT_DIR=$(echo $BASE_DIR/outputs/)
FASTQ=$(echo $BASE_DIR/fastq_in_use1/)

srun Rscript $BASE_DIR"/cellecta1_barcode_extractions_scRNA_NB.R" --sample=$FASTQ
