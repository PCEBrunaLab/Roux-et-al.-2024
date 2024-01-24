#! /usr/bin/bash
#SBATCH --job=bwa_alignment
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=30:00:00
#SBATCH --partition=compute
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J63_DNA_sequencing/logs/%x.%j-.out
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J63_DNA_sequencing/logs/%x.%j-.err

#SBATCH --mail-user=sian.hamer@icr.ac.uk
#SBATCH --mail-type=ALL

# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate hotspot

# Set working directory
BASE_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J63_DNA_sequencing
cd $BASE_DIR

srun Rscript $BASE_DIR"/bwa_alignment.R"
