#! /usr/bin/bash
#SBATCH --job=cellranger_multi
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=18:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J60_singlecell/logs/%x.%j_.err
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J60_singlecell/logs/%x.%j_.out

#SBATCH --mail-user=sian.hamer@icr.ac.uk
#SBATCH --mail-type=ALL

# Set working directory
BASE_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J60_singlecell
cd $BASE_DIR

# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate hotspot

#Set path to cellranger
export PATH=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/cellranger/apps/cellranger-7.1.0:$PATH

## submit a cellranger job
cellranger multi --id=Neuro_A --csv=config_A.csv
