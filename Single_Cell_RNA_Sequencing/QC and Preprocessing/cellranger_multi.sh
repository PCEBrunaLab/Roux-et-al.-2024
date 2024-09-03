#! /usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=18:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute


# Set working directory
BASE_DIR=...
cd $BASE_DIR

# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate ...

#Set path to cellranger
export PATH=...cellranger-7.1.0:$PATH

## submit a cellranger job
cellranger multi --id=Neuro_A --csv=config_A.csv
