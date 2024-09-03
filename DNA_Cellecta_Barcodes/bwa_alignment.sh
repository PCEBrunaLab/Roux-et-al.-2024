#! /usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=30:00:00

# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate ...

# Set working directory
BASE_DIR=...
cd $BASE_DIR

srun Rscript $BASE_DIR"/bwa_alignment.R"
