#!/bin/bash
#SBATCH --job-name=GenomeIndex
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12


# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init

SCRATCH_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES


# Change to the job directory
cd $SCRATCH_DIR

# Changing directories to where the genome files are located
cd /data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/GenomeReference

# Running BWA
module load BWA/0.7.17  

# Inex the reference file
bwa index -a bwtsw 	Homo_sapiens_assembly38.fasta


