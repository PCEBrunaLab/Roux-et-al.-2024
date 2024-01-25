#!/bin/bash
#SBATCH --job-name=QC
#SBATCH --ntasks=1
#SBATCH --partition=compute
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

# Changing directories to where the fastq files are located
cd /data/scratch/DMP/DUDMP/PAEDCANC/asadr/WGS/clean

# Running FASTQC
module load FastQC/0.11.9
mkdir -p fastqc
# Find all of the fastq files
find /data/scratch/DMP/DUDMP/PAEDCANC/asadr/WGS/clean -name "*.gz" -exec fastqc -o fastqc/ {} \;

# run multiqc to compile individual fastqc files, this helps visualization of fastqc reports
conda activate multiqc1.9
#mkdir -p multiqc
multiqc /data/scratch/DMP/DUDMP/PAEDCANC/asadr/WGS/clean/ -o multiqc
