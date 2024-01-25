#!/bin/bash
#SBATCH --job-name=GATK_bundle
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --output=GATK_bundle_%j.out
#SBATCH --error=GATK_bundle_%j.err

module load wget
module load gzip

# Define the working directory
work_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES
GATK_bundle=${work_dir}/GenomeReference

# Change to the GATK_bundle directory
cd ${GATK_bundle}

# Download GATK resource bundle
wget --user=gsapubftp-anonymous --password= --recursive --no-parent --no-directories --no-check-certificate ftp://ftp.broadinstitute.org/bundle/hg38/

# Decompress the files
gunzip *.gz
