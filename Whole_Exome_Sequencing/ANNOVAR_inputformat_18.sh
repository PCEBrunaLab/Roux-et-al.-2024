#!/bin/bash
#SBATCH --job-name=ANNOVAR_inputformat
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load Annovar
module load ANNOVAR/20200607

# Paths
Annovar_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/Annovar
Filter_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/FilterCriteria

# Change to the directory
cd $Annovar_dir

# Convert VCF to ANNOVAR input format
for VCF_FILE in $Filter_dir/*_Targeted.vcf.gz; do
    OUTPUT_NAME=$(basename $VCF_FILE ".vcf.gz")
    convert2annovar.pl -format vcf4 $VCF_FILE > $Annovar_dir/${OUTPUT_NAME}.avinput

done