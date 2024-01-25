#!/bin/bash
#SBATCH --job-name=ANNOVAR_to_MAF_conversion
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

# Convert ANNOVAR output to MAF format
for MULTIANNOT_FILE in $Annovar_dir/*.multianno.csv; do
    OUTPUT_PREFIX=$(basename $MULTIANNOT_FILE ".multianno.csv")
    
    # Convert the ANNOVAR output to MAF
    annovarToMaf.pl -i ${Annovar_dir}/${OUTPUT_PREFIX}.multianno.csv \
    -o ${Annovar_dir}/${OUTPUT_PREFIX}.maf
done
