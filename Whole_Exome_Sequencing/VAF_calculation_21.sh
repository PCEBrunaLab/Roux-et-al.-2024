#!/bin/bash
#SBATCH --job-name=VAF_calculation
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load necessary modules
module load BCFtools/1.11

# Paths
Filter_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/FilterCriteria

# Change to the job directory
cd $Filter_dir

for VCF_FILE in $Filter_dir/*_Targeted.vcf.gz; do
    # Extract base name of the VCF file for output naming
    BASENAME=$(basename "$VCF_FILE" _Targeted.vcf.gz)
    TEMP_FILE="${Filter_dir}/${BASENAME}_temp.txt"
    OUTPUT_FILE="${Filter_dir}/${BASENAME}_VAF.txt"
  
    # Extract INFO and FORMAT fields to the temporary file
    bcftools query -f '%CHROM:%POS:%REF:%ALT[\t%AD\t%DP]\n' $VCF_FILE > $TEMP_FILE

    # Calculate VAF using awk and write to the output file
    awk -F'\t' '{ split($2, depth, ","); vaf = depth[2] / $3; print $1 "\t" vaf; }' $TEMP_FILE > $OUTPUT_FILE

    # Optionally, remove the temporary file
    rm $TEMP_FILE

done




