#!/bin/bash
#SBATCH --job-name=Mutect_FilterCriteria
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load required modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init

# Activate GATK environment
conda activate gatk4.1.9.0

# Define paths
reference_genome="/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/GenomeReference/Homo_sapiens_assembly38.fasta"
Filter_dir="/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/VCF_Files"

# Loop through all filtered VCF files and apply criteria
for FILTERED_VCF_FILE in $Filter_dir/*_PassOnly.vcf.gz
do
    # Extract base name of the VCF file
    BASENAME_FILTERED=$(basename "$FILTERED_VCF_FILE" _PassOnly.vcf.gz)  
  
    # Use GATK's VariantFiltration for applying filter criteria
    gatk VariantFiltration \
        -R ${reference_genome} \
        -V ${FILTERED_VCF_FILE} \
        --filter-expression "DP < 30 || TLOD < 6.3 || AF < 0.01 || MBQ < 25 || MMQ < 40 || SBF > 0.9" \
        --filter-name "DepthOrFreqFilter" \
        -O ${Filter_dir}/${BASENAME_FILTERED}_CriteriaFiltered.vcf.gz

    # Extract variants that pass the new filters
    gatk SelectVariants \
        -R ${reference_genome} \
        -V ${Filter_dir}/${BASENAME_FILTERED}_CriteriaFiltered.vcf.gz \
        --exclude-filtered \
        -O ${Filter_dir}/${BASENAME_FILTERED}_FinalFiltered.vcf.gz
done

