#!/bin/bash
#SBATCH --job-name=TARGET_BED_SureSelect_Strelka
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
conda activate gatk4.1.9.0

# Define paths
reference_genome=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/GenomeReference/Homo_sapiens_assembly38.fasta
Filter_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/mutect
TARGET_BED=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/GenomeReference/S07604514_Regions.bed

# Loop through all directories matching *_mutect_run in the Filter_dir
for DIR in $Filter_dir/*_mutect_run; do
    # Assuming the directory and file exist, construct the file path
    FILTERED_VCF_FILE="${DIR}/results/variants/*_FinalFiltered.vcf.gz"
    BASENAME_FILTERED=$(basename "$DIR")

    # Intersect VCF with BED
    gatk SelectVariants \
        -R ${reference_genome} \
        -V ${FILTERED_VCF_FILE} \
        -L ${TARGET_BED} \
        -O "${DIR}/${BASENAME_FILTERED}_Targeted.vcf.gz"
done
