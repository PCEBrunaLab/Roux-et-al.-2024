#!/bin/bash
#SBATCH --job-name=Mutect_SelectVariants
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init

# Paths
vcf_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/VCF_Files
reference_genome=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/GenomeReference/Homo_sapiens_assembly38.fasta
Filter_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/VCF_Files

# Change to the job directory
cd $vcf_dir

# Load GATK module
conda activate gatk4.1.9.0

# Loop over all Filtered VCF files
for FILTERED_VCF_FILE in $Filter_dir/*_Filtered.vcf.gz
do
  # Get the base name of the filtered vcf file (without the _Filtered.vcf.gz extension)
  BASENAME_FILTERED=$(basename "$FILTERED_VCF_FILE" _Filtered.vcf.gz)
  
  # Select variants with "PASS" using GATK's SelectVariants
  gatk SelectVariants \
    -R ${reference_genome} \
    -V ${FILTERED_VCF_FILE} \
    --select-type-to-include SNP \
    --select-type-to-include INDEL \
    --exclude-filtered \
    -O ${Filter_dir}/${BASENAME_FILTERED}_PassOnly.vcf.gz
  
done