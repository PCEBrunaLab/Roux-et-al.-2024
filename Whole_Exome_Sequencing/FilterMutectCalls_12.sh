#!/bin/bash
#SBATCH --job-name=FilterMutectCalls
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
vcf_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/ReadGroup
reference_genome=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/GenomeReference/Homo_sapiens_assembly38.fasta
Filter_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/VCF_Files
SNPEFF=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/SNPEFF

# Change to the job directory
cd $vcf_dir

# Load GATK module
conda activate gatk4.1.9.0


# Loop over all VCF files
for VCF_FILE in $vcf_dir/*.vcf.gz
do
  # Get the base name of the vcf file (without the .vcf.gz extension)
  BASENAME=$(basename "$VCF_FILE" .vcf.gz)
  
  # Filter variants with FilterMutectCalls
  gatk FilterMutectCalls \
    -R ${reference_genome} \
    -V ${VCF_FILE} \
    -O ${Filter_dir}/${BASENAME}_Filtered.vcf.gz  
  
done

