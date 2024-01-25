#!/bin/bash
#SBATCH --job-name=MuTect2
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
ReadGroup_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/ReadGroup
reference_genome=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/GenomeReference/Homo_sapiens_assembly38.fasta
germline_resource=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/Germline/af-only-gnomad.hg38.vcf

# Change to the job directory
cd $ReadGroup_dir

# Load GATK module
conda activate gatk4.1.9.0

# Iterate over all bam files in the ReadGroup directory
for BAM_FILE in $ReadGroup_dir/*recalibrated.bam; do
  # Get the base name of the bam file (without the .bam extension)
  BASENAME=$(basename "$BAM_FILE" .bam)

  # Call variants with MuTect2
  gatk Mutect2 \
    -R ${reference_genome} \
    -I ${BAM_FILE} \
    --germline-resource ${germline_resource} \
    -O ${ReadGroup_dir}/${BASENAME}.vcf.gz

done

