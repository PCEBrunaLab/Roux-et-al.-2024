#!/bin/bash
#SBATCH --job-name=GermlineFilter_Strelka
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8000
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Load necessary modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init
module load BCFtools/1.11

# Paths
ReadGroup_dir=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/ReadGroup/
germline_resource=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/Germline/af-only-gnomad.hg38.vcf.gz

# Change to the job directory
cd $VARIANTS_DIR

# Process the VCF files
for VCF_FILE in $ReadGroup_dir/results/variants/variants.vcf.gz; do
    OUTPUT_VCF=${VCF_FILE%.vcf.gz}_filtered_germline.vcf.gz

    # Exclude known germline variants from Strelka2 VCF
    bcftools isec -w 1 -C -Oz -o $OUTPUT_VCF $VCF_FILE $germline_resource

    # Index the filtered VCF
    tabix -p vcf $OUTPUT_VCF
done
