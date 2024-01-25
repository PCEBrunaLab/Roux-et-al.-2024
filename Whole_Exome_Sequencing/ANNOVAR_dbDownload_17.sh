#!/bin/bash
#SBATCH --job-name=ANNOVAR_dbDownload
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

# Change to the directory
cd $Annovar_dir

# Download the refGene database for hg38
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar knownGene humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad_genome humandb/
annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
