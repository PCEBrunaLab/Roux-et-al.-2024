#!/bin/bash
#SBATCH --job-name=GenomeIndex
#SBATCH --partition=smp
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=16000
#SBATCH --cpus-per-task=12
#SBATCH --time=2:00:00


# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda init

SCRATCH_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq


# Change to the job directory
cd $SCRATCH_DIR

# Changing directories to where the fastq files are located
cd /data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/GenomeReference

# Running STAR
conda activate star2.7.6a

# Run STAR

STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/GenomeReference \
--genomeFastaFiles /data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/GenomeReference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /data/scratch/DMP/DUDMP/PAEDCANC/asadr/BulkRNAseq/rawdata/Mapping/GenomeReference/Homo_sapiens.GRCh38.109.gtf \
--sjdbOverhang 99 \

