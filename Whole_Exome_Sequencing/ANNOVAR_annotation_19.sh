#!/bin/bash
#SBATCH --job-name=ANNOVAR_annotation
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
humandb_path=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/WES/Annovar/humandb 

# Annotate variants
for AVINPUT_FILE in $Annovar_dir/*.avinput; do
    OUTPUT_PREFIX=$(basename $AVINPUT_FILE ".avinput")
    table_annovar.pl $AVINPUT_FILE $humandb_path -buildver hg38 \
    -out ${Annovar_dir}/${OUTPUT_PREFIX} \
    -protocol refGene,knownGene, \
    -operation g,g, \
    -nastring . \
    -csvout
done
