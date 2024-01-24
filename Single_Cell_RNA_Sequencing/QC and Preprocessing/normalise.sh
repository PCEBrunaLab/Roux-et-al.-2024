#! /usr/bin/bash
#SBATCH --job=Neuro_norm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --partition=compute
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J60_singlecell/logs/%x.%j.err
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J60_singlecell/logs/%x.%j.out

# Set working directory
BASE_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J60_singlecell
cd $BASE_DIR

# load R
# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate hotspot

## loop over samples and extract CMO summary information
## per sample
SRCDIR=$(echo $BASE_DIR/src)
LOGDIR=$(echo $BASE_DIR/logs)
OUTDIR=$(echo $BASE_DIR/SCE.dir)

# define input and output files
SCE=$(echo $BASE_DIR"/SCE.dir/Neuro-A_edited_SCE.RDS,"$BASE_DIR"/SCE.dir/Neuro-B_edited_SCE.RDS,"$BASE_DIR"/SCE.dir/Neuro-C_edited_SCE.RDS,"$BASE_DIR"/SCE.dir/Neuro-S_edited_SCE.RDS")
CALLS=$(echo $BASE_DIR"/SCE.dir/Neuro_A_demuxed.tsv,"$BASE_DIR"/SCE.dir/Neuro_B_demuxed.tsv,"$BASE_DIR"/SCE.dir/Neuro_C_demuxed.tsv,"$BASE_DIR"/SCE.dir/Neuro_S_demuxed.tsv")
OUTPUT=$(echo $BASE_DIR"/SCE.dir/Neuro_SCE-norm.RDS")

srun Rscript $SRCDIR"/norm_counts.R" --SCE=$SCE --breaks --calls=$CALLS --output=$OUTPUT

## also normalise the cell lines separately

