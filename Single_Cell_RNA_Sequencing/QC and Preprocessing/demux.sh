#! /usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=8000

### demultiplex with CellPlex CMOs

# Set working directory
BASE_DIR=...
cd $BASE_DIR

# load R and conda environment
# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate ..

#Conduct per multiplexed pool
SRCDIR=$(echo $BASE_DIR/src)
LOGDIR=$(echo $BASE_DIR/logs)
OUTDIR=$(echo $BASE_DIR/SCE.dir)
QUANT=0.99 # quantile threshold to declare positive CMO-assigned cells

SNAME=$BASE_DIR/SCE.dir/Neuro-A_SCE.RDS
OUTFILE=$OUTDIR/"Neuro_A_demuxed.tsv"

srun Rscript $SRCDIR"/demux.R" --SCE=$OUTDIR/"Neuro-A_SCE.RDS" --ncores=4 --quantile=$QUANT --output=$OUTFILE 
