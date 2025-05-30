#! /usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=8000

# Set working directory
BASE_DIR=...
cd $BASE_DIR

# load R
# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate ...

## per sample
SRCDIR=$(echo $BASE_DIR/src)
LOGDIR=$(echo $BASE_DIR/logs)
OUTDIR=$(echo $BASE_DIR/SCE.dir)

SCE=$OUTDIR/"Neuro-A_SCE.RDS"

SNAME=Neuro_A
OUTFILE=$(echo $OUTDIR/Neuro_A)
ID=Neuro_A

srun Rscript $SRCDIR"/cmo_stats.R" --SCE=$SCE --id=$ID --output=$OUTFILE 
