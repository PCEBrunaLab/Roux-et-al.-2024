#! /usr/bin/bash
#SBATCH --job=emptyDrops
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=12000
#SBATCH --partition=smp
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J60_singlecell/logs/%x.%j-.out
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J60_singlecell/logs/%x.%j-.err

#SBATCH --mail-user=sian.hamer@icr.ac.uk
#SBATCH --mail-type=ALL

# Set working directory
BASE_DIR=/data/scratch/DMP/DUDMP/PAEDCANC/shamer/J60_singlecell
cd $BASE_DIR

# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate hotspot

## loop over sample and run emptyDrops + cell filtering on each sample separately - returns an SCE object
## per sample
SRCDIR=$(echo $BASE_DIR/src/)
OUTDIR=$(echo $BASE_DIR/SCE.dir/)
UMI=500 # threshold for background empty droplets - this will depend heavily on the sequencing depth of the samples

sampdir=$BASE_DIR/Neuro_A/

H5FILE=$(find $sampdir/outs/multi/count/ -name "*matrix.h5")
OUTFILE=$OUTDIR/"Neuro-A_SCE.RDS"

srun Rscript $SRCDIR"/test_emptyDrops.R" --h5file=$H5FILE --id=$sampdir --ncores=12 --umithreshold=$UMI --output=$OUTFILE