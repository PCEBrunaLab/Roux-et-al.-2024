#! /usr/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=12000

# Set working directory
BASE_DIR=...
cd $BASE_DIR

# Activate mamba env
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
mamba activate ...

## loop over sample and run emptyDrops + cell filtering on each sample separately - returns an SCE object
## per sample
SRCDIR=$(echo $BASE_DIR/src/)
OUTDIR=$(echo $BASE_DIR/SCE.dir/)
UMI=500 # threshold for background empty droplets - this will depend heavily on the sequencing depth of the samples

#conduct per Multiplexed pool
sampdir=$BASE_DIR/Neuro_A/

H5FILE=$(find $sampdir/outs/multi/count/ -name "*matrix.h5")
OUTFILE=$OUTDIR/"Neuro-A_SCE.RDS"

srun Rscript $SRCDIR"/test_emptyDrops.R" --h5file=$H5FILE --id=$sampdir --ncores=12 --umithreshold=$UMI --output=$OUTFILE
