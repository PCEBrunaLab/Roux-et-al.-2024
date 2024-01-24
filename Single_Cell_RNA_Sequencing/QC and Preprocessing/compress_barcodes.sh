#! /usr/bin/bash

## Take the unique observed barcodes
BASE_DIR="/nfs/research/marioni/mdmorgan/Neuroblastoma_SingleCell"
SRC=$(echo $BASE_DIR"/src")

INPREF=("Neuro_A" "Neuro_B" "Neuro_C" "Neuro_S")
#INPREF=("Neuro_C")
OUTDIR=$BASE_DIR"/cellbarcodes.dir"


for SAMP in ${INPREF[@]};
do
    JNAME=$SAMP"_extract_barcodes"
    ERROR=$BASE_DIR"/logs/"$SAMP"_unique_barcodes.err"
    OUTLOG=$BASE_DIR"/logs/"$SAMP"_unique_barcodes.out"

    BCFILE=$(find $OUTDIR -name "$SAMP*" | grep "bc")
    OUTFILE=$OUTDIR/$SAMP"_unique.txt.gz"

    JOB="bsub -q standard -M 40000 -R "rusage[mem=38000]" -e $ERROR -o $OUTLOG -J $JNAME -T 1500 python $SRC/compress_barcodes.py --file $BCFILE --output $OUTFILE"

    echo $JOB
    eval $JOB
done
