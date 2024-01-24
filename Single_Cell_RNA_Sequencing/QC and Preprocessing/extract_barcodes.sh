#! /usr/bin/bash

## process Cellecta barcode libraries
BASE_DIR="/nfs/research/marioni/mdmorgan/Neuroblastoma_SingleCell"
SRC=$(echo $BASE_DIR"/src")

INPREF=("Neuro_A" "Neuro_B" "Neuro_C" "Neuro_S")
#INPREF=("Neuro_C")
OUTDIR=$BASE_DIR"/cellbarcodes.dir"

## define a joining function
function join_by {
    local IFS="$1"
    shift
    echo "$*"
}

for SAMP in ${INPREF[@]};
do
    JNAME=$SAMP"_extract_barcodes"
    ERROR=$BASE_DIR"/logs/"$SAMP"_extract_barcodes.err"
    OUTLOG=$BASE_DIR"/logs/"$SAMP"_extract_barcodes.out"

    BAMDIR=$SAMP"/SC_MULTI_CS/SC_MULTI_CORE/STRUCTIFY_PER_SAMPLE_OUTS/fork0/chnk0/files"

    BAMFIND=()
    while IFS= read -r -d $'\0'; do
	BAMFIND+=("$REPLY")
    done < <(find $BAMDIR"/" -name "*bam" -print0 )

    BAMFILES=$(join_by , "${BAMFIND[@]}")
    echo $BAMFILES

    OUTFILE=$OUTDIR/$SAMP"_bc.txt.gz"

    JOB="bsub -q standard -M 40000 -R "rusage[mem=38000]" -e $ERROR -o $OUTLOG -J $JNAME -T 1500 python $SRC/extract_barcodes.py --bam $BAMFILES --output $OUTFILE"

    echo $JOB
    eval $JOB    
done
