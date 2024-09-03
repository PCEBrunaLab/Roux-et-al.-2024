#! /usr/bin/bash

module load r-4.1.1-gcc-9.3.0-jkdw35f

BASEDIR=...
SRC=$(echo $BASEDIR"/src")

INDIR=$(echo $BASEDIR"/barcodes")
ALLFILES=$(find $INDIR -name "*_barcode_freq.txt.gz" | cut -d "/" -f 8 | cut -d "_" -f 1 | sort | uniq)
echo $ALLFILES
ARRAY=( $(for x in $ALLFILES; do echo $x; done) )
LEN=${#ARRAY[@]}

for FILE in ${ARRAY[@]}; do

    INFILE=$(echo $INDIR"/"$FILE"_barcode_freq.txt.gz")
    OUTPREF=$(echo $BASEDIR"/corrected/Neuro_$FILE")
    ERROR=$(echo $BASEDIR"/logs/Neuro_$FILE-Error_correction.err")
    OUTLOG=$(echo $BASEDIR"/logs/Neuro_$FILE-Error_correction.out")

    JOB="bsub -q standard -M 28000 -R "rusage[mem=26000]" -T 1400 -J 'Error_correct-$FILE' -e $ERROR -o $OUTLOG  Rscript $SRC/barcode_error_correct.R --file $INFILE --output $OUTPREF"

    echo $JOB
    eval $JOB
done

