#! /usr/bin/bash

## process Cellecta barcode libraries

BASE_DIR="/nfs/research/marioni/mdmorgan/Neuroblastoma_SingleCell"
SRC=$(echo $BASE_DIR"/src")
FASTDIR=$(echo $BASE_DIR"/J60_pools_1_2_3_run1315/fastq_path/J60/")
WHITELIST=$(echo $BASE_DIR"/barcodes/Whitelist_14nt.txt")
WHITE30=$(echo $BASE_DIR"/barcodes/Whitelist_30nt.txt")

INPREF=("A" "B" "C")
#INPREF=("A")

## define a joining functino
function join_by {
    local IFS="$1"
    shift
    echo "$*"
}

for SAMP in ${INPREF[@]};
do
    JNAME="Neuro_"$SAMP"_Cellecta_process"
    ERROR=$BASE_DIR"/logs/Neuro_"$SAMP"_Cellecta_process.err"
    OUTLOG=$BASE_DIR"/logs/Neuro_"$SAMP"_Cellecta_process.out"

    # Index 1 files
    I1FIND=()
    while IFS= read -r -d $'\0'; do
	  I1FIND+=("$REPLY")
    done < <(find $FASTDIR"/" -name "Cellecta_"$SAMP"*I1*fastq.gz" -print0)

    I1FILES=$(join_by , "${I1FIND[@]}")

    # Index 2 files
    I2FIND=()
    while IFS= read -r -d $'\0'; do
	  I2FIND+=("$REPLY")
    done < <(find $FASTDIR"/" -name "Cellecta_"$SAMP"*I2*fastq.gz" -print0)

    I2FILES=$(join_by , "${I2FIND[@]}")

    # Read 1 files
    R1FIND=()
    while IFS= read -r -d $'\0'; do
	R1FIND+=("$REPLY")
    done < <(find $FASTDIR"/" -name "Cellecta_"$SAMP"*R1*fastq.gz" -print0)

    R1FILES=$(join_by , "${R1FIND[@]}")

    # Read 2 files
    R2FIND=()
    while IFS= read -r -d $'\0'; do
	  R2FIND+=("$REPLY")
    done < <(find $FASTDIR"/" -name "Cellecta_"$SAMP"*R2*fastq.gz" -print0)

    R2FILES=$(join_by , "${R2FIND[@]}")
    
    OUTFILE=$(echo $BASE_DIR"/barcodes/"$SAMP"_barcode_freq.txt.gz")

    JOB="bsub -q bigmem -M 62000 -R "rusage[mem=60000]" -T 120 -n 24 -J $JNAME -e $ERROR -o $OUTLOG python $SRC/process_sc_barcodes.py --I1 $I1FILES --I2 $I2FILES --R1 $R1FILES --R2 $R2FILES --whitelist-14nt $WHITELIST --whitelist-30nt $WHITE30 --output $OUTFILE --processes 24 --chunk-size 100000"

    echo $JOB
    eval $JOB
done
    
	  
