#BWA Index Generation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required.packages <- c("tidyverse", "data.table", "qdapRegex")

if (!all(required.packages %in% installed.packages()[,"Package"])){
  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]])
}

library(tidyverse)
library(data.table)
library(qdapRegex)

#Read in csv files for all 14bp and 30bp Cellecta barcode sequence within pool used.
#These can be found in Cellecta-SEQ-CloneTracker-XP-10M-Barcodes.xlsx file for associated pool
bc.14 <- readLines("bc.14.csv")
bc.30 <- readLines("bc.30.csv")

#Collate all possible sequences of bc.14 and bc.30 with TGGT as spacer between
test <- apply(expand.grid(bc.14, bc.30), 1, paste, collapse="TGGT")
barcodes <- as.data.frame(test)
write.csv(barcodes, "full_barcode.csv")

#Seperately manually add barcode names to each barcode; barcode 1 - n name associated with each barcode
#Re-upload named.csv files to scratch fpr downstream processing
listed_barcodes <- read.csv("full_barcodes_named.csv", header = T)

#Collate all possible sequences with flanking sequence for more accurate alignment
for (i in 1:nrow(listed_barcodes)) {
  write(paste0('>', listed_barcodes$name[i]),"barcodes_full.fasta",append=TRUE,sep="")
  write(paste0(listed_barcodes$sequence[i]), "barcodes_full.fasta", append = TRUE, sep = "")
}

#Generate index in terminal using following command:
bwa index barcodes_full.fasta

#The two files .ann and .pac must then be altered to add all the barcode sequences to the end, otherwise these will not be used as references for alignment during bwa aligner.
cellecta.ann <- readLines("barcodes_full.fasta.ann", warn=FALSE)
cellecta.pac <- readLines("barcodes_full.fasta.pac", warn=FALSE)

#To add barcodes to these files we must add first line >sequence1 followed by line 2 with sequence1. 
## Generate a .txt file with this information to add to the .ann and .pac
#Alter pool 1 .ann file
for (i in 1:nrow(listed_barcodes_pool_1)) {
  write(paste0(">", listed_barcodes_pool_1$sequence[i]), "barcodes_full_pool_1.fasta.ann", append = TRUE)
  write(listed_barcodes_pool_1$sequence[i], "barcodes_full_pool_1.fasta.ann", append = TRUE)
}

#Alter pool 1 .pac file
for (i in 1:nrow(listed_barcodes_pool_1)) {
  write(paste0(">", listed_barcodes_pool_1$sequence[i]), "barcodes_full_pool_1.fasta.pac", append = TRUE)
  write(listed_barcodes_pool_1$sequence[i], "barcodes_full_pool_1.fasta.pac", append = TRUE)
}


