#! /usr/bin/env Rscript

## compute some summaries over the CMOs for each sample
#load libraries
library(SingleCellExperiment)
library(DropletUtils)
library(scran)
library(scater)
library(Matrix)
library(MatrixGenerics)
library(BiocParallel)
library(biomaRt)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the SCE file")

parser <- add_option(parser, c("-i", "--id"), type="character",
                     help="Sample ID")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output combined SCE object")

opt <- parse_args(parser)

message("Reading SCE file ", opt$SCE)
in.sce <- readRDS(opt$SCE)
colnames(in.sce) <- paste(colData(in.sce)$Sample, colData(in.sce)$Barcode, sep="_")

cmo.feats <- rownames(in.sce)[grepl(rownames(in.sce), pattern="CMO[0-9]+")]
message("Found ", length(cmo.feats), " CMOs")
cmo.counts <- counts(in.sce[cmo.feats, ])

cmo.sum <- data.frame("CMO"=cmo.feats, "Sample"=opt$id,
	              "Counts"=rowSums2(cmo.counts)) # proportion of cells in each cmo
cell.sum <- data.frame("CellID"=colnames(in.sce), "Sample"=opt$id,
                       "Counts"=colSums2(cmo.counts)) # proportion of cmos in each cell

cmo.ofile <- paste0(opt$output, "_CMOsummary.tsv")
message("Writing CMO summary to ", cmo.ofile)
write.table(cmo.sum, file=cmo.ofile, sep="\t", quote=FALSE, row.names=FALSE)

cell.ofile <- paste0(opt$output, "_Cellsummary.tsv")
message("Writing cell summary to ", cell.ofile)
write.table(cell.sum, file=cell.ofile, sep="\t", quote=FALSE, row.names=FALSE)

message("All done")




