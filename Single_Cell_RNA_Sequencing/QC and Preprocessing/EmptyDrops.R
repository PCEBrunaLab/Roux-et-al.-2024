#!/usr/bin/env Rscript

# This script takes in all of the 10X experiments calls cells and normalizes them all together
# Hopefully this will help to mitigate against any major batch effects between the different experiments.
# Due to the size of the data I need to run the cell calling and filtering on each sample separately.

#install packages
required.packages <- c("devtools", "SingleCellExperiment", "DropletUtils", "scran", "Matrix", 
                       "MatrixGenerics", "BiocParallel", "biomaRt", "optparse", "scater", "Cairo")

if (!all(required.packages %in% installed.packages()[,"Package"])){
  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]])
}

#load libraries
library(devtools)
library(SingleCellExperiment)
library(DropletUtils)
library(scran)
library(scater)
library(Matrix)
library(MatrixGenerics)
library(BiocParallel)
library(biomaRt)
library(optparse)

start <- Sys.time()

parser <- OptionParser()
parser <- add_option(parser, c("-a", "--h5file"), type="character",
       	             help="The path to the HDF5 file")

parser <- add_option(parser, c("-i", "--id"), type="character",
                     help="Sample ID")

parser <- add_option(parser, c("-n", "--ncores"), type="numeric",
                     help="The number of cores to use for parallelised steps")

parser <- add_option(parser, c("-u", "--umithreshold"), type="numeric",
                     help="UMI threshold to define the background distribution of empty droplets")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output combined SCE object")

opt <- parse_args(parser)

message(paste0("Running cell calling on sample: ", opt$id))
# rather than relying on the CellRanger to call cells I'll use emptyDrops and remove poor quality cells at a later point.
mcparam <- MulticoreParam(workers=opt$ncores)
register(mcparam)

message("Reading from HDF5 file")
in.sce <- read10xCounts(samples=opt$h5file, sample.names=opt$id, BPPARAM=mcparam)

message("Running EmptyDrops to estimate cells with ", opt$ncores, " cores")
intest_calls <- emptyDrops(in.sce, assay.type="counts", niters=20000, ignore=4999, BPPARAM=mcparam,
                           lower=opt$umithreshold, retain=Inf)

message("Identifying cells at 1% FDR")
# identify cells using a FDR 1%
sig_cells <- intest_calls$FDR <= 0.01 & !is.na(intest_calls$FDR)

# subset the called cells - the numbers are very similar to what CellRanger gives
sce.drops <- in.sce[, sig_cells]
message("Keeping ", sum(sig_cells), " non-empty droplet barcodes")

######################
## Single-cell QC
######################
## Remove very low complexity cell libraries, i.e. ones that express < 1000 genes
lib.sizes <- MatrixGenerics::colSums2(counts(sce.drops))
nonzero.genes <- MatrixGenerics::colSums2(counts(sce.drops) > 0)

# remove cells with < 1000 genes
#sum(nonzero.genes > 1000)

message("Selecting cells with > 1000 UMIs")
sce.drops <- sce.drops[, nonzero.genes > 1000]

## remove cells with excessively high mitochondrial content
biomaRt.connection <- useMart("ensembl", "hsapiens_gene_ensembl")
gene.df <- getBM(attributes = c("ensembl_gene_id", "chromosome_name"), filters="ensembl_gene_id", 
                 values=rownames(sce.drops), mart=biomaRt.connection)
rownames(gene.df) <- gene.df$ensembl_gene_id

mt.counts <- counts(sce.drops)[which(rownames(sce.drops) %in% gene.df$ensembl_gene_id[gene.df$chromosome_name == "MT"]), ]
mt.fraction <- MatrixGenerics::colSums2(mt.counts)/MatrixGenerics::colSums2(counts(sce.drops))

message("Filtering cells with high MT content")
# fit a median-centred, MAD variance model normal to dervive p-values for outlier proportions
# this isn't the best calibrated model, so might be a little liberal
mt.p <- pnorm(mt.fraction, mean=median(mt.fraction), sd=mad(mt.fraction)*5, lower.tail=FALSE)
mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method="fdr") < 0.1)])
print(mt.lim)

# cells with a high mitochondrial fraction > `mt.lim` are removed as outliers.
message(paste0("Removing ", sum(mt.fraction >= mt.lim), " cells for high MT content"))
sce.drops <- sce.drops[, mt.fraction < mt.lim]

message(paste("Writing unnormalised SCE object to file:", opt$output)) 
# Now we have fairly decent cells I'll proceed to the normalisation using deconvolution-estimated size factors
saveRDS(sce.drops, file=opt$output) 
# Define filtered file name
filt.file <- gsub(opt$output, pattern="_SCE\\.RDS", replacement="_filterInfo.tsv") 
message(paste("Saving filtering info to file:", filt.file)) 
filt.df <- data.frame("N.NotEmpty"=sum(sig_cells), 
                      "NonZeroGenes"=length(nonzero.genes), "Gt1kUMIs"=sum(nonzero.genes > 1000)) 
write.table(filt.df, file=filt.file, sep="\t", quote=FALSE, row.names=FALSE) 
message("All done.") 
message(paste("Time taken to complete:", Sys.time() - start)) 