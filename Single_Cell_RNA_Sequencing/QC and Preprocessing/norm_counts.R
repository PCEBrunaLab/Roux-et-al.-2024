#! /usr/bin/env Rscript

## Given an input dataset of multiple batches - normalize across batches

library(DropletUtils)
library(SingleCellExperiment)
library(scran)
library(scater)
library(batchelor)
library(irlba)
library(optparse)
library(Matrix)
library(umap)
library(BiocParallel)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="A comma-separated list of SCE file paths")

parser <- add_option(parser, c("-c", "--calls"), type="character",
                     help="A comma-separated list of doublet calls")

parser <- add_option(parser, c("-b", "--breaks"), action="store_true", default=TRUE,
                     help="If set then the data will be broken into arbitrary chunks to speed up size factor estimation")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files denoised data and combined SCE object")

opt <- parse_args(parser)

sce.file.list <- unlist(strsplit(opt$SCE, split=",", fixed=TRUE))
print(sce.file.list)
sce.names <- gsub(sce.file.list, pattern="(\\S+)/(Neuro_[A-Z])(edited_SCE\\.RDS)", replacement="\\2")
print(sce.names)
names(sce.file.list) <- sce.names

sce.list <- list()
for(x in seq_along(sce.names)){
  x.nom <- sce.names[x]
  x.sce <- readRDS(sce.file.list[[x.nom]])
  colnames(x.sce) <- paste(colData(x.sce)$Sample, colData(x.sce)$Barcode, sep="_")
  sce.list[[x.nom]] <- x.sce
}

message("Combining ", length(sce.list), " SCE files together")
big.sce <- do.call(cbind, sce.list)

call.file.list <- unlist(strsplit(opt$calls, split=",", fixed=TRUE))
call.names <- gsub(call.file.list, pattern="(\\S+)/(Neuro_[A-Z])(_demuxed\\.tsv)", replacement="\\2")
names(call.file.list) <- call.names

call.list <- list()
for(x in seq_along(call.names)){
  x.name <- call.names[x]
  x.call <- read.table(call.file.list[[x.name]], sep="\t", header=TRUE)
  x.call$ID <- paste(x.name, x.call$ID, sep="_")
  call.list[[x.name]] <- x.call
}

call.df <- do.call(rbind.data.frame, call.list)
rownames(call.df) <- call.df$ID

# take the intersection of singlets and called barcodes and make the SCE object
keep.bcs <- intersect(call.df$ID[call.df$Class %in% c("Singlet")], colnames(big.sce))

message("Keeping ", length(keep.bcs), " singlet cells")
big.sce <- big.sce[, colnames(big.sce) %in% keep.bcs]
counts(big.sce) <- as(counts(big.sce), "dgCMatrix")

n.cells <- ncol(big.sce)
n.genes <- nrow(big.sce)
message(paste0("Computing gene expression sparsity over ", n.cells, " droplets."))
#gene_sparsity <- rowSums(counts(big.sce) < 1)/n.cells
#keep_genes <- gene_sparsity < 0.99
keep_genes <- rowMeans(counts(big.sce)) > 0.01
genes <- rownames(big.sce)

message(paste0("Using ", sum(keep_genes), " genes for size factor estimation"))

if(isFALSE(opt$breaks)){
    # try to find the smallest integer that leads to whole number division
    poss.ints <- c(3:50)
    mods <- ncol(big.sce) %% poss.ints
    if(sum(mods == 0) < 1){
      warning("No integer values <50 that creates an integer devisor - setting to 5 chunks")
      n.chunks <- 5
    } else{
      n.chunks <- min(poss.ints[which(mods == 0)])
    }
    
    max.n <- floor(ncol(big.sce)/n.chunks)
    d <- seq_along(c(1:ncol(big.sce)))
    max.points <- unlist(lapply(split(d, ceiling(d/max.n)), max))
    sce.list <- list()
    for(q in seq_along(max.points)){
        q.max <- max.points[q]
	if(q == 1){
            x.sce <- big.sce[, 1:q.max]
	} else{
	x.sce <- big.sce[, (old.max+1):q.max]
	}
	old.max <- q.max

	# set cluster size to 5% of data
       	cluster.size <- ceiling(ncol(x.sce) * 0.1)
        message(paste0("Cluster size set to ", cluster.size))
	print(class(counts(x.sce)))
        clusters <- quickCluster(x.sce, min.size=cluster.size,
                                 subset.row=keep_genes,
                                 method="igraph")
        max.size <- floor(cluster.size/2)

        # change the window size in 50% increments
	size.inc <- ceiling(max.size * 0.5)
	message(paste0("Estimating size factors using ", size.inc, " cell increments"))
    	# how are there so many negative size factors?
    	x.sce <- computeSumFactors(x.sce,
                                   max.cluster.size=max.size,
                                   positive=TRUE,
			           subset.row=keep_genes,
				   BPPARAM=MulticoreParam(workers=2),
                                   assay.type='counts', clusters=clusters)

        neg.sf <- sum(sizeFactors(x.sce) < 0)
	if(neg.sf > 0){
           message(paste0(neg.sf, " negative size factors estimated - consider a higher sparsity threshold"))
        }

        message("Normalising single-cell expression values")
        x.sce <- logNormCounts(x.sce)
	colData(x.sce)$Batch <- as.character(q)
	saveRDS(x.sce, file=gsub(opt$output, pattern="\\.RDS", replacement=paste0("Batch", q, ".RDS")))
	sce.list[[paste0(q)]] <- x.sce
	sink(file="/dev/null")
	rm(list=c("x.sce", "clusters"))
	gc()
	sink(file=NULL)
    }
    big.sce <- do.call(cbind, sce.list)

    message(paste0("Performing multi-batch normalisation across ", length(unique(big.sce$Batch)), " batches"))
    big.sce <- multiBatchNorm(big.sce, batch=big.sce$Batch, assay.type="counts",
	                      min.mean=0.01, normalize.all=TRUE, preserve.single=TRUE)
    } else{
    # set cluster size to 5% of data
    cluster.size <- ceiling(ncol(big.sce) * 0.01)
    message(paste0("Cluster size set to ", cluster.size))

    clusters <- quickCluster(big.sce, min.size=cluster.size,
                             subset.row=keep_genes,
         		     BPPARAM=MulticoreParam(workers=12),
                             method="igraph")
    max.size <- floor(cluster.size/2)

    # change the window size in 50% increments
    size.inc <- ceiling(max.size * 0.5)
    message(paste0("Estimating size factors using ", size.inc, " cell increments"))
    # how are there so many negative size factors?
    big.sce <- computeSumFactors(big.sce,
                             max.cluster.size=max.size,
                             positive=TRUE,
			     subset.row=keep_genes,
                             #BPPARAM=MulticoreParam(workers=12),
                             assay.type='counts', clusters=clusters)
  
    # the scater function is normalise, not normalize
    # I appreciate the British spelling, but this might cause a clash
    # with igraph for those not familiar
    # count the number of negative size factors
    neg.sf <- sum(sizeFactors(big.sce) < 0)
    if(neg.sf > 0){
        message(paste0(neg.sf, " negative size factors estimated - consider a higher sparsity threshold"))
    }

    message("Normalising single-cell expression values")
    big.sce <- logNormCounts(big.sce)
}

big.sce <- big.sce[!is.na(rownames(rowData(big.sce))), ]
#rownames(big.sce) <- rowData(big.sce)$ensembl_gene_id
rownames(big.sce) <- genes

saveRDS(big.sce, file=opt$output)
message("All done")