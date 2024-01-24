#! /usr/bin/env Rscript

# Use the joint distributions of counts over all of the CMOs to assign cells
# to particular samples


required.packages <- c("devtools", "SingleCellExperiment", "DropletUtils", "scran", "Matrix", 
                       "MatrixGenerics", "BiocParallel", "biomaRt", "optparse", "scater", "Cairo", "fitdistrplus")

if (!all(required.packages %in% installed.packages()[,"Package"])){
  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]])
}

#load libraries
library(SingleCellExperiment)
library(DropletUtils)
library(fitdistrplus)
library(scran)
library(scater)
library(Matrix)
library(MatrixGenerics)
library(BiocParallel)
library(optparse)

start <- Sys.time()

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the SCE object file")

parser <- add_option(parser, c("-n", "--ncores"), type="numeric",
                     help="The number of cores to use for parallelised steps")

parser <- add_option(parser, c("-q", "--quantile"), type="numeric",
                     help="The quantile threshold to define the true signal of CMO tag counts")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output CMO assignment for each single cell")

opt <- parse_args(parser)


mcparam <- MulticoreParam(workers=opt$ncores)
register(mcparam)

message(paste0("Reading counts for sample: ", opt$id))
in.sce <- readRDS(opt$SCE)

if(is.null(colnames(in.sce))){
  message("Setting barcode names on cells")
  colnames(in.sce) <- colData(in.sce)$Barcode
}

cmo.tags <- rownames(in.sce)[rowData(in.sce)[, "Type"] == "Multiplexing Capture"]
message("Subsetting to ", length(cmo.tags), " CMOs")
in.sce <- in.sce[cmo.tags, ]

used.cmo <- rowSums2(counts(in.sce)) > floor(ncol(in.sce)/6)
print(rowSums2(counts(in.sce)))
message("Dropping ", sum(!used.cmo), " unused CMOs")
in.sce <- in.sce[used.cmo, ]

# remove cells with insufficient CMO tag coverage
keep.cells <- colSums2(counts(in.sce)) > 0

message("Dropping ", sum(!keep.cells), " cells with poor CMO coverage")
in.sce <- in.sce[, keep.cells]
print(in.sce)

message("Normalising CMO counts")
norm_matrix <- function(X){
  as(apply(X, 2, function(P) P/sum(P)), "dgCMatrix")
}

# CPM normalise the HTO counts
norm.cmo <- norm_matrix(counts(in.sce))
cmo.tags <- rownames(in.sce)

print(sum(is.na(norm.cmo)))

# define the expected number of clusters - this should be based on the number of unique tags
exp_clusters <- function(C){
  n <- length(C)
  combos <- nrow(expand.grid(c(1:n), c(1:n)))
  n + ((combos - n)/2)
}

exp.k <- exp_clusters(rownames(in.sce))
print(exp.k)

# do a k-means cluster on the cells in feature space. Need to define a threshold for negative cells.
cell_class <- list()

exp.k <- length(cmo.tags)
message("Performing k-means clustering with ", exp.k, " clusters")
k.means <- kmeans(x=t(as.matrix(norm.cmo)), centers=exp.k)
# fit a negative binomial model to the remaining unnormalised HTO counts for each HTO barcode
quantile_class <- list()

message("Assignming CMOs: ", length(cmo.tags), " uniqe CMOs found")
for(x in seq_along(cmo.tags)){
    x.tag <- cmo.tags[x]

    xtag.counts <- counts(in.sce)[x.tag, ]
    print(x.tag)
    print(head(xtag.counts))

    # fit a background negative binomial, excluding the cells for a given cluster with the highest
    # average counts for that HTO tag

    max.center <- max(k.means$centers[, x.tag])
    max.clust <- which(k.means$centers[, x.tag] == max.center)
    cmo.cells <- k.means$cluster == max.clust

    message("Fitting background negative binomial distribution for ", x.tag)
    # exclude top 0.5% of cells as there are +ve cells that are mis-classified by the crude k-means
    not.cluster.counts <- xtag.counts[!cmo.cells]
    top.05.counts <- not.cluster.counts[order(not.cluster.counts, decreasing=TRUE)]
    top.05 <- floor(length(top.05.counts) * 0.005)
    
    x.negbin <- fitdist(top.05.counts[top.05:length(top.05.counts)], "nbinom")
    x.99.q <- as.numeric(quantile(x.negbin, probs=opt$quantile)$quantiles)

    # classify each cells that is >= the 99th quantile as "positive"
    pos.cells <- as.numeric(xtag.counts >= x.99.q)
    names(pos.cells) <- names(xtag.counts)
    message("Found " , sum(pos.cells), " cells in the ", opt$quantile, " quantile of counts distribution")
    quantile_class[[x.tag]] <- pos.cells
    print(sum(pos.cells))
 }

# classify each cell into the relevant HTO
# if more than one then call it a multiplet

message("Calling cell - CMO assignments")
cmo.class <- do.call(cbind.data.frame,
                     quantile_class)
print(head(cmo.class))
cmo.sums <- rowSums(cmo.class)
print(head(names(cmo.sums)))
singlets <- names(cmo.sums)[cmo.sums == 1]
multiplets <- names(cmo.sums)[cmo.sums > 1]
dropout <- names(cmo.sums)[cmo.sums == 0]
single.class <- apply(cmo.class[singlets, ], 1, FUN=function(C) colnames(cmo.class)[C == 1])
multi.class <- apply(cmo.class[multiplets, ], 1, FUN=function(Q) paste(colnames(cmo.class)[Q == 1], collapse="."))
  
all.cell.class <- data.frame("ID"=c(singlets, multiplets, dropout),
	       	  	     "Class"=c(rep("Singlet", length(singlets)),
			               rep("Multiplet", length(multiplets)),
				       rep("Dropout", length(dropout))),
                             "CMO"=c(single.class[singlets], multi.class[multiplets], rep("None", length(dropout))))

message("Writing CMO assigment to ", opt$output)
write.table(all.cell.class,
            file=opt$output,
            sep="\t", row.names=FALSE, quote=FALSE)

print( Sys.time() - start )


