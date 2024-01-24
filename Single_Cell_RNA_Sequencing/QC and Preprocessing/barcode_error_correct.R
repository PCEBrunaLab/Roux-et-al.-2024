#! /usr/bin/env Rscript

## Perform error correction on barcodes
## This uses locality-sensitive hashing to reduce the number of pair-wise distance calculations that need to be made
## For error correction we only need to know _which_ barcodes are within an edit distance of `t` to do the grouping.

library(optparse)
library(stringdist)
library(Matrix)
library(igraph)
#library(MatrixExtra)


collapseBarcodes <- function(count.vec, threshold=0.5, parent=NULL){
    ## collapse down all barcodes to the most frequent (or parent if provided)
    ## only collapse if the counts are <= threshold of the parent/most frequent node
    
    if(is.null(parent)){
        parent.node <- names(count.vec)[which(count.vec == max(count.vec))]
    } else{
        parent.node <- parent
    }
    
    p.thresh <- ceiling(max(count.vec) * threshold)
    
    # need the mapping of barcode to collapsed barcode
    bc.map <- data.frame("Parent"=parent.node, "BC"=parent.node)
    rownames(bc.map) <- bc.map$BC 
    
    barcodes <- setdiff(names(count.vec), parent.node)
    for(x in seq_along(barcodes)){
        BC <- barcodes[x]
        if(count.vec[BC] <= p.thresh){
            bc.map <- rbind.data.frame(bc.map, data.frame("Parent"=parent.node, "BC"=BC))
        } else{
            bc.map <- rbind.data.frame(bc.map, data.frame("Parent"=BC, "BC"=BC))
        }
    }
    
    rownames(bc.map) <- NULL
    return(bc.map)
}


makeKmers <- function(k){
  # return all k^4 kmers using ATCG
  nts <- c("A", "T", "C", "G")
  kmer.grid <- expand.grid(nts, nts)

  for(x in seq_len(k-2)){
    kmer.grid <- expand.grid(apply(kmer.grid, 1, FUN=paste, collapse=""), nts)
  }
  kmer.grid <- unique(apply(kmer.grid, 1, FUN=paste, collapse=""))
  return(kmer.grid)
}


makeCharMatrix <- function(barcodes, kmers){
    # make a sparse characteristic matrix of the barcodes and kmers
    char.mat <- as(sapply(kmers, FUN=function(KMER) sapply(barcodes, FUN=function(BC) grepl(x=BC, pattern=KMER))) + 0, "dgCMatrix")
    dimnames(char.mat) <- list(barcodes, kmers)
    # remove Kmers that don't appear
    char.mat <- char.mat[, colSums(char.mat) > 0, drop=FALSE]

    return(char.mat)
}


sieveOfSundaram <- function(n){
    # The sieve of Sundaram is a simple deterministic algorithm for finding all the prime numbers up to a specified integer
    # strictly we should include 2, but we aren't interested in a small prime number here
    max.n <- ceiling((n+1)/2)
    L <- rep(TRUE, max.n)
    names(L) <- seq_len(max.n)
    
    for(i in seq_len(length(L))){
        for(j in seq(i, length(L))){
            if((i + j + (2*i*j)) <= max.n){
                L[i + j + (2*i*j)] <- FALSE
            }
        }
    }
    
    primes <- (2 * as.numeric(names(L[L]))) + 1
    
    return(max(primes[which(primes <= n)])) # just get the largest prime number
}


makeSigMatrix <- function(n.hash, barcodes, char.mat, prime.n){
  # make the hash signature matrix - output should be hashes by barcodes
  # running sample this many times on lots of barcodes is very slow
  # instead we need to generate N hash functions based on (a * x + 1) %% prime <- where prime is the nearest prime
  # number to the number of barcodes

  hash.facs <- sample(seq(1, 2*n.hash, 2)) # a series of integer scalars for the hash functions <- this is the first N odd numbers
  hash.perms <- sapply(seq_len(n.hash), FUN=function(HASH){
				                      (((hash.facs[HASH] * seq_len(ncol(char.mat))) + 1) %% prime.n) + 1
                                            }) # random hash functions
  
  #hash.perms <- sapply(seq_len(n.hash), FUN=function(HASH) sample(colnames(char.mat))) # permute the k-shingles

  sig.mat <- t(sapply(seq_len(n.hash), FUN=function(HASH) {	
    sapply(seq_len(nrow(char.mat)), FUN=function(BCS){
        ifelse(is.infinite(min(which(char.mat[BCS, hash.perms[, HASH]] == 1))), Inf, min(which(char.mat[BCS, hash.perms[, HASH]] == 1)))
        })
    }))

  dimnames(sig.mat) <- list(seq_len(n.hash), rownames(char.mat))
  return(sig.mat)
}


hashToBuckets <- function(r, b, n.hash, sig.mat, seed=42, w=1){
    # As we have an r-length vector for each column of the signature matrix, can we use a dot-product with the 
    # same random vector to map onto a scalar value which can then be assigned to a bucket? How about a random
    # gaussian vector?
    set.seed(42)
    rand.vec <- abs(rnorm(r, mean=0, sd=1))
    bands <- split(seq_len(n.hash), ceiling(seq_len(n.hash)/r))
    hash.buckets <- list()

    for(j in seq_along(bands)){
        j.band <- bands[[j]]
    	j.start <- min(j.band)
    	j.end <- max(j.band)
    	j.hashes <- apply(sig.mat[c(j.start:j.end), ], 2, FUN=function(HX){
            if(!all(is.infinite(HX))){
                floor((((HX %*% rand.vec)[, 1]) + runif(1, 0, w))/w)
            } else{
            NA
            }
	}) # hash each signature column
        names(j.hashes) <- colnames(sig.mat)
    	j.hashes <- j.hashes[!is.na(j.hashes)] # need to remove Infinite and NA sets
	j.hashes <- j.hashes[!is.infinite(j.hashes)] # need to remove Infinite and NA sets
    
        
        j.buckets <- seq_len(max(j.hashes) + 1) - 1
    	j.hash.list <- list()
    
	for(k in seq_along(j.buckets)){
            j.hash.list[[paste0(k)]] <- sig.mat[c(j.start:j.end), names(j.hashes[which(j.hashes == (j.buckets[k]+1))]), drop=FALSE]
        }
    
	hash.buckets[[paste0(j)]] <- j.hash.list
    }

    return(hash.buckets)
}


findShinglePairs <- function(hash.buckets, r){
    # check if a pair of shingles is similar in multiple band hashes, i.e. in the same bucket in multiple bands
    cand.pair.shingles <- c()

    for(j in seq_along(hash.buckets)){
        j.hash.buck <- hash.buckets[[j]]
	for(q in seq_along(j.hash.buck)){
            if(ncol(j.hash.buck[[q]]) > 1){
                # do these agree across 90% of entries?
            	comp.grid <- expand.grid(colnames(j.hash.buck[[q]]), colnames(j.hash.buck[[q]]))
            	comp.grid <- comp.grid[!apply(comp.grid, 1, FUN=function(PX) any(duplicated(PX))), ]
            	n.half <- nrow(comp.grid)/2
            	comp.grid <- comp.grid[c(1:n.half), ] # don't need redundant rows
            
		is.pair <- apply(comp.grid, 1, FUN=function(SIG){
		    sum((j.hash.buck[[q]])[, SIG[1]] == (j.hash.buck[[q]])[, SIG[2]]) > 0
                })
                if(sum(is.pair) > 0){
                    q.pairs <- apply(comp.grid[is.pair, ], 1, paste, collapse="_")
                    cand.pair.shingles <- unique(c(cand.pair.shingles, q.pairs))
            	}
            }
        }
    }
    
    return(unique(cand.pair.shingles))
}


computeEditDistances <- function(barcodes, cand.pair.shingles, char.mat, thresh){
    ## loop over the candidate pairs of k-shingles
    ## extract the barcodes with a 1 in these k-shingles and compute the edit distance between them?
    edit.mat <- as(matrix(0L, ncol=length(barcodes), nrow=length(barcodes)), "dgCMatrix")
    dimnames(edit.mat) <- list(barcodes, barcodes)
    
    for(k in seq_along(cand.pair.shingles)){
        k.pair <- unlist(strsplit(cand.pair.shingles[k], split="_", fixed=TRUE))
        k.mat <- char.mat[k.pair, , drop=FALSE]
	k.mat <- k.mat[rowSums(k.mat) > 0, , drop=FALSE]
    	k.dist <- stringdistmatrix(rownames(k.mat), rownames(k.mat), method="hamming", useBytes=TRUE)
    	dimnames(k.dist) <- list(rownames(k.mat), rownames(k.mat))
    	edit.mat[rownames(k.dist), rownames(k.dist)] <- (k.dist <= thresh) + 0
    }

    return(edit.mat)
}
   

parser <- OptionParser()
parser <- add_option(parser, c("-f", "--file"), type="character",
                     help="File containing barcode frequencies for a single sample")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output of file of corrected barcode frequencies")

opt <- parse_args(parser)

message("Reading barcode file ", opt$file)
bcs.df <- read.table(opt$file, sep="\t", header=TRUE, stringsAsFactors=FALSE)
   
message("Error correcting 14nt barcodes")
correct.14nt <- list()
nts <- c("A", "T", "C", "G")

### LSH isn't good with short barcodes - they are short enough and small enough number to compute distances deterministically though

i.bc14 <- unique(bcs.df$BC.14)
message("Found ", length(i.bc14), " barcodes")
i.bc14.counts <- sapply(i.bc14, FUN=function(BX) sum(bcs.df[bcs.df$BC.14 %in% BX, ]$count))

# stringdistmatrix fails for large samples size
message("Computing distances between ", length(i.bc14), " 14nt barcodes")
i.14.dist <- stringdistmatrix(i.bc14, i.bc14, method="hamming", useBytes=TRUE)
dimnames(i.14.dist) <- list(i.bc14, i.bc14)

# make the adj matrix by cutting > threshold
i.14.adj <- as((i.14.dist < 1.4) + 0, "dgCMatrix")

message("Constructing distinct graphs")
i.14.graphs <- decompose.graph(graph_from_adjacency_matrix(i.14.adj)) # pull out the disjoint graphs into a list of graphs

message("Collapsing similar barcode counts")
i.map <- do.call(rbind.data.frame, lapply(i.14.graphs,
                                          FUN=function(IX){
                                              collapseBarcodes(i.bc14.counts[names(i.bc14.counts) %in% names(V(IX))])
                                             }))

bc14.map.df <- merge(bcs.df, i.map, by.y=c("BC"), by.x=c("BC.14"))
rownames(bc14.map.df) <- NULL

message("Error correcting 30nt barcodes")
i.bc30 <- unique(bcs.df$BC.30)
message("Found ", length(i.bc30), " barcodes")
i.bc30.counts <- sapply(i.bc30, FUN=function(BX) sum(bcs.df[bcs.df$BC.30 %in% BX, ]$count))

# stringdistmatrix fails for large samples size
message("Computing distances between ", length(i.bc30), " 30nt barcodes")
i.30.dist <- stringdistmatrix(i.bc30, i.bc30, method="hamming", useBytes=TRUE)
dimnames(i.30.dist) <- list(i.bc30, i.bc30)

# make the adj matrix by cutting > threshold
i.30.adj <- as((i.30.dist < 3) + 0, "dgCMatrix")

message("Constructing distinct graphs")
i.30.graphs <- decompose.graph(graph_from_adjacency_matrix(i.30.adj)) # pull out the disjoint graphs into a list of graphs

message("Collapsing similar barcode counts")
i.map <- do.call(rbind.data.frame, lapply(i.30.graphs,
                                          FUN=function(IX){
                                              collapseBarcodes(i.bc30.counts[names(i.bc30.counts) %in% names(V(IX))])
                                             }))

bc30.map.df <- merge(bcs.df, i.map, by.y=c("BC"), by.x=c("BC.30"))
rownames(bc30.map.df) <- NULL

ofile.14nt <- paste0(opt$output, "_14ntBarcodes.txt")
ofile.30nt <- paste0(opt$output, "_30ntBarcodes.txt")

message("Writing error corrected 14nt barcodes to ", ofile.14nt)
write.table(bc14.map.df, file=ofile.14nt, sep="\t", row.names=FALSE, quote=FALSE)

message("Writing error corrected 30nt barcodes to ", ofile.30nt)
write.table(bc30.map.df, file=ofile.30nt, sep="\t", row.names=FALSE, quote=FALSE)

message("Error correcting cell barcodes")
i.bccb <- unique(bcs.df$CB)
message("Found ", length(i.bccb), " cell barcodes")
i.bccb.counts <- sapply(i.bccb, FUN=function(BX) sum(bcs.df[bcs.df$CB %in% BX, ]$count))

# stringdistmatrix fails for large samples size
message("Computing distances between ", length(i.bccb), " cell barcodes")
# this needs to be slightly different

i.cb.dist <- stringdistmatrix(i.bccb, i.bccb, method="hamming", useBytes=TRUE)
dimnames(i.cb.dist) <- list(i.bccb, i.bccb)

# make the adj matrix by cutting > threshold
i.cb.adj <- as((i.cb.dist < 1) + 0, "dgCMatrix")

message("Constructing distinct graphs")
i.cb.graphs <- decompose.graph(graph_from_adjacency_matrix(i.cb.adj)) # pull out the disjoint graphs into a list of graphs

message("Collapsing similar barcode counts") # the counts for any given barcode should be fairly small no?
i.map <- do.call(rbind.data.frame, lapply(i.cb.graphs,
                                          FUN=function(IX){
                                              collapseBarcodes(i.bccb.counts[names(i.bccb.counts) %in% names(V(IX))])
                                             }))

bccb.map.df <- merge(bcs.df, i.map, by.y=c("BC"), by.x=c("CB"))
rownames(bccb.map.df) <- NULL

ofile.cbnt <- paste0(opt$output, "_CellBarcodes.txt")
message("Writing error corrected cell barcodes to ", ofile.cbnt)
write.table(bccb.map.df, file=ofile.cbnt, sep="\t", row.names=FALSE, quote=FALSE)

print(head(bccb.map.df))

## error correct UMIs _within_ each cell barcode?
message("Error correcting UMIs within each cell barcode")
new.unique.cb <- unique(bccb.map.df$CB)

i.bccb <- unique(bcs.df$CB)
message("Found ", length(i.bccb), " corrected cell barcodes")

umi.corr.list <- list()
for(x in seq_along(i.bccb)){
    x.bcsdf <- bccb.map.df[bccb.map.df$CB %in% i.bccb[x], ]
    x.umi <- unique(x.bcsdf$UMI)
    message("Found ", length(x.umi), " UMIs")
    
    i.umi.counts <- sapply(x.umi, FUN=function(BX) sum(x.bcsdf[x.bcsdf$UMI %in% BX, ]$count))

    # stringdistmatrix fails for large samples size
    message("Computing distances between ", length(x.umi), " UMIs")
    i.umi.dist <- stringdistmatrix(x.umi, x.umi, method="hamming", useBytes=TRUE)
    dimnames(i.umi.dist) <- list(x.umi, x.umi)

    # make the adj matrix by cutting > threshold
    i.umi.adj <- as((i.umi.dist < 3) + 0, "dgCMatrix")

    message("Constructing distinct graphs")
    i.umi.graphs <- decompose.graph(graph_from_adjacency_matrix(i.umi.adj)) # pull out the disjoint graphs into a list of graphs

    message("Collapsing similar barcode counts")
    i.map <- do.call(rbind.data.frame, lapply(i.umi.graphs,
                                              FUN=function(IX){
                                                  collapseBarcodes(i.umi.counts[names(i.umi.counts) %in% names(V(IX))])
                                                 }))

    umi.map.df <- merge(x.bcsdf, i.map, by.y=c("BC"), by.x=c("UMI"))
    rownames(umi.map.df) <- NULL
    umi.map.df$CB <- i.bccb[x]
    umi.corr.list[[i.bccb[x]]] <- umi.map.df
}

umi.all.df <- do.call(rbind.data.frame, umi.corr.list)
ofile.umi <- paste0(opt$output, "_UMIs.txt")
message("Writing error corrected UMIs to ", ofile.umi)
write.table(umi.all.df, file=ofile.umi, sep="\t", row.names=FALSE, quote=FALSE)

message("All done")
