#Compute Nearest Gene Per Peaks

#Initial R Env with correct installation of packages - ArchR_Renv
renv::activate()

#Set Working Directory and Arch R parameters
setwd("~/Desktop/sn_clones/")

addArchRThreads(threads = 4) 
addArchRGenome("hg38")

#Load in data files
project <- readRDS("project_clones_working/Save-ArchR-Project.rds")
peak_set <- getPeakSet(project)

#Get Gene Annotations
gene_annot <- getGenes(project)
# separate in + and - strand to define TSS
# TSS are the start coordinate for genes on positive strand and end coordinate for negative strand
gene_annot_TSS <- gene_annot
end(gene_annot_TSS[which(strand(gene_annot_TSS) == "+")]) <- start(gene_annot_TSS[which(strand(gene_annot_TSS) == "+")]) 
start(gene_annot_TSS[which(strand(gene_annot_TSS) == "-")]) <- end(gene_annot_TSS[which(strand(gene_annot_TSS) == "-")]) 

#Nearest gene for each peak
nearest_gene <- as.data.frame(GenomicRanges::nearest(peak_set, gene_annot, select = "all", ignore.strand = FALSE))
colnames(nearest_gene) <- c("peak", "gene") # replace column names
nearest_gene$peak <- paste0(seqnames(peak_set), "_", start(peak_set), "_", end(peak_set))[nearest_gene$peak] # replace by peak name
nearest_gene$gene <- gene_annot$symbol[nearest_gene$gene] # replace by gene name
nearest_gene2 <- nearest_gene %>% group_by(peak) %>% summarise(gene = paste0(gene, collapse=", ")) %>% as.data.frame() # collapse when several genes associated to peak
saveRDS(nearest_gene2, "datafiles/nearest_gene_per_peak.RDS")

#Distance to nearest TSS for each peak
dist_TSS <- as.data.frame(GenomicRanges::distanceToNearest(peak_set, gene_annot_TSS, ignore.strand = TRUE)) # ignore strand because we work with 1 bp for gene TSS
colnames(dist_TSS) <- c("peak", "gene", "dist") # replace column names
dist_TSS$peak <- paste0(seqnames(peak_set), "_", start(peak_set), "_", end(peak_set))[dist_TSS$peak] # replace by peak name
dist_TSS$gene <- gene_annot$symbol[dist_TSS$gene] # replace by gene name
saveRDS(dist_TSS, "datafiles/dist_nearest_gene_per_peak.RDS")

