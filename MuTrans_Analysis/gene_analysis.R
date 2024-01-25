## Gene set enrichment analysis (GSEA) on MuTrans data

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)
library(reticulate)

## [ Anndata to Seurat ] ----

## Load python packages
use_condaenv("anndata")
sc <- import("scanpy")
hdf5 <- import("hdf5plugin")

## Untreated data
ut.adata <- sc$read_h5ad("ut_mutrans.h5ad")

ut.counts <- t(as.matrix(ut.adata$layers["logcounts"]))
colnames(ut.counts) <- ut.adata$obs_names$to_list()
rownames(ut.counts) <- ut.adata$var_names$to_list()
ut.counts <- Matrix::Matrix(ut.counts, sparse = TRUE)

ut.seurat <- CreateSeuratObject(ut.counts)
ut.seurat <- AddMetaData(ut.seurat, ut.adata$obs)

saveRDS(ut.seurat, "ut_seurat_mutrans.rds")

## Cisplatin data
cis.adata <- sc$read_h5ad("cis_mutrans.h5ad")

cis.counts <- t(as.matrix(cis.adata$layers["logcounts"]))
colnames(cis.counts) <- cis.adata$obs_names$to_list()
rownames(cis.counts) <- cis.adata$var_names$to_list()
cis.counts <- Matrix::Matrix(cis.counts, sparse = TRUE)

cis.seurat <- CreateSeuratObject(cis.counts)
cis.seurat <- AddMetaData(cis.seurat, cis.adata$obs)

saveRDS(cis.seurat, "cis_seurat_mutrans.rds")

## JQ1 data
jq1.adata <- sc$read_h5ad("jq1_mutrans.h5ad")

jq1.counts <- t(as.matrix(jq1.adata$layers["logcounts"]))
colnames(jq1.counts) <- jq1.adata$obs_names$to_list()
rownames(jq1.counts) <- jq1.adata$var_names$to_list()
jq1.counts <- Matrix::Matrix(jq1.counts, sparse = TRUE)

jq1.seurat <- CreateSeuratObject(jq1.counts)
jq1.seurat <- AddMetaData(jq1.seurat, jq1.adata$obs)

saveRDS(jq1.seurat, "jq1_seurat_mutrans.rds")

## [ Markers ] ----

## Untreated data
ut.seurat <- readRDS("ut_seurat_mutrans.rds")

Idents(ut.seurat) <- ut.seurat$attractor

ut.wilcox.markers <- FindAllMarkers(ut.seurat, assay = "RNA", slot = "counts",
                                    only.pos = TRUE, densify = TRUE, min.pct = 0.05,
                                    logfc.threshold = 0.2, test.use = "wilcox")

ut.wilcox.list <- lapply(unique(ut.wilcox.markers$cluster), function(x){
  tmp.df <- ut.wilcox.markers[ut.wilcox.markers$cluster==x & ut.wilcox.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})

names(ut.wilcox.list) <- unique(ut.wilcox.markers$cluster)

require(clusterProfiler)
require(org.Hs.eg.db)

## Gene Ontology biological processes
ut.wilcox.GO.BP.list <- lapply(ut.wilcox.list, function(x){
  enrichGO(x$gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

require(msigdbr)

## C2 curated gene sets
C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])

ut.wilcox.C2.list <- lapply(ut.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

## Hallmark gene sets
H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

ut.wilcox.H.list <- lapply(ut.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

require(writexl)

ut.wilcox.lists <- list("ut.wilcox.GO.BP.list" = ut.wilcox.GO.BP.list,
                        "ut.wilcox.C2.list" = ut.wilcox.C2.list,
                        "ut.wilcox.H.list" = ut.wilcox.H.list,
                        "ut.wilcox.list" = ut.wilcox.list)

ut.wilcox.names <- sub("list", "sheets", names(ut.wilcox.lists))

for(i in 1:length(ut.wilcox.lists)){
  wilcox.list <- ut.wilcox.lists[[i]]
  list.sheets <- lapply(wilcox.list, as.data.frame)
  assign(ut.wilcox.names[[i]], list.sheets)
}

write_xlsx(ut.wilcox.GO.BP.sheets, "ut_markers_BP_ontology.xlsx")
write_xlsx(ut.wilcox.C2.sheets, "ut_markers_C2_geneset.xlsx")
write_xlsx(ut.wilcox.H.sheets, "ut_markers_Hallmarks.xlsx")
write_xlsx(ut.wilcox.sheets, "ut_markers_wilcox.xlsx")

## cisplatin data
cis.seurat <- readRDS("cis_seurat_mutrans.rds")

Idents(cis.seurat) <- cis.seurat$attractor

cis.wilcox.markers <- FindAllMarkers(cis.seurat, assay = "RNA", slot = "counts",
                                     only.pos = TRUE, densify = TRUE, min.pct = 0.05,
                                     logfc.threshold = 0.2, test.use = "wilcox")

cis.wilcox.list <- lapply(unique(cis.wilcox.markers$cluster), function(x){
  tmp.df <- cis.wilcox.markers[cis.wilcox.markers$cluster==x & cis.wilcox.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})

names(cis.wilcox.list) <- unique(cis.wilcox.markers$cluster)

## Gene Ontology biological processes
cis.wilcox.GO.BP.list <- lapply(cis.wilcox.list, function(x){
  enrichGO(x$gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

## C2 curated gene sets
C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])

cis.wilcox.C2.list <- lapply(cis.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

## Hallmark gene sets
H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

cis.wilcox.H.list <- lapply(cis.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

cis.wilcox.lists <- list("cis.wilcox.GO.BP.list" = cis.wilcox.GO.BP.list,
                         "cis.wilcox.C2.list" = cis.wilcox.C2.list,
                         "cis.wilcox.H.list" = cis.wilcox.H.list,
                         "cis.wilcox.list" = cis.wilcox.list)

cis.wilcox.names <- sub("list", "sheets", names(cis.wilcox.lists))

for(i in 1:length(cis.wilcox.lists)){
  wilcox.list <- cis.wilcox.lists[[i]]
  list.sheets <- lapply(wilcox.list, as.data.frame)
  assign(cis.wilcox.names[[i]], list.sheets)
}

write_xlsx(cis.wilcox.GO.BP.sheets, "cis_markers_BP_ontology.xlsx")
write_xlsx(cis.wilcox.C2.sheets, "cis_markers_C2_geneset.xlsx")
write_xlsx(cis.wilcox.H.sheets, "cis_markers_Hallmarks.xlsx")
write_xlsx(cis.wilcox.sheets, "cis_markers_wilcox.xlsx")

## jq1 data
jq1.seurat <- readRDS("jq1_seurat_mutrans.rds")

Idents(jq1.seurat) <- jq1.seurat$attractor

jq1.wilcox.markers <- FindAllMarkers(jq1.seurat, assay = "RNA", slot = "counts",
                                     only.pos = TRUE, densify = TRUE, min.pct = 0.05,
                                     logfc.threshold = 0.2, test.use = "wilcox")

jq1.wilcox.list <- lapply(unique(jq1.wilcox.markers$cluster), function(x){
  tmp.df <- jq1.wilcox.markers[jq1.wilcox.markers$cluster==x & jq1.wilcox.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})

names(jq1.wilcox.list) <- unique(jq1.wilcox.markers$cluster)

## Gene Ontology biological processes
jq1.wilcox.GO.BP.list <- lapply(jq1.wilcox.list, function(x){
  enrichGO(x$gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

## C2 curated gene sets
C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])

jq1.wilcox.C2.list <- lapply(jq1.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

## Hallmark gene sets
H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

jq1.wilcox.H.list <- lapply(jq1.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

jq1.wilcox.lists <- list("jq1.wilcox.GO.BP.list" = jq1.wilcox.GO.BP.list,
                         "jq1.wilcox.C2.list" = jq1.wilcox.C2.list,
                         "jq1.wilcox.H.list" = jq1.wilcox.H.list,
                         "jq1.wilcox.list" = jq1.wilcox.list)

jq1.wilcox.names <- sub("list", "sheets", names(jq1.wilcox.lists))

for(i in 1:length(jq1.wilcox.lists)){
  wilcox.list <- jq1.wilcox.lists[[i]]
  list.sheets <- lapply(wilcox.list, as.data.frame)
  assign(jq1.wilcox.names[[i]], list.sheets)
}

write_xlsx(jq1.wilcox.GO.BP.sheets, "jq1_markers_BP_ontology.xlsx")
write_xlsx(jq1.wilcox.C2.sheets, "jq1_markers_C2_geneset.xlsx")
write_xlsx(jq1.wilcox.H.sheets, "jq1_markers_Hallmarks.xlsx")
write_xlsx(jq1.wilcox.sheets, "jq1_markers_wilcox.xlsx")
