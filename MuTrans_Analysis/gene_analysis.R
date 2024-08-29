## Gene analysis on MuTrans data

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)
library(reticulate)

setwd("C:/Users/echen/OneDrive - The Institute of Cancer Research/Documents/nb sc analysis 2 updated/mutrans")

## [ Anndata to Seurat ] ----

## Load python packages
use_condaenv("anndata")
sc <- import("scanpy")
hdf5 <- import("hdf5plugin")

## Untreated data
ut.adata <- sc$read_h5ad("MuTrans-release/Example/data/ut_mutrans.h5ad")

ut.counts <- t(as.matrix(ut.adata$layers["raw_counts"]))
colnames(ut.counts) <- ut.adata$obs_names$to_list()
rownames(ut.counts) <- ut.adata$var_names$to_list()
ut.counts <- Matrix::Matrix(ut.counts, sparse = TRUE)

ut.seurat <- CreateSeuratObject(ut.counts)
ut.seurat <- AddMetaData(ut.seurat, ut.adata$obs)

saveRDS(ut.seurat, "../data/ut_seurat_mutrans.rds")

## Cisplatin entry data
entry.adata <- sc$read_h5ad("MuTrans-release/Example/data/entry_mutrans.h5ad")

entry.counts <- t(as.matrix(entry.adata$layers["raw_counts"]))
colnames(entry.counts) <- entry.adata$obs_names$to_list()
rownames(entry.counts) <- entry.adata$var_names$to_list()
entry.counts <- Matrix::Matrix(entry.counts, sparse = TRUE)

entry.seurat <- CreateSeuratObject(entry.counts)
entry.seurat <- AddMetaData(entry.seurat, entry.adata$obs)

saveRDS(entry.seurat, "../data/entry_seurat_mutrans.rds")

## Cisplatin exit data (cisplatin -> recovery 1 week -> recovery 4 weeks)
exitv2.adata <- sc$read_h5ad("MuTrans-release/Example/data/exit_v2_mutrans.h5ad")

exitv2.counts <- t(as.matrix(exitv2.adata$layers["raw_counts"]))
colnames(exitv2.counts) <- exitv2.adata$obs_names$to_list()
rownames(exitv2.counts) <- exitv2.adata$var_names$to_list()
exitv2.counts <- Matrix::Matrix(exitv2.counts, sparse = TRUE)

exitv2.seurat <- CreateSeuratObject(exitv2.counts)
exitv2.seurat <- AddMetaData(exitv2.seurat, exitv2.adata$obs)

saveRDS(exitv2.seurat, "../data/exit_v2_seurat_mutrans.rds")

## [ UT markers ] ----

ut.seurat <- readRDS("../data/ut_seurat_mutrans.rds")

Idents(ut.seurat) <- ut.seurat$attractor

ut.seurat <- NormalizeData(ut.seurat)
ut.seurat <- FindVariableFeatures(ut.seurat, nfeatures = 1000)
ut.seurat <- ScaleData(ut.seurat)

ut.wilcox.markers <- FindAllMarkers(ut.seurat, assay = "RNA", slot = "data",
                                    only.pos = TRUE, densify = TRUE, min.pct = 0.05,
                                    logfc.threshold = 0.2, test.use = "wilcox")

dir.create("markers", recursive = TRUE)

saveRDS(ut.wilcox.markers, "markers/ut_markers_wilcox.rds")

ut.wilcox.list <- lapply(unique(ut.wilcox.markers$cluster), function(x){
  tmp.df <- ut.wilcox.markers[ut.wilcox.markers$cluster==x & ut.wilcox.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})

names(ut.wilcox.list) <- unique(ut.wilcox.markers$cluster)

saveRDS(ut.wilcox.list, "markers/ut_markers_wilcox_list.rds")

table(ut.wilcox.markers$cluster)
## 0    1    2    3    4 
## 777 1150  888  597  700 

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

saveRDS(ut.wilcox.GO.BP.list, "markers/ut_markers_BP_ontology.rds")
saveRDS(ut.wilcox.C2.list, "markers/ut_markers_C2_geneset.rds")
saveRDS(ut.wilcox.H.list, "markers/ut_markers_Hallmarks.rds")

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

write_xlsx(ut.wilcox.GO.BP.sheets, "markers/ut_markers_BP_ontology.xlsx")
write_xlsx(ut.wilcox.C2.sheets, "markers/ut_markers_C2_geneset.xlsx")
write_xlsx(ut.wilcox.H.sheets, "markers/ut_markers_Hallmarks.xlsx")
write_xlsx(ut.wilcox.sheets, "markers/ut_markers_wilcox.xlsx")

## [ Entry markers ] ----

entry.seurat <- readRDS("../data/entry_seurat_mutrans.rds")

Idents(entry.seurat) <- entry.seurat$attractor

entry.seurat <- NormalizeData(entry.seurat)
entry.seurat <- FindVariableFeatures(entry.seurat, nfeatures = 1000)
entry.seurat <- ScaleData(entry.seurat)

entry.wilcox.markers <- FindAllMarkers(entry.seurat, assay = "RNA", slot = "data",
                                       only.pos = TRUE, densify = TRUE, min.pct = 0.05,
                                       logfc.threshold = 0.2, test.use = "wilcox")

saveRDS(entry.wilcox.markers, "markers/entry_markers_wilcox.rds")

entry.wilcox.list <- lapply(unique(entry.wilcox.markers$cluster), function(x){
  tmp.df <- entry.wilcox.markers[entry.wilcox.markers$cluster==x & entry.wilcox.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})

names(entry.wilcox.list) <- unique(entry.wilcox.markers$cluster)

saveRDS(entry.wilcox.list, "markers/entry_markers_wilcox_list.rds")

table(entry.wilcox.markers$cluster)
## 0    1    2    3    4    5 
## 956  748  742  774  629 1051 

require(clusterProfiler)
require(org.Hs.eg.db)

## Gene Ontology biological processes
entry.wilcox.GO.BP.list <- lapply(entry.wilcox.list, function(x){
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

entry.wilcox.C2.list <- lapply(entry.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

## Hallmark gene sets
H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

entry.wilcox.H.list <- lapply(entry.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

saveRDS(entry.wilcox.GO.BP.list, "markers/entry_markers_BP_ontology.rds")
saveRDS(entry.wilcox.C2.list, "markers/entry_markers_C2_geneset.rds")
saveRDS(entry.wilcox.H.list, "markers/entry_markers_Hallmarks.rds")

require(writexl)

entry.wilcox.lists <- list("entry.wilcox.GO.BP.list" = entry.wilcox.GO.BP.list,
                           "entry.wilcox.C2.list" = entry.wilcox.C2.list,
                           "entry.wilcox.H.list" = entry.wilcox.H.list,
                           "entry.wilcox.list" = entry.wilcox.list)

entry.wilcox.names <- sub("list", "sheets", names(entry.wilcox.lists))

for(i in 1:length(entry.wilcox.lists)){
  wilcox.list <- entry.wilcox.lists[[i]]
  list.sheets <- lapply(wilcox.list, as.data.frame)
  assign(entry.wilcox.names[[i]], list.sheets)
}

write_xlsx(entry.wilcox.GO.BP.sheets, "markers/entry_markers_BP_ontology.xlsx")
write_xlsx(entry.wilcox.C2.sheets, "markers/entry_markers_C2_geneset.xlsx")
write_xlsx(entry.wilcox.H.sheets, "markers/entry_markers_Hallmarks.xlsx")
write_xlsx(entry.wilcox.sheets, "markers/entry_markers_wilcox.xlsx")

## [ Exit markers ]----

exitv2.seurat <- readRDS("../data/exit_v2_seurat_mutrans.rds")

Idents(exitv2.seurat) <- exitv2.seurat$attractor

exitv2.seurat <- NormalizeData(exitv2.seurat)
exitv2.seurat <- FindVariableFeatures(exitv2.seurat, nfeatures = 1000)
exitv2.seurat <- ScaleData(exitv2.seurat)

exitv2.wilcox.markers <- FindAllMarkers(exitv2.seurat, assay = "RNA", slot = "data",
                                        only.pos = TRUE, densify = TRUE, min.pct = 0.05,
                                        logfc.threshold = 0.2, test.use = "wilcox")

saveRDS(exitv2.wilcox.markers, "markers/exit_v2_markers_wilcox.rds")

exitv2.wilcox.list <- lapply(unique(exitv2.wilcox.markers$cluster), function(x){
  tmp.df <- exitv2.wilcox.markers[exitv2.wilcox.markers$cluster==x & exitv2.wilcox.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})

names(exitv2.wilcox.list) <- unique(exitv2.wilcox.markers$cluster)

saveRDS(exitv2.wilcox.list, "markers/exit_v2_markers_wilcox_list.rds")

table(exitv2.wilcox.markers$cluster)
## 0   1   2   3   4 
## 736 845 744 766 824

require(clusterProfiler)
require(org.Hs.eg.db)

## Gene Ontology biological processes
exitv2.wilcox.GO.BP.list <- lapply(exitv2.wilcox.list, function(x){
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

exitv2.wilcox.C2.list <- lapply(exitv2.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

## Hallmark gene sets
H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

exitv2.wilcox.H.list <- lapply(exitv2.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

saveRDS(exitv2.wilcox.GO.BP.list, "markers/exit_v2_markers_BP_ontology.rds")
saveRDS(exitv2.wilcox.C2.list, "markers/exit_v2_markers_C2_geneset.rds")
saveRDS(exitv2.wilcox.H.list, "markers/exit_v2_markers_Hallmarks.rds")

require(writexl)

exitv2.wilcox.lists <- list("exitv2.wilcox.GO.BP.list" = exitv2.wilcox.GO.BP.list,                       
                            "exitv2.wilcox.C2.list" = exitv2.wilcox.C2.list,
                            "exitv2.wilcox.H.list" = exitv2.wilcox.H.list,
                            "exitv2.wilcox.list" = exitv2.wilcox.list)

exitv2.wilcox.names <- sub("list", "sheets", names(exitv2.wilcox.lists))

for(i in 1:length(exitv2.wilcox.lists)){
  wilcox.list <- exitv2.wilcox.lists[[i]]
  list.sheets <- lapply(wilcox.list, as.data.frame)
  assign(exitv2.wilcox.names[[i]], list.sheets)
}

write_xlsx(exitv2.wilcox.GO.BP.sheets, "markers/exit_v2_markers_BP_ontology.xlsx")
write_xlsx(exitv2.wilcox.C2.sheets, "markers/exit_v2_markers_C2_geneset.xlsx")
write_xlsx(exitv2.wilcox.H.sheets, "markers/exit_v2_markers_Hallmarks.xlsx")
write_xlsx(exitv2.wilcox.sheets, "markers/exit_v2_markers_wilcox.xlsx")

