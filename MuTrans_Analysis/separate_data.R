## Separate data by condition

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)

nb.amt.seurat <- readRDS("nb_seurat_AMT.rds")

barcode.df <- read.csv("barcode_df.csv", row.names = 1)

## [ Separate data ] ----

nb.ut.seurat <- subset(x = nb.amt.seurat, subset = Condition == "Untreated" | Condition == "Untreated_rec")
nb.cis.seurat <- subset(x = nb.amt.seurat, subset = Condition == "Cisplatin" | Condition == "Cisplatin_rec")
nb.jq1.seurat <- subset(x = nb.amt.seurat, subset = Condition == "JQ1" | Condition == "JQ1_rec")

DefaultAssay(nb.ut.seurat) <- "RNA"
DefaultAssay(nb.cis.seurat) <- "RNA"
DefaultAssay(nb.jq1.seurat) <- "RNA"

nb.ut.seurat <- SCTransform(nb.ut.seurat, method = "glmGamPoi",
                            vst.flavor = "v2",
                            vars.to.regress = c("Seurat.Cycle.Score", "subsets_mt_percent"), 
                            verbose = TRUE,
                            do.scale = TRUE,
                            do.center = TRUE,
                            variable.features.n = 1000)

nb.cis.seurat <- SCTransform(nb.cis.seurat, method = "glmGamPoi",
                             vst.flavor = "v2",
                             vars.to.regress = c("Seurat.Cycle.Score", "subsets_mt_percent"), 
                             verbose = TRUE,
                             do.scale = TRUE,
                             do.center = TRUE,
                             variable.features.n = 1000)

nb.jq1.seurat <- SCTransform(nb.jq1.seurat, method = "glmGamPoi",
                             vst.flavor = "v2",
                             vars.to.regress = c("Seurat.Cycle.Score", "subsets_mt_percent"), 
                             verbose = TRUE,
                             do.scale = TRUE,
                             do.center = TRUE,
                             variable.features.n = 1000)

nb.ut.sce <- as.SingleCellExperiment(nb.ut.seurat)
nb.cis.sce <- as.SingleCellExperiment(nb.cis.seurat)
nb.jq1.sce <- as.SingleCellExperiment(nb.jq1.seurat)

## [ Add barcodes ] ----

ut.filt.sce <- nb.ut.sce[, colnames(nb.ut.sce) %in% rownames(barcode.df)]
ut.filt.sce$Full.BCS <- barcode.df[colnames(ut.filt.sce), ]$Full.BCS
saveRDS(ut.filt.sce, "ut_sce.rds")

cis.filt.sce <- nb.cis.sce[, colnames(nb.cis.sce) %in% rownames(barcode.df)]
cis.filt.sce$Full.BCS <- barcode.df[colnames(cis.filt.sce), ]$Full.BCS
saveRDS(cis.filt.sce, "cis_sce.rds")

jq1.filt.sce <- nb.jq1.sce[, colnames(nb.jq1.sce) %in% rownames(barcode.df)]
jq1.filt.sce$Full.BCS <- barcode.df[colnames(jq1.filt.sce), ]$Full.BCS
saveRDS(jq1.filt.sce, "jq1_sce.rds")
