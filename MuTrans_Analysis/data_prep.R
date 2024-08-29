## MuTrans Data Prep

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)
library(scCustomize)

setwd("Z:/echen/Updated NB analysis 2")

## [ Separate data ] ----

## Read in final Seurat object
nb.seurat <- readRDS("data/nb_seurat.rds")

## Split objcect by "Condition"
seurat.list <- SplitObject(nb.seurat, split.by = "Condition")
list2env(seurat.list, envir = globalenv())

## List objects for DTP entry/exit datasets
ut.list <- list(POT, Untreated)
entry.list <- list(POT , `Cisplatin(1)_ON`, `Cisplatin(2)_ON`)
exit.list <- list(`Cisplatin(1)_ON`, `Cisplatin(2)_ON`, Cisplatin_1weeksOFF, Cisplatin_4weeksOFF)

ut.seurat <- Merge_Seurat_List(ut.list)
entry.seurat <- Merge_Seurat_List(entry.list)
exit.seurat <- Merge_Seurat_List(exit.list)

DefaultAssay(ut.seurat) <- "RNA"
ut.seurat <- JoinLayers(ut.seurat)
ut.seurat <- NormalizeData(ut.seurat)
ut.seurat[["RNA"]] <- as(ut.seurat[["RNA"]], "Assay")

saveRDS(ut.seurat, "data/nb_ut_seurat.rds")

DefaultAssay(entry.seurat) <- "RNA"
entry.seurat <- JoinLayers(entry.seurat)
entry.seurat <- NormalizeData(entry.seurat)
entry.seurat[["RNA"]] <- as(entry.seurat[["RNA"]], "Assay")

saveRDS(entry.seurat, "data/nb_entry_seurat.rds")

DefaultAssay(exit.seurat) <- "RNA"
exit.seurat <- JoinLayers(exit.seurat)
exit.seurat <- NormalizeData(exit.seurat)
exit.seurat[["RNA"]] <- as(exit.seurat[["RNA"]], "Assay")

saveRDS(exit.seurat, "data/nb_exit_seurat.rds")

## [ Seurat to AnnData ] ----

ut.seurat <- readRDS("data/nb_ut_seurat.rds")
SaveH5Seurat(ut.seurat, filename = "data/nb_ut_seurat.h5Seurat", overwrite = T, verbose = T)
Convert("data/nb_ut_seurat.h5Seurat", dest = "h5ad", assay = "RNA", overwrite = T)

entry.seurat <- readRDS("data/nb_entry_seurat.rds")
SaveH5Seurat(entry.seurat, filename = "data/nb_entry_seurat.h5Seurat", overwrite = T, verbose = T)
Convert("data/nb_entry_seurat.h5Seurat", dest = "h5ad", assay = "RNA", overwrite = T)

exit.seurat <- readRDS("data/nb_exit_seurat.rds")
SaveH5Seurat(exit.seurat, filename = "data/nb_exit_seurat.h5Seurat", overwrite = T, verbose = T)
Convert("data/nb_exit_seurat.h5Seurat", dest = "h5ad", assay = "RNA", overwrite = T)


