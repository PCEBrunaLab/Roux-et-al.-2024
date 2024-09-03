## scRNA-seq analysis with Seurat

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## Any missing packages:

required.packages <- c("scater", "scran", "scuttle", "Seurat", "patchwork", 
                       "biomaRt", "reshape2", "EnsDb.Hsapiens.v86", "ggsankey", 
                       "tricycle", "hues", "glmGamPoi", "sctransform")

if (!all(required.packages %in% installed.packages()[,"Package"])){
  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]])
}

## [ Load dependencies ] ----

library(scater)
library(scran)
library(scuttle)
library(Seurat)
library(patchwork)
library(biomaRt)
library(reshape2)
library(hues)

nb.sce <- readRDS("nb_sce.RDS")

## Define lincRNA set: snapshot of the hg38 / GRCh38 version of the ensembl database (v86, 2017)
## https://bioconductor.org/packages/3.16/data/annotation/
require(EnsDb.Hsapiens.v86)
ens.86.genes <- genes(EnsDb.Hsapiens.v86)
linc.genes <- ens.86.genes$gene_id[ens.86.genes$gene_biotype=="lincRNA"]

sample.cols <- iwanthue(length(unique(nb.sce$sample_id)))
names(sample.cols) <- unique(nb.sce$sample_id)

## [ Data prep ] ----
## For cell lines and PDO models only
#ensembl <- useMart("ensembl", "hsapiens_gene_ensembl")
#ens.bm <- getBM(attributes = c("ensembl_gene_id", 
#                               "chromosome_name",
#                               "start_position",
#                               "end_position",
#                               "external_gene_name", 
#                               "hgnc_symbol", 
#                               "description"), 
#                mart = ensembl, 
#                filters = "ensembl_gene_id", 
#                values = rownames(nb.sce))
ens.bm <- read.csv("ensembl_biomart.csv", row.names = 1)

## Keep track of duplicates
dup.ensg <- ens.bm$ensembl_gene_id[duplicated(ens.bm$ensembl_gene_id)]
dup.ensg <- setNames(ens.bm$external_gene_name[ens.bm$ensembl_gene_id %in% dup.ensg],
                     ens.bm$ensembl_gene_id[ens.bm$ensembl_gene_id %in% dup.ensg])

ens.filt.bm <- ens.bm
ens.filt.bm$description <- gsub("(^.+) \\[Source.+", "\\1", ens.filt.bm$description)

ens.filt.bm <- ens.filt.bm[ens.filt.bm$chromosome_name %in% c(1:22, "X", "Y"), ]

## Remove anything not named
ens.filt.bm <- ens.filt.bm[!ens.filt.bm$external_gene_name=="", ]

## Keep track of everything that didn't get annotated
ens.removed.bm <- ens.bm[!ens.bm$ensembl_gene_id %in% ens.filt.bm$ensembl_gene_id, ]
nrow(ens.removed.bm)

length(grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE))
# 11008 / 11172 are novel transcripts or novel transcripts antisense to a gene

grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE, invert = TRUE)
## The rest is a mixed bag: ~160 genes either novel, mitochondrially encoded, non-coding, or just empty description & name

## Any duplicates?
sort(ens.filt.bm$external_gene_name[duplicated(ens.filt.bm$external_gene_name)])

#ARMCX5-GPRASP2
#CCL3L1 - part of chemokine cluster on chr17
#DNAJC9-AS1 - heatshock protein homolog -AS1 is readthrough
#ELFN2 - postsynaptic cell adhesion molecule in group III mGluRs
#GOLGA8M - something to do with the Golgi membrane
#GPR84-AS1 - a g-protein coupled receptor complex member readthrough
#MATR3 - nuclear matrix protein          
#MKKS - centrosomal shuttling protein
#NPIPA9 - nuclear pore complex related
#PRICKLE2-AS1 - homolog of prickle readthrough
#RAET1E-AS1 - MHC class I related genes member readthrough
#RGS5 - hypoxia-induced gene that is involved in the induction of endothelial apoptosis
#SPATA13 - spermatogenesis associated
#TBCE - tubulin folding cofactor E          
#TMSB15B - thymosin family member (there's a lot)
#TNFRSF10A-DT - TNF-receptor superfamily member; DT = divergent transcript (?)

#LINC00484
#LINC00486
#LINC00595
#LINC01115
#LINC01238
#LINC01605
#LINC03021
#LINC03023
#LINC03025

## Keep uniques
ens.filt.bm <- ens.filt.bm[!duplicated(ens.filt.bm$external_gene_name), ]
nrow(ens.filt.bm)

mt.genes <- ens.bm$ensembl_gene_id[grep("^MT-", ens.bm$external_gene_name)]
ribo.genes <- ens.bm$ensembl_gene_id[grep("^RP[SL]", ens.bm$external_gene_name)]
linc.genes <- linc.genes[linc.genes %in% ens.bm$ensembl_gene_id]



## For PDX models
#ENSEMBL Data Preparation
require(EnsDb.Hsapiens.v86)
ens.86.genes <- genes(EnsDb.Hsapiens.v86)
linc.genes <- ens.86.genes$gene_id[ens.86.genes$gene_biotype=="lincRNA"]

#We aligned these samples to both mouse and human genome
ensembl.h <- useMart("ensembl", "hsapiens_gene_ensembl")
ensembl.m <- useMart("ensembl", "mmusculus_gene_ensembl")

# Separate human and mouse gene IDs based on prefixes
human_gene_ids <- rownames(nb.sce)[grepl("^GRCh38", rownames(nb.sce))]
mouse_gene_ids <- rownames(nb.sce)[grepl("^GRCm39", rownames(nb.sce))]

# Remove prefixes
human_gene_ids_filter <- gsub(human_gene_ids, pattern = "(^GRCh38_)(S*)", replacement = "\\2")
mouse_gene_ids_filter <- gsub(mouse_gene_ids, pattern = "(^GRCm39_)(S*)", replacement = "\\2")

# Fetch human gene information
ens.h.bm <- getBM(attributes = c("ensembl_gene_id", 
                                   "chromosome_name",
                                   "start_position",
                                   "end_position",
                                   "external_gene_name", 
                                   "hgnc_symbol", 
                                   "description"), 
                    mart = ensembl.h, 
                    filters = "ensembl_gene_id", 
                    values = human_gene_ids_filter)

# Fetch mouse gene information
ens.m.bm <- getBM(attributes = c("ensembl_gene_id", 
                                   "chromosome_name",
                                   "start_position",
                                   "end_position",
                                   "external_gene_name", 
                                   "mgi_symbol", 
                                   "description"), 
                    mart = ensembl.m, 
                    filters = "ensembl_gene_id", 
                    values = mouse_gene_ids_filter)

#We have slight differences in column names, remove different column
ens.h.bm$hgnc_symbol <- NULL
ens.m.bm$mgi_symbol <- NULL

#Assign column to show which genome the gene has come from
ens.h.bm$genome <- "human"
ens.m.bm$genome <- "mouse"

# Combine the results from human and mouse
combined_ens.bm <- rbind(ens.h.bm, ens.m.bm)

#Save ensembl.csv for reuse to save on connection issues with ensembl
write.csv(combined_ens.bm, "datafiles/ensembl_pdx_mouse_and_human_biomart.csv")
gc()

#Identify duplicated gene names
dup.ensg <- combined_ens.bm$ensembl_gene_id[duplicated(combined_ens.bm$ensembl_gene_id)] 
dup.ensg <- setNames(combined_ens.bm$external_gene_name[combined_ens.bm$ensembl_gene_id %in% dup.ensg],
                     combined_ens.bm$ensembl_gene_id[combined_ens.bm$ensembl_gene_id %in% dup.ensg])
#NCMAP-DT, EIF1B-AS1, LINC00595 duplicated human
#AA536875 duplicated mouse

#Created filtered ens.bm and keep only useful information
ens.filt.bm <- combined_ens.bm
ens.filt.bm$description <- gsub("(^.+) \\[Source.+", "\\1", ens.filt.bm$description)
ens.filt.bm <- ens.filt.bm[ens.filt.bm$chromosome_name %in% c(1:22, "X", "Y"), ]

#Remove anything not named
ens.filt.bm <- ens.filt.bm[!ens.filt.bm$external_gene_name=="", ]

#Keep track of everything that didn't get annotated
ens.removed.bm <- combined_ens.bm[!combined_ens.bm$ensembl_gene_id %in% ens.filt.bm$ensembl_gene_id, ]
nrow(ens.removed.bm) #13037 transcripts not named

length(grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE))
#12763 are novel transcripts or novel transcripts antisense to a gene

grep("[Nn]ovel transcript", ens.removed.bm$description, value = TRUE, invert = TRUE)
#The rest is a mixed bag: ~275 genes either novel, mitochondrially encoded, non-coding, or just empty description & name

#Any duplicates?
sort(ens.filt.bm$external_gene_name[duplicated(ens.filt.bm$external_gene_name)])

#Keep uniques
ens.filt.bm <- ens.filt.bm[!duplicated(ens.filt.bm$external_gene_name), ]
nrow(ens.filt.bm) #59062

#Identify mitochondrial, ribosomal and linc genes
mt.genes <- ens.h.bm$ensembl_gene_id[grep("^MT-", ens.h.bm$external_gene_name)]
ribo.genes <- ens.h.bm$ensembl_gene_id[grep("^RP[SL]", ens.h.bm$external_gene_name)]
linc.genes <- linc.genes[linc.genes %in% ens.h.bm$ensembl_gene_id]

#Remove prefix from sce objects
rownames(nb.sce) <- gsub(rownames(nb.sce), pattern = "(GRC[hm3839]*)_(ENS.*)", replacement = "\\2")


## Basic QC
nb.sce <- addPerCellQCMetrics(nb.sce, flatten = TRUE, subsets = list(mt = mt.genes, linc = linc.genes, ribo = ribo.genes))

nb.sce <- nb.sce[ens.filt.bm$ensembl_gene_id,]
rownames(nb.sce) <- ens.filt.bm$external_gene_name
rownames(ens.filt.bm) <- ens.filt.bm$external_gene_name

rowData(nb.sce) <- ens.filt.bm

metadata.df <-  as.data.frame(colData(nb.sce))[, c("Sample", "Barcode", "Batch", "Class", "sample_id",
                                                   "description", "Rec", "Condition", "detected", "total", 
                                                   "subsets_mt_percent", "subsets_ribo_percent", "subsets_linc_percent")]

nb.seurat <- CreateSeuratObject(counts = assay(nb.sce, "counts"),
                                assay = "RNA",
                                meta.data = metadata.df)

## [ QC ] ----

det.sce.gg <- plotColData(nb.sce, y = "detected", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "Genes / cell", title = "Detected genes per cell") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

umi.sce.gg <- plotColData(nb.sce, y = "total", x = "Condition", colour_by = "Condition") + 
  labs(x = element_blank(), y = "UMI / cell", title = "Detected UMI per cell") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5) +
  scale_y_continuous(trans = "log10", labels = scales::comma) +
  annotation_logticks(sides = "l", outside = TRUE) + coord_cartesian(clip = "off")

mt.sce.gg <- plotColData(nb.sce, y = "subsets_mt_percent", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "Mitochondrial\nexpression % of total", 
       title = "Proportion of mitochondrial genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

rb.sce.gg <- plotColData(nb.sce, y = "subsets_ribo_percent", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "Ribosomal\nexpression % of total", 
       title = "Proportion of ribosomal genes\nin total cell expression") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

det.sce.gg + umi.sce.gg + mt.sce.gg + rb.sce.gg + plot_layout(nrow = 2) & theme(aspect.ratio = 0.5)

## Examine whether there are any genes that dominate expression in cells - large % of expression occupied by single gene
counts.assay <- counts(nb.sce)
counts.assay@x <- counts.assay@x/rep.int(colSums(counts.assay), diff(counts.assay@p))
top.expr <- order(Matrix::rowSums(counts.assay), decreasing = TRUE)[20:1]
top.expr <- as.matrix(t(counts.assay[top.expr, ]))
top.expr <- as.data.frame(top.expr)
ord.idx <- colnames(top.expr)
top.expr <- melt(top.expr, variable.name = "Gene", value.name = "Prop")
top.expr$Gene <- factor(top.expr$Gene, levels = ord.idx)

ggplot(top.expr,
       mapping = aes(x = Gene, y = Prop * 100)) +
  geom_boxplot() +
  theme_bw() + theme(panel.grid.minor.x = element_blank(),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.major.y = element_line(linetype = "dotted")) +
  coord_flip() +
  labs(y = "% total expression", x = "Gene", title = "Percentage expression of a single gene\nper total cell expression")

rm(counts.assay)
gc()

malat1.sce.gg <- plotExpression(nb.sce, 
               features = "MALAT1", 
               exprs_values = "logcounts",
               x = "Condition", 
               colour_by = "Condition") +
  labs(x = element_blank(), y = "Log2 Expression") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

lincrna.sce.gg <- plotColData(nb.sce, y = "subsets_linc_percent", x = "Condition", colour_by = "Condition") +
  labs(x = element_blank(), y = "lincRNA\nexpression % of total") +
  scale_colour_manual(values = group.cols) +
  theme(legend.position = "none", axis.text.x = element_text(size = 6)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.9, linewidth = 0.5)

malat1.sce.gg + lincrna.sce.gg + plot_layout(nrow = 2) & theme(aspect.ratio = 0.5)

## Cell Cycle prediction
nb.cycle.seurat <- NormalizeData(nb.seurat)
nb.cycle.seurat <- CellCycleScoring(nb.cycle.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

nb.sce$Seurat.Phase <- nb.cycle.seurat$Phase
nb.sce$Seurat.S <- nb.cycle.seurat$S.Score
nb.sce$Seurat.G2M <- nb.cycle.seurat$G2M.Score

nb.seurat$Seurat.Phase <- nb.cycle.seurat$Phase
nb.seurat$Seurat.S <- nb.cycle.seurat$S.Score
nb.seurat$Seurat.G2M <- nb.cycle.seurat$G2M.Score

## [ Filtering ] ----

## Doublet detection already done
table(nb.sce$Class)

sample.cols <- iwanthue(length(names(table(nb.sce$sample_id))))
names(sample.cols) <- names(table(nb.sce$sample_id))

## Setting the expression threshold over 1000 doesn't remove any cells
detected.filt <- colnames(nb.sce)[nb.sce$detected > 1000]
## Genes have to have at least 5 reads 
rowcounts.filt <- rownames(nb.sce)[Matrix::rowSums(counts(nb.sce)) > 5]
## No MT genes in the data, they've been filtered out already
mt.genes.filt <- grep("^MT-", rownames(nb.sce), invert = TRUE, value = TRUE)
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

nb.filt.sce <- nb.sce[genes.filt, detected.filt]

## ~ 6000 genes lost
dim(nb.sce) - dim(nb.filt.sce)

mito.filt <- colnames(nb.filt.sce)[nb.filt.sce$subsets_mt_percent < 10]
ribo.filt <- colnames(nb.filt.sce)[nb.filt.sce$subsets_ribo_percent > 5]
qc.filt <- intersect(mito.filt, ribo.filt)

nb.filt.sce <- nb.filt.sce[,qc.filt]

## ~ 2600 cells lost
dim(nb.sce) - dim(nb.filt.sce)

## The same filtering in Seurat
rowcounts.filt <- rownames(nb.seurat)[Matrix::rowSums(GetAssayData(nb.seurat, slot = "counts", assay = "RNA")) > 5]
mt.genes.filt <- grep("^MT-", rownames(nb.seurat), invert = TRUE, value = TRUE)
genes.filt <- intersect(mt.genes.filt, rowcounts.filt)

nb.seurat$nFeature_RNA == nb.seurat$detected

nb.filt.seurat <- subset(nb.seurat, 
                         cells = WhichCells(nb.seurat, 
                                            expression = subsets_mt_percent < 10 &
                                              subsets_ribo_percent > 5 &
                                              detected > 1000),
                         features = genes.filt)

## The dimensions are equal
dim(nb.filt.seurat) == dim(nb.filt.sce)

## [ Normalisation ] ----

DefaultAssay(nb.filt.seurat) <- "RNA"
nb.filt.seurat$Seurat.Cycle.Score <- nb.filt.seurat$Seurat.S - nb.filt.seurat$Seurat.G2M

## Normalise the data regressing out the cell cycle and mitochondrial proportion
## and calculate the top 1000 hypervariable genes 
nb.filt.seurat <- SCTransform(nb.filt.seurat, method = "glmGamPoi",
                              vst.flavor = "v2",
                              vars.to.regress = c("Seurat.Cycle.Score", "subsets_mt_percent"), 
                              verbose = TRUE,
                              do.scale = TRUE,
                              do.center = TRUE,
                              variable.features.n = 1000)

## Do PCA and look at the first 50 dims
nb.filt.seurat <- RunPCA(nb.filt.seurat, verbose = FALSE, npcs = 100, 
                         features = nb.filt.seurat@assays$SCT@var.features)
ElbowPlot(object = nb.filt.seurat, ndims = 50, reduction = "pca")

## Run Harmony to do batch correction
require(harmony)
nb.harmony.seurat <- RunHarmony(nb.filt.seurat, 
                                group.by.vars = "Batch",
                                assay.use = "SCT", 
                                reduction.save = "Harmony",
                                theta = 2, 
                                max.iter.harmony = 50, 
                                plot_convergence = FALSE)

nb.harmony.seurat <- RunUMAP(nb.harmony.seurat, reduction = "Harmony", dims = 1:30)
nb.harmony.seurat <- FindNeighbors(nb.harmony.seurat, reduction = "Harmony", dims = 1:30)

## Seurat can generate more or fewer clusters based on a resolution parameter, this does a param search across a range
cluster.search <- function(seurat, from = 0.2, to = 2, by = 0.2){
  tmp.res <- lapply(seq(from = from, to = to, by = by), function(x){
    tmp.clust <- FetchData(FindClusters(seurat, verbose = TRUE, resolution = x), 
                           c(paste0("SCT_snn_res.", x), "seurat_clusters"))
    colnames(tmp.clust) <- c(paste0("snn_res.", x), paste0("seurat_clusters.", x))
    return(tmp.clust)
  })
  tmp.res <- do.call("cbind", tmp.res)
  return(AddMetaData(seurat, metadata = tmp.res))
}

nb.harmony.seurat <- cluster.search(nb.harmony.seurat)

## Custom ggplot theme
umap.theme <- function() {
  require(ggplot2)
  return(theme_bw() +
           theme(aspect.ratio = 1, 
                 panel.grid = element_blank(),
                 panel.border = element_blank(),
                 axis.line = element_line(colour = "#16161D", linewidth = 0.8),
                 axis.ticks = element_line(colour = "#16161D", linewidth = 0.8),
                 legend.title = element_blank()))
}

## Extended Data Figure 3e
dim1.gg <- DimPlot(nb.harmony.seurat, group.by = "Condition", order = TRUE) +
  scale_colour_manual(values = group.cols) +
  umap.theme() + labs(title = "Conditions") +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 5), title = "Condition"))

## Figure 3a
dim2.gg <- DimPlot(nb.harmony.seurat,
                   group.by = "seurat_clusters.0.2", order = TRUE) +
  scale_colour_manual(values = iwanthue(length(unique(nb.harmony.seurat$seurat_clusters.0.2)))) +
  umap.theme() + labs(title = "Clusters resolution = 0.2") +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 5), title = "Cluster Number"))

dim3.gg <- DimPlot(nb.harmony.seurat, group.by = "Condition", split.by = "sample_id", order = TRUE) +
  scale_colour_manual(values = group.cols) +
  umap.theme() + labs(title = "Samples") + theme(legend.position = "bottom") +
  facet_wrap(~sample_id, nrow = 3) +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 5), title = "Condition"))

require(patchwork)
dim1.gg | dim2.gg

dim3.gg

saveRDS(nb.harmony.seurat, "nb_seurat_harmony_dimred.rds")

cool.warm.pal <- colorRampPalette(c("#1d4877", "#1b8a5a", "#fbb021", "#f68838", "#ee3e32"))

feat1.gg <- FeaturePlot(nb.harmony.seurat,
            features = "detected",
            order = TRUE,
            min.cutoff = "q5", max.cutoff = "q90") +
  scale_colour_gradientn(colours = cool.warm.pal(100)) +
  umap.theme() +
  labs(title = "Detected Features")

feat2.gg <- FeaturePlot(nb.harmony.seurat,
                        features = "subsets_mt_percent",
                        order = TRUE,
                        min.cutoff = "q5", max.cutoff = "q90") +
  scale_colour_gradientn(colours = cool.warm.pal(100)) +
  umap.theme() +
  labs(title = "% Mitochondrial genes / cell")

feat3.gg <- FeaturePlot(nb.harmony.seurat,
                        features = "subsets_ribo_percent",
                        order = TRUE,
                        min.cutoff = "q5", max.cutoff = "q90") +
  scale_colour_gradientn(colours = cool.warm.pal(100)) +
  umap.theme() +
  labs(title = "% Ribosomal genes / cell")

feat4.gg <- FeaturePlot(nb.harmony.seurat,
                        features = "Seurat.S",
                        order = TRUE,
                        min.cutoff = "q5", max.cutoff = "q90") +
  scale_colour_gradientn(colours = cool.warm.pal(100)) +
  umap.theme() +
  labs(title = "Seurat Cycle S-phase")

feat5.gg <- FeaturePlot(nb.harmony.seurat,
                        features = "Seurat.G2M",
                        order = TRUE,
                        min.cutoff = "q5", max.cutoff = "q90") +
  scale_colour_gradientn(colours = cool.warm.pal(100)) +
  umap.theme() +
  labs(title = "Seurat Cycle G2M-phase")

feat6.gg <- FeaturePlot(nb.harmony.seurat,
                        features = "subsets_linc_percent",
                        order = TRUE,
                        min.cutoff = "q5", max.cutoff = "q90") +
  scale_colour_gradientn(colours = pals::ocean.haline(100)) +
  umap.theme() +
  labs(title = "% lincRNA genes / cell")

dim4.gg <- DimPlot(nb.harmony.seurat,
                    group.by = "Seurat.Phase", order = TRUE) +
  scale_colour_manual(values = iwanthue(length(unique(nb.harmony.seurat$Seurat.Phase)))) +
  umap.theme() + labs(title = "Seurat Cycle Phases")

feat1.gg + feat2.gg + feat3.gg +
  feat4.gg + feat5.gg + dim4.gg +
  plot_layout(nrow = 2, ncol = 3)

dim2.gg <- DimPlot(nb.harmony.seurat,
                   group.by = "seurat_clusters.0.2", order = TRUE, label = TRUE, repel = TRUE) +
  scale_colour_manual(values = iwanthue(length(unique(nb.harmony.seurat$seurat_clusters.0.2)))) +
  umap.theme() + labs(title = "10x Harmony Clusters res=0.2")

dim2.gg$layers[[2]]$aes_params$bg.colour <- "white"
dim2.gg$layers[[2]]$aes_params$direction <- "x"
dim2.gg$layers[[2]]$geom_params$max.overlaps <- 100

dim1.gg + dim2.gg + feat3.gg + feat6.gg + plot_layout(nrow = 2, ncol = 2)

## [ Markers ] ----

nb.harmony.seurat$Test_Clusters <- nb.harmony.seurat$seurat_clusters.0.2
Idents(nb.harmony.seurat) <- nb.harmony.seurat$Test_Clusters

nb.clust.wilcox.markers <- FindAllMarkers(nb.harmony.seurat,
                                          assay = "SCT",
                                          slot = "data",
                                          only.pos = TRUE,
                                          densify = TRUE,
                                          min.pct = 0.05,
                                          logfc.threshold = 0.2,
                                          test.use = "wilcox")

nb.clust.wilcox.list <- lapply(unique(nb.clust.wilcox.markers$cluster), function(x){
  tmp.df <- nb.clust.wilcox.markers[nb.clust.wilcox.markers$cluster==x & nb.clust.wilcox.markers$p_val_adj < 0.05, ]
  tmp.df[order(tmp.df$avg_log2FC, decreasing = TRUE),]
})

names(nb.clust.wilcox.list) <- unique(nb.clust.wilcox.markers$cluster)

require(clusterProfiler)
require(org.Hs.eg.db)
nb.clust.wilcox.GO.BP.list <- lapply(nb.clust.wilcox.list, function(x){
  enrichGO(x$gene,
           OrgDb = "org.Hs.eg.db",
           keyType = "SYMBOL",
           ont = "BP", 
           readable = TRUE, 
           pAdjustMethod = "BH")
})

require(msigdbr)
C2.msigsdb <- msigdbr(species = "Homo sapiens", category = "C2")
C2.df <- as.data.frame(C2.msigsdb[,c("gs_name", "gene_symbol")])

nb.clust.wilcox.C2.list <- lapply(nb.clust.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = C2.df)
})

H.msigsdb <- msigdbr(species = "Homo sapiens", category = "H")
H.df <- as.data.frame(H.msigsdb[,c("gs_name", "gene_symbol")])

nb.clust.wilcox.H.list <- lapply(nb.clust.wilcox.list, function(x){
  enricher(x$gene,
           pAdjustMethod = "BH", TERM2GENE = H.df)
})

saveRDS(nb.clust.wilcox.GO.BP.list, "nb_seurat_markers_BP_ontology.rds")
saveRDS(nb.clust.wilcox.C2.list, "nb_seurat_markers_C2_geneset.rds")
saveRDS(nb.clust.wilcox.H.list, "nb_seurat_markers_Hallmarks.rds")
saveRDS(nb.clust.wilcox.list, "nb_seurat_markers_wilcox.rds")

require(writexl)
nb.clust.wilcox.lists <- list("nb.clust.wilcox.GO.BP.list" = nb.clust.wilcox.GO.BP.list,
                              "nb.clust.wilcox.C2.list" = nb.clust.wilcox.C2.list,
                              "nb.clust.wilcox.H.list" = nb.clust.wilcox.H.list,
                              "nb.clust.wilcox.list" = nb.clust.wilcox.list)

nb.clust.wilcox.names <- sub("list", "sheets", names(nb.clust.wilcox.lists))

for(i in 1:length(nb.clust.wilcox.lists)){
  wilcox.list <- cnb.clust.wilcox.lists[[i]]
  list.sheets <- lapply(wilcox.list, as.data.frame)
  assign(nb.clust.wilcox.names[[i]], list.sheets)
}

write_xlsx(nb.clust.wilcox.GO.BP.sheets, "nb_seurat_markers_BP_ontology.xlsx")
write_xlsx(nb.clust.wilcox.C2.sheets, "nb_seurat_markers_C2_geneset.xlsx")
write_xlsx(nb.clust.wilcox.H.sheets, "nb_seurat_markers_Hallmarks.xlsx")
write_xlsx(nb.clust.wilcox.sheets, "nb_seurat_markers_wilcox.xlsx")

## [ ADRN/MES signatures ] ----

require(readxl)
vanGroningen.df <- as.data.frame(read_xlsx(path = "vanGroningen_2017.xlsx", col_names = FALSE))
colnames(vanGroningen.df) <- c("Gene", "Signature")

vanGroningen.sig <- lapply(unique(vanGroningen.df$Signature), function(x){
  vanGroningen.df[vanGroningen.df$Signature==x, "Gene"]
})

names(vanGroningen.sig) <- unique(vanGroningen.df$Signature)
vanGroningen.sig

nb.harmony.seurat <- AddModuleScore(nb.harmony.seurat, features = vanGroningen.sig, assay = "SCT", seed = 12345, 
                                    name = c("MES.Sig", "ADRN.Sig"))
colnames(nb.harmony.seurat@meta.data) <- gsub("\\.Sig[1-9]$", "\\.Sig", colnames(nb.harmony.seurat@meta.data))

require(ggbeeswarm)
vln1.gg <- VlnPlot(nb.harmony.seurat, 
        features = c("MES.Sig"), group.by = "Test_Clusters") &
  umap.theme() + theme(aspect.ratio = 0.5)

vln1.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln1.gg$layers[[2]] <- NULL

vln2.gg <- VlnPlot(nb.harmony.seurat, 
                   features = c("ADRN.Sig"), group.by = "Test_Clusters") &
  umap.theme() + theme(aspect.ratio = 0.5)

vln2.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln2.gg$layers[[2]] <- NULL

vln1.gg | vln2.gg

dim5.gg <- DimPlot(nb.harmony.seurat,
                   group.by = "Test_Clusters", order = TRUE, label = TRUE, repel = TRUE) +
  umap.theme() + labs(title = "Working Clusters")

dim5.gg

patch.design <- 
  "
  12
  13
"

dim5.gg + vln1.gg + vln2.gg + plot_layout(design = patch.design)

## Extended Data Figure 3f
mes.gg <- FeaturePlot(nb.harmony.seurat,
                      features = "MES.Sig",
                      order = TRUE,
                      min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.haline(100)) +
  umap.theme()

## Extended Data Figure 3f
adrn.gg <- FeaturePlot(nb.harmony.seurat,
                       features = "ADRN.Sig",
                       order = TRUE,
                       min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.matter(100)) +
  umap.theme()

dim5.gg + mes.gg + adrn.gg + plot_layout(design = patch.design)

## [ AMT signature ] ----

## AMT score = MES score - ADRN score
nb.harmony.seurat$AMT.Sig <- nb.harmony.seurat$MES.Sig - nb.harmony.seurat$ADRN.Sig

sig.gg <- FeaturePlot(nb.harmony.seurat,
                      features = "AMT.Sig",
                      order = TRUE,
                      min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(title = "AMT score = MES expression - ADRN expression", colour = "Expression") +
  umap.theme()

sig.gg + theme(text = element_text(size = 15),
               legend.key.size = unit(1, "cm"))

## [ AMT score ] ----

## Calculate an AMT score by testing for a difference between the mean expression of MES and ADRN gene sets
nb.counts <- nb.harmony.seurat@assays$SCT$data

t.stats <- matrix(data = NA, nrow = ncol(nb.counts), ncol = 2)
rownames(t.stats) <- colnames(nb.counts)
colnames(t.stats) <- c("t", "p.value")

nb.counts.mes <- nb.counts[rownames(nb.counts) %in% vanGroningen.sig$MES, ]
nb.counts.adrn <- nb.counts[rownames(nb.counts) %in% vanGroningen.sig$ADRN, ]

for(cell in 1:ncol(nb.counts)) {
  result <- t.test(nb.counts.mes[, cell], nb.counts.adrn[, cell])
  t.stats[cell, 1] <- result$statistic
  t.stats[cell, 2] <- result$p.value
}
saveRDS(t.stats, "AMT_t_statistics.rds")

t.stats.df <- as.data.frame(t.stats)

require(dplyr)
amt.states <- t.stats.df %>%
  mutate(ADRN.MES = case_when(
    t > 0 & p.value < 0.05 ~ "MES",
    t > 0 & p.value >= 0.05 | t < 0 & p.value >= 0.05 ~ "intermediate",
    t < 0 & p.value < 0.05 ~ "ADRN"
  )) %>%
  arrange(ADRN.MES)

## Extended Data Figure 3g
amt.states.gg <-
  amt.states %>%
  ggplot(aes(ADRN.MES, fill = ADRN.MES)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values = c("purple","grey","orange")) +
  scale_x_discrete(labels = c("ADRN", "Intermediate", "MES")) +
  labs(x = "",
       y = "Frequency",
       title = "3 states",
       fill = "Annotation") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

amt.states.gg + theme(text = element_text(size = 15))

amt.states <- amt.states[match(colnames(nb.harmony.seurat), rownames(amt.states)), ]

nb.harmony.seurat$AMT.score <- amt.states$t
nb.harmony.seurat$AMT.state <- amt.states$ADRN.MES

## Extended Data Figure 3f
amt.gg <- FeaturePlot(nb.harmony.seurat,
                      features = "AMT.score",
                      min.cutoff = -15, max.cutoff = 15) +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(title = "AMT score", colour = "Score") +
  umap.theme()

amt.gg + theme(text = element_text(size = 15),
                legend.key.size = unit(1, "cm"))

amt.cols <- c("ADRN" = "purple",
              "intermediate" = "grey",
              "MES" = "orange")

## Figure 3c
adrn.mes.gg <- DimPlot(nb.harmony.seurat,
                       group.by = "AMT.state") +
  scale_colour_manual(values = amt.cols, labels = c("Adrenergic", "Intermediate", "Mesenchymal")) +
  umap.theme() +
  labs(title = "AMT state")

adrn.mes.gg +
  theme(text = element_text(size = 15)) +
  guides(colour = guide_legend(override.aes = list(size = 5), title = "State"))

saveRDS(nb.harmony.seurat, "nb_seurat_AMT.rds")

nb.amt.sce <- as.SingleCellExperiment(nb.harmony.seurat)
saveRDS(nb.amt.sce, "nb_sce_AMT.rds")


## [ Dyer Signature ] ----
Dyer signature colours as taken from paper
dyer.cols <- c("#27296f", "#7e171a", "#009169")


