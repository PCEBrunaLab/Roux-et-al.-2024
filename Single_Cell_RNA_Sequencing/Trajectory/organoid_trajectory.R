## NB organoids trajectory analysis

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)
library(SeuratDisk)
library(scater)
library(scCustomize)
library(TSCAN)
library(slingshot)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(ggVennDiagram)
library(scattermore)
library(patchwork)
library(pheatmap)
library(hues)
library(viridis)
library(readxl)

setwd("~/my_rds/NB PDO analysis")

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

identity.cols <- c("SYMP-Dyer cycling1" = "#0077BB",
                   "ADRN stress" = "#882255",
                   "ADRN Myc" = "#BB22EE",
                   "SYMP-Dyer cycling2" = "#225555",
                   "MES-Dyer recovery1" = "#CCAA11",
                   "SYMP-Dyer cycling3" = "#332288",
                   "MES-Dyer recovery2" = "#EE7733",
                   "ADRN p53" = "#EE3377",
                   "MES-Dyer recovery3" = "#CC3311")

condition.cols <- c("NB039_cisplatin" = "#DDAA33",
                    "NB039_cisplatin_recovery" = "#9A6A20",
                    "NB039_untreated" = "#43ACFF",
                    "NB067_cisplatin" = "#BBCC33",
                    "NB067_cisplatin_recovery" = "#758020",
                    "NB067_untreated" = "#004488")

## [ Organoid signature plots ] ----

organoid.seurat <- readRDS("data/nb_seurat_AMT_non-batch_organoids.rds")

## Cluster colours
#cluster.cols <- iwanthue(length(unique(organoid.seurat$seurat_clusters)))
#saveRDS(cluster.cols, "data/organoid_cluster_cols.rds")
cluster.cols <- readRDS("data/organoid_cluster_cols.rds")

num.gg <- DimPlot(organoid.seurat,
                  group.by = "seurat_clusters", order = TRUE) +
  umap.theme() + labs(title = "Cluster number") +
  scale_colour_manual(values = cluster.cols)

num.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/organoid_umap_cluster_number.pdf", width = 5.8, height = 5.8)
ggsave("plots/organoid_umap_cluster_number.png", width = 5.8, height = 5.8)

adrn.gg <- FeaturePlot(organoid.seurat,
                       features = "ADRN.Sig",
                       order = TRUE,
                       min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.matter(100)) +
  labs(title = "vanGroningen\nADRN signature", colour = "Expression") +
  umap.theme()

mes.gg <- FeaturePlot(organoid.seurat,
                      features = "MES.Sig",
                      order = TRUE,
                      min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.haline(100)) +
  labs(title = "vanGroningen\nMES signature", colour = "Expression") +
  umap.theme()

adrn.gg + mes.gg &
  theme(text = element_text(size = 15),
        legend.key.size = unit(1, "cm"))

ggsave("plots/signatures/organoid_vG_ADRN_MES_umaps.pdf", width = 5.8, height = 8.3)
ggsave("plots/signatures/organoid_vG_ADRN_MES_umaps.png", width = 5.8, height = 8.3)

vln1.gg <- VlnPlot(organoid.seurat,
                   features = "MES.Sig", group.by = "seurat_clusters") +
  scale_fill_manual(values = cluster.cols) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln1.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln1.gg$layers[[2]] <- NULL

vln2.gg <- VlnPlot(organoid.seurat,
                   features = "ADRN.Sig", group.by = "seurat_clusters") +
  scale_fill_manual(values = cluster.cols) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln2.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln2.gg$layers[[2]] <- NULL

vln3.gg <- VlnPlot(organoid.seurat,
                   features = "AMT.score", group.by = "seurat_clusters") +
  scale_fill_manual(values = cluster.cols) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln3.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln3.gg$layers[[2]] <- NULL

(vln3.gg | (vln1.gg / vln2.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/organoid_vG_ADRN_MES_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/organoid_vG_ADRN_MES_violins.png", width = 8.3, height = 8.3)

amt.gg <- DimPlot(organoid.seurat,
                  group.by = "AMT.state", order = TRUE) +
  umap.theme() + labs(title = "AMT state") +
  scale_colour_manual(values = c("#990099", "#F37735", "lightgrey"),
                      labels = c("ADRN", "MES", "Intermediate"))

amt.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/organoid_umap_amt_state.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/organoid_umap_amt_state.png", width = 5.8, height = 5.8)

organoid.seurat$Cluster_Name <- recode(as.character(organoid.seurat$seurat_clusters),
                                  "0" = "SYMP-Dyer cycling1",
                                  "1" = "ADRN stress",
                                  "2" = "ADRN Myc",
                                  "3" = "SYMP-Dyer cycling2",
                                  "4" = "MES-Dyer recovery1",
                                  "5" = "SYMP-Dyer cycling3",
                                  "6" = "MES-Dyer recovery2",
                                  "7" = "ADRN p53",
                                  "8" = "MES-Dyer recovery3")

cluster.gg <- DimPlot(organoid.seurat,
                      group.by = "Cluster_Name", order = TRUE) +
  umap.theme() + labs(title = "Initial Cluster Annotation") +
  scale_colour_manual(values = identity.cols)

cluster.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/organoid_umap_initial_clusters.pdf", width = 5.8, height = 5.8)
ggsave("plots/organoid_umap_initial_clusters.png", width = 5.8, height = 5.8)

## [ Organoid coexpression ] ----

## Which cells coexpress van Groningen ADRN and MES signatures
coexprs = function(x){
  if (x[["MES.Sig"]] > 0.1 & x[["ADRN.Sig"]] > 0.1) {
    result <- 1
  } else {
    result <- 0
  }
  return(result)
}

metadata.df <- subset(organoid.seurat@meta.data, select = c("MES.Sig", "ADRN.Sig"))

coexprs.df <- apply(metadata.df, 1, coexprs)
organoid.seurat@meta.data <- cbind(organoid.seurat@meta.data, Coexprs.Sig = coexprs.df)

coexprs.gg <- DimPlot(organoid.seurat,
                      group.by = "Coexprs.Sig", order = TRUE) +
  umap.theme() + labs(title = "ADRN and MES Coexpression\n(MES > 0.1, ADRN > 0.1)") +
  scale_colour_manual(values = c("#4477AA", "#EE6677"), labels = c("FALSE", "TRUE")) + 
  labs(colour = "Coexpression")

coexprs.gg

ggsave("plots/coexpression/organoid_umap_coexpress.pdf", width = 5.8, height = 5.8)
ggsave("plots/coexpression/organoid_umap_coexpress.png", width = 5.8, height = 5.8)

saveRDS(organoid.seurat, "data/organoid_seurat_start.rds")

## [ Separated data ] ----

## Split object by organoid model
organoid.seurat <- readRDS("data/organoid_seurat_start.rds")

seurat.list <- SplitObject(organoid.seurat, split.by = "Model")
list2env(seurat.list, envir = globalenv())

saveRDS(NB039, "data/nb039_seurat.rds")
saveRDS(NB067, "data/nb067_seurat.rds")

## NB039
nb039.seurat <- readRDS("data/nb039_seurat.rds")

nb039.seurat <- RunPCA(nb039.seurat, verbose = FALSE, npcs = 100, 
                       features = nb039.seurat@assays$SCT@var.features)
ElbowPlot(object = nb039.seurat, ndims = 50, reduction = "pca")

nb039.seurat <- RunUMAP(nb039.seurat, reduction = "pca", dims = 1:30)
nb039.seurat <- FindNeighbors(nb039.seurat, reduction = "pca", dims = 1:30)
nb039.seurat <- FindClusters(nb039.seurat, resolution = 0.3)

DimPlot(nb039.seurat,
        group.by = "seurat_clusters", order = TRUE) +
  umap.theme() + labs(title = "Cluster number") +
  scale_colour_manual(values = iwanthue(length(unique(nb039.seurat$seurat_clusters)))) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/nb039_umap_cluster_number.pdf", width = 5.8, height = 5.8)
ggsave("plots/nb039_umap_cluster_number.png", width = 5.8, height = 5.8)

DimPlot(nb039.seurat,
        group.by = "Cluster_Name", order = TRUE) +
  umap.theme() + labs(title = "Initial Cluster Annotation") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/nb039_umap_initial_clusters.pdf", width = 5.8, height = 5.8)
ggsave("plots/nb039_umap_initial_clusters.png", width = 5.8, height = 5.8)

saveRDS(nb039.seurat, "data/nb039_seurat_start.rds")

nb039.sce <- as.SingleCellExperiment(nb039.seurat)
reducedDim(nb039.sce, "PCA") <- Embeddings(nb039.seurat, "pca")
reducedDim(nb039.sce, "UMAP") <- Embeddings(nb039.seurat, "umap")
saveRDS(nb039.sce, "data/nb039_sce_start.rds")

## NB067
nb067.seurat <- readRDS("data/nb067_seurat.rds")

nb067.seurat <- RunPCA(nb067.seurat, verbose = FALSE, npcs = 100, 
                       features = nb067.seurat@assays$SCT@var.features)
ElbowPlot(object = nb067.seurat, ndims = 50, reduction = "pca")

nb067.seurat <- RunUMAP(nb067.seurat, reduction = "pca", dims = 1:30)
nb067.seurat <- FindNeighbors(nb067.seurat, reduction = "pca", dims = 1:30)
nb067.seurat <- FindClusters(nb067.seurat, resolution = 0.3)

DimPlot(nb067.seurat,
        group.by = "seurat_clusters", order = TRUE) +
  umap.theme() + labs(title = "Cluster number") +
  scale_colour_manual(values = iwanthue(length(unique(nb067.seurat$seurat_clusters)))) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/nb067_umap_cluster_number.pdf", width = 5.8, height = 5.8)
ggsave("plots/nb067_umap_cluster_number.png", width = 5.8, height = 5.8)

DimPlot(nb067.seurat,
        group.by = "Cluster_Name", order = TRUE) +
  umap.theme() + labs(title = "Initial Cluster Annotation") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/nb067_umap_initial_clusters.pdf", width = 5.8, height = 5.8)
ggsave("plots/nb067_umap_initial_clusters.png", width = 5.8, height = 5.8)

saveRDS(nb067.seurat, "data/nb067_seurat_start.rds")

nb067.sce <- as.SingleCellExperiment(nb067.seurat)
reducedDim(nb067.sce, "PCA") <- Embeddings(nb067.seurat, "pca")
reducedDim(nb067.sce, "UMAP") <- Embeddings(nb067.seurat, "umap")
saveRDS(nb067.sce, "data/nb067_sce_start.rds")

## [ Combined trajectory inference ] ----

organoid.sce <- readRDS("data/organoid_sce_start.rds")

organoid.slingshot.sce <- slingshot(organoid.sce, reducedDim = "PCA",
                                    clusterLabels = "seurat_clusters", start.clus = 0,
                                    stretch = 0)

saveRDS(organoid.slingshot.sce, "data/organoid_sce_slingshot_no_stretch.rds")

organoid.slingshot.sce <- readRDS("data/organoid_sce_slingshot_no_stretch.rds")

## Extract matrix of pseudotime values along each lineage
slingshot.lineages <- slingPseudotime(organoid.slingshot.sce)
head(slingshot.lineages)
## 4 lineages

## Single pseudotime for all cells
shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)

umap.gg <- plotUMAP(organoid.slingshot.sce, colour_by = I(shared.pseudotime)) +
  scale_colour_viridis() +
  labs(colour = "Pseudotime")
slingshot.embedded.all <- embedCurves(organoid.slingshot.sce, "UMAP")
slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
for (lineage in slingshot.embedded.all) {
  slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
  umap.gg <- umap.gg + geom_path(data = slingshot.embedded.all,
                                 aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2)
}

saveRDS(slingshot.embedded.all, "data/organoid_slingshot_embedded_all_no_stretch.rds")

umap.gg + theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/organoid_umap_slingshot_no_stretch.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/organoid_umap_slingshot_no_stretch.png", width = 8.3, height = 5.8)

## Save as SDS
slingshot.embedded.all.sds <- embedCurves(organoid.slingshot.sce, "UMAP")
slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)
saveRDS(slingshot.embedded.all.sds, "data/organoid_slingshot_sds_embedded_all_no_stretch.rds")

## Plot pseudotime values for cells in each lineage
no.cols <- 2
pseudotime <- slingPseudotime(organoid.slingshot.sce)
names <- colnames(pseudotime)
no.rows <- ceiling(length(names)/no.cols)
pal <- viridis(100)
pdf("plots/slingshot/organoid_slingshot_pseudotime_lineages_no_stretch.pdf", width = 8.3, height = 8.3)
png("plots/slingshot/organoid_slingshot_pseudotime_lineages_no_stretch.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
par(mfrow = c(no.rows, no.cols))
for (i in names){
  cols <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(organoid.slingshot.sce, "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
  lines(slingshot.embedded.all.sds,
        lwd = 1, col = "black", type = "lineages", cex = 1)
}
dev.off()

## Embed lineages one by one
slingshot.embedded <- embedCurves(organoid.slingshot.sce, "UMAP")

lineage.names <- colnames(slingshot.embedded)
umap.clust.gg <- plotUMAP(organoid.slingshot.sce, colour_by = "Cluster_Name") +
  umap.theme() + scale_colour_manual(values = identity.cols) +
  theme(text = element_text(size = 15)) + guides(colour = guide_legend(override.aes = list(size = 5)))
umap.amt.gg <- plotUMAP(organoid.slingshot.sce, colour_by = "AMT.score") +
  umap.theme() + scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))

for (i in 1:length(lineage.names)){
  embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
  embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
  saveRDS(embedded.lineage, paste0("data/organoid_slingshot_embedded_no_stretch_", lineage.names[[i]], ".rds"))
  
  plotUMAP(organoid.slingshot.sce, colour_by = paste0("slingPseudotime_", i)) +
    geom_path(data = embedded.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]]) + umap.theme() +
    theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/organoid_slingshot_no_stretch_", lineage.names[[i]], ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/organoid_slingshot_no_stretch_", lineage.names[[i]], ".png"), width = 8.3, height = 5.8)
  
  umap.clust.gg +
    geom_path(data = embedded.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/organoid_slingshot_no_stretch_", lineage.names[[i]], "_cluster.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/organoid_slingshot_no_stretch_", lineage.names[[i]], "_cluster.png"), width = 8.3, height = 5.8)
  
  umap.amt.gg +
    geom_path(data = embedded.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/organoid_slingshot_no_stretch_", lineage.names[[i]], "_AMT.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/organoid_slingshot_no_stretch_", lineage.names[[i]], "_AMT.png"), width = 8.3, height = 5.8)
}

## Lineages of interest for downstream analysis:
## Lineage 1 = SYMP-SYMP
## Lineage 2 = SYMP-MES
## Lineage 3 = SYMP-ADRN
## Lineage 4 = SYMP-SYMP

amt2.lineage <- readRDS("data/organoid_slingshot_embedded_no_stretch_Lineage2.rds")

amt2.gg <- umap.clust.gg +
  geom_path(data = amt2.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
  labs(title = "Lineage 2:\nSYMP to MES", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt2.gg

ggsave("plots/slingshot/annotated/organoid_umap_slingshot_amt2.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/organoid_umap_slingshot_amt2.png", width = 8.3, height = 5.8)

## Correlation between AMT score and pseudotime
organoid.slingshot.sce <- readRDS("data/organoid_sce_slingshot_no_stretch.rds")

slingshot.df <- as.data.frame(colData(organoid.slingshot.sce)[, c("slingPseudotime_2", "AMT.score")])

lin2.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_2, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 35, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "Lineage 2",
       colour = "AMT score")

lin2.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/annotated/organoid_scatter_amt2.pdf", width = 5.8, height = 5.8)
ggsave("plots/slingshot/annotated/organoid_scatter_amt2.png", width = 5.8, height = 5.8)

## [ Separated trajectory inference ] ----

## NB039
nb039.sce <- readRDS("data/nb039_sce_start.rds")

nb039.slingshot.sce <- slingshot(nb039.sce, reducedDim = "PCA",
                                 clusterLabels = "seurat_clusters", start.clus = 0,
                                 stretch = 0)

saveRDS(nb039.slingshot.sce, "data/nb039_sce_slingshot_no_stretch.rds")

nb039.slingshot.sce <- readRDS("data/nb039_sce_slingshot_no_stretch.rds")

## Extract matrix of pseudotime values along each lineage
slingshot.lineages <- slingPseudotime(nb039.slingshot.sce)
head(slingshot.lineages)
## 3 lineages

## Single pseudotime for all cells
shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)

umap.gg <- plotUMAP(nb039.slingshot.sce, colour_by = I(shared.pseudotime)) +
  scale_colour_viridis() +
  labs(colour = "Pseudotime")
slingshot.embedded.all <- embedCurves(nb039.slingshot.sce, "UMAP")
slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
for (lineage in slingshot.embedded.all) {
  slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
  umap.gg <- umap.gg + geom_path(data = slingshot.embedded.all,
                                 aes(x = umap_1, y = umap_2), linewidth = 1.2)
}

saveRDS(slingshot.embedded.all, "data/nb039_slingshot_embedded_all_no_stretch.rds")

umap.gg + theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/nb039_umap_slingshot_no_stretch.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/nb039_umap_slingshot_no_stretch.png", width = 8.3, height = 5.8)

## Save as SDS
slingshot.embedded.all.sds <- embedCurves(nb039.slingshot.sce, "UMAP")
slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)
saveRDS(slingshot.embedded.all.sds, "data/nb039_slingshot_sds_embedded_all_no_stretch.rds")

## Plot pseudotime values for cells in each lineage
no.cols <- 2
pseudotime <- slingPseudotime(nb039.slingshot.sce)
names <- colnames(pseudotime)
no.rows <- ceiling(length(names)/no.cols)
pal <- viridis(100)
pdf("plots/slingshot/nb039_slingshot_pseudotime_lineages_no_stretch.pdf", width = 8.3, height = 8.3)
png("plots/slingshot/nb039_slingshot_pseudotime_lineages_no_stretch.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
par(mfrow = c(no.rows, no.cols))
for (i in names){
  cols <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(nb039.slingshot.sce, "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
  lines(slingshot.embedded.all.sds,
        lwd = 1, col = "black", type = "lineages", cex = 1)
}
dev.off()

## Embed lineages one by one
slingshot.embedded <- embedCurves(nb039.slingshot.sce, "UMAP")

lineage.names <- colnames(slingshot.embedded)
umap.clust.gg <- plotUMAP(nb039.slingshot.sce, colour_by = "Cluster_Name") +
  umap.theme() + scale_colour_manual(values = identity.cols) +
  theme(text = element_text(size = 15)) + guides(colour = guide_legend(override.aes = list(size = 5)))
umap.amt.gg <- plotUMAP(nb039.slingshot.sce, colour_by = "AMT.score") +
  umap.theme() + scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))

for (i in 1:length(lineage.names)){
  embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
  embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
  saveRDS(embedded.lineage, paste0("data/nb039_slingshot_embedded_no_stretch_", lineage.names[[i]], ".rds"))
  
  plotUMAP(nb039.slingshot.sce, colour_by = paste0("slingPseudotime_", i)) +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]]) + umap.theme() +
    theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb039_slingshot_no_stretch_", lineage.names[[i]], ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb039_slingshot_no_stretch_", lineage.names[[i]], ".png"), width = 8.3, height = 5.8)
  
  umap.clust.gg +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb039_slingshot_no_stretch_", lineage.names[[i]], "_cluster.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb039_slingshot_no_stretch_", lineage.names[[i]], "_cluster.png"), width = 8.3, height = 5.8)
  
  umap.amt.gg +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb039_slingshot_no_stretch_", lineage.names[[i]], "_AMT.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb039_slingshot_no_stretch_", lineage.names[[i]], "_AMT.png"), width = 8.3, height = 5.8)
}

## Lineages of interest for downstream analysis:
## Lineage 1 = SYMP-MES
## Lineage 2 = SYMP-MES?
## Lineage 3 = SYMP-SYMP

amt1.lineage <- readRDS("data/nb039_slingshot_embedded_no_stretch_Lineage1.rds")

amt1.gg <- umap.clust.gg +
  geom_path(data = amt1.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
  labs(title = "Lineage 1:\nSYMP to MES", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt1.gg

ggsave("plots/slingshot/annotated/nb039_umap_slingshot_amt1.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/nb039_umap_slingshot_amt1.png", width = 8.3, height = 5.8)

amt2.lineage <- readRDS("data/nb039_slingshot_embedded_no_stretch_Lineage2.rds")

amt2.gg <- umap.clust.gg +
  geom_path(data = amt2.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
  labs(title = "Lineage 2:\nSYMP to MES?", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt2.gg

ggsave("plots/slingshot/annotated/nb039_umap_slingshot_amt2.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/nb039_umap_slingshot_amt2.png", width = 8.3, height = 5.8)

## Correlation between AMT score and pseudotime
nb039.slingshot.sce <- readRDS("data/nb039_sce_slingshot_no_stretch.rds")

slingshot.df <- as.data.frame(colData(nb039.slingshot.sce)[, c("slingPseudotime_1", "slingPseudotime_2", "AMT.score")])

lin1.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_1, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 50, label.y = -15) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "NB039 Lineage 1",
       colour = "AMT score")

lin1.gg + 
  theme(text = element_text(size = 20, colour = "black"), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/annotated/nb039_scatter_amt1_square.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/nb039_scatter_amt1_square.png", width = 8.3, height = 5.8)

lin2.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_2, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 50, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "Lineage 2",
       colour = "AMT score")

lin2.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/annotated/nb039_scatter_amt2.pdf", width = 5.8, height = 5.8)
ggsave("plots/slingshot/annotated/nb039_scatter_amt2.png", width = 5.8, height = 5.8)

## NB067
nb067.sce <- readRDS("data/nb067_sce_start.rds")

nb067.slingshot.sce <- slingshot(nb067.sce, reducedDim = "PCA",
                                 clusterLabels = "seurat_clusters", start.clus = 0,
                                 stretch = 0)

saveRDS(nb067.slingshot.sce, "data/nb067_sce_slingshot_no_stretch.rds")

nb067.slingshot.sce <- readRDS("data/nb067_sce_slingshot_no_stretch.rds")

## Extract matrix of pseudotime values along each lineage
slingshot.lineages <- slingPseudotime(nb067.slingshot.sce)
head(slingshot.lineages)
## 2 lineages

## Single pseudotime for all cells
shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)

umap.gg <- plotUMAP(nb067.slingshot.sce, colour_by = I(shared.pseudotime)) +
  scale_colour_viridis() +
  labs(colour = "Pseudotime")
slingshot.embedded.all <- embedCurves(nb067.slingshot.sce, "UMAP")
slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
for (lineage in slingshot.embedded.all) {
  slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
  umap.gg <- umap.gg + geom_path(data = slingshot.embedded.all,
                                 aes(x = umap_1, y = umap_2), linewidth = 1.2)
}

saveRDS(slingshot.embedded.all, "data/nb067_slingshot_embedded_all_no_stretch.rds")

umap.gg + theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/nb067_umap_slingshot_no_stretch.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/nb067_umap_slingshot_no_stretch.png", width = 8.3, height = 5.8)

## Save as SDS
slingshot.embedded.all.sds <- embedCurves(nb067.slingshot.sce, "UMAP")
slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)
saveRDS(slingshot.embedded.all.sds, "data/nb067_slingshot_sds_embedded_all_no_stretch.rds")

## Plot pseudotime values for cells in each lineage
no.cols <- 2
pseudotime <- slingPseudotime(nb067.slingshot.sce)
names <- colnames(pseudotime)
no.rows <- ceiling(length(names)/no.cols)
pal <- viridis(100)
pdf("plots/slingshot/nb067_slingshot_pseudotime_lineages_no_stretch.pdf", width = 8.3, height = 5.8)
png("plots/slingshot/nb067_slingshot_pseudotime_lineages_no_stretch.png", width = 8.3, height = 5.8,
    units = "in", res = 200)
par(mfrow = c(no.rows, no.cols))
for (i in names){
  cols <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(nb067.slingshot.sce, "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
  lines(slingshot.embedded.all.sds,
        lwd = 1, col = "black", type = "lineages", cex = 1)
}
dev.off()

## Embed lineages one by one
slingshot.embedded <- embedCurves(nb067.slingshot.sce, "UMAP")

lineage.names <- colnames(slingshot.embedded)
umap.clust.gg <- plotUMAP(nb067.slingshot.sce, colour_by = "Cluster_Name") +
  umap.theme() + scale_colour_manual(values = identity.cols) +
  theme(text = element_text(size = 15)) + guides(colour = guide_legend(override.aes = list(size = 5)))
umap.amt.gg <- plotUMAP(nb067.slingshot.sce, colour_by = "AMT.score") +
  umap.theme() + scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))

for (i in 1:length(lineage.names)){
  embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
  embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
  saveRDS(embedded.lineage, paste0("data/nb067_slingshot_embedded_no_stretch_", lineage.names[[i]], ".rds"))
  
  plotUMAP(nb067.slingshot.sce, colour_by = paste0("slingPseudotime_", i)) +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]]) + umap.theme() +
    theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb067_slingshot_no_stretch_", lineage.names[[i]], ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb067_slingshot_no_stretch_", lineage.names[[i]], ".png"), width = 8.3, height = 5.8)
  
  umap.clust.gg +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb067_slingshot_no_stretch_", lineage.names[[i]], "_cluster.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb067_slingshot_no_stretch_", lineage.names[[i]], "_cluster.png"), width = 8.3, height = 5.8)
  
  umap.amt.gg +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb067_slingshot_no_stretch_", lineage.names[[i]], "_AMT.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/nb067_slingshot_no_stretch_", lineage.names[[i]], "_AMT.png"), width = 8.3, height = 5.8)
}

## Lineages of interest for downstream analysis:
## Lineage 1 = SYMP-MES
## Lineage 2 = SYMP-ADRN

amt1.lineage <- readRDS("data/nb067_slingshot_embedded_no_stretch_Lineage1.rds")

amt1.gg <- umap.clust.gg +
  geom_path(data = amt1.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
  labs(title = "Lineage 1:\nSYMP to MES", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt1.gg

ggsave("plots/slingshot/annotated/nb067_umap_slingshot_amt1.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/nb067_umap_slingshot_amt1.png", width = 8.3, height = 5.8)

## Correlation between AMT score and pseudotime
nb067.slingshot.sce <- readRDS("data/nb067_sce_slingshot_no_stretch.rds")

slingshot.df <- as.data.frame(colData(nb067.slingshot.sce)[, c("slingPseudotime_1", "AMT.score")])

lin1.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_1, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 40, label.y = -14) +
  theme_bw() +
  theme(aspect.ratio = 1,        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "NB067 Lineage 1",
       colour = "AMT score")

lin1.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/annotated/nb067_scatter_amt1_square.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/nb067_scatter_amt1_square.png", width = 8.3, height = 5.8)

## [ Combined gene dynamics ] ----

organoid.slingshot.sce <- readRDS("data/organoid_sce_slingshot_no_stretch.rds")

## Linear model fitting for lineage 2
lineage2 <- testPseudotime(organoid.slingshot.sce, pseudotime = organoid.slingshot.sce$slingPseudotime_2)
lineage2.genes <- lineage2[order(lineage2$p.value), ]
head(lineage2.genes, 10)

dir.create("data/organoid_dynamic_genes/combined/full_lists", recursive = TRUE)

saveRDS(lineage2.genes, "data/organoid_dynamic_genes/combined/full_lists/lineage2_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 2
lineage2.genes.down <- as.data.frame(lineage2.genes[lineage2.genes$logFC < 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage2.genes.down, 10)

saveRDS(lineage2.genes.down, "data/organoid_dynamic_genes/combined/lineage2_down_genes.rds")
write.csv(lineage2.genes.down,"data/organoid_dynamic_genes/combined/lineage2_down_genes.csv", row.names = FALSE)

lineage2.genes.down <- readRDS("data/organoid_dynamic_genes/combined/lineage2_down_genes.rds")

lin2.down.gg <- plotExpression(organoid.slingshot.sce, features = head(lineage2.genes.down$gene, 10), x = "slingPseudotime_2", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 2: downregulated")

lin2.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/organoid_amt2_down_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/organoid_amt2_down_top10.png", width = 11.7, height = 8.3)

on.lineage2 <- !is.na(organoid.slingshot.sce$slingPseudotime_2)
lin2.down.heatmap <- plotHeatmap(organoid.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage2.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin2.down.heatmap$gtable$grobs[[6]]$children$GRID.text.7558$label <- "Pseudotime"
lin2.down.heatmap$gtable$grobs[[6]]$children$GRID.text.7562$label <- "Cluster Name"

lin2.down.heatmap

pdf("plots/gene_dynamics/organoid_amt2_down_top50.pdf", width = 8.3, height = 8.3)
lin2.down.heatmap
dev.off()

png("plots/gene_dynamics/organoid_amt2_down_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin2.down.heatmap
dev.off()

# Genes that increase in expression with increasing pseudotime along lineage 2
lineage2.genes.up <- as.data.frame(lineage2.genes[lineage2.genes$logFC > 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage2.genes.up, 10)

saveRDS(lineage2.genes.up, "data/organoid_dynamic_genes/combined/lineage2_up_genes.rds")
write.csv(lineage2.genes.up,"data/organoid_dynamic_genes/combined/lineage2_up_genes.csv", row.names = FALSE)

lineage2.genes.up <- readRDS("data/organoid_dynamic_genes/combined/lineage2_up_genes.rds")

lin2.up.gg <- plotExpression(organoid.slingshot.sce, features = head(lineage2.genes.up$gene, 10), x = "slingPseudotime_2", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 2: upregulated")

lin2.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/organoid_amt2_up_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/organoid_amt2_up_top10.png", width = 11.7, height = 8.3)

lin2.up.gg + lin2.down.gg + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Lineage 2") &
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/organoid_amt2_lineage_expression.pdf", width = 16.5, height = 8.3)
ggsave("plots/gene_dynamics/organoid_amt2_lineage_expression.png", width = 16.5, height = 8.3)

lin2.up.heatmap <- plotHeatmap(organoid.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage2.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin2.up.heatmap$gtable$grobs[[6]]$children$GRID.text.12131$label <- "Pseudotime"
lin2.up.heatmap$gtable$grobs[[6]]$children$GRID.text.12135$label <- "Cluster Name"

lin2.up.heatmap

pdf("plots/gene_dynamics/organoid_amt2_up_top50.pdf", width = 8.3, height = 8.3)
lin2.up.heatmap
dev.off()

png("plots/gene_dynamics/organoid_amt2_up_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin2.up.heatmap
dev.off()

## [ Separated gene dynamics ] ----

## NB039
nb039.slingshot.sce <- readRDS("data/nb039_sce_slingshot_no_stretch.rds")

## Linear model fitting for lineage 1
lineage1 <- testPseudotime(nb039.slingshot.sce, pseudotime = nb039.slingshot.sce$slingPseudotime_1)
lineage1.genes <- lineage1[order(lineage1$p.value), ]
head(lineage1.genes, 10)

dir.create("data/organoid_dynamic_genes/nb039/full_lists", recursive = TRUE)

saveRDS(lineage1.genes, "data/organoid_dynamic_genes/nb039/full_lists/lineage1_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 1
lineage1.genes.down <- as.data.frame(lineage1.genes[lineage1.genes$logFC < 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage1.genes.down, 10)

saveRDS(lineage1.genes.down, "data/organoid_dynamic_genes/nb039/lineage1_down_genes.rds")
write.csv(lineage1.genes.down,"data/organoid_dynamic_genes/nb039/lineage1_down_genes.csv", row.names = FALSE)

lineage1.genes.down <- readRDS("data/organoid_dynamic_genes/nb039/lineage1_down_genes.rds")

lin1.down.gg <- plotExpression(nb039.slingshot.sce, features = head(lineage1.genes.down$gene, 10), x = "slingPseudotime_1", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 1: downregulated")

lin1.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/nb039_amt1_down_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/nb039_amt1_down_top10.png", width = 11.7, height = 8.3)

on.lineage1 <- !is.na(nb039.slingshot.sce$slingPseudotime_1)
lin1.down.heatmap <- plotHeatmap(nb039.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage1.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin1.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin1.down.heatmap$gtable$grobs[[6]]$children$GRID.text.14631$label <- "Pseudotime"
lin1.down.heatmap$gtable$grobs[[6]]$children$GRID.text.14635$label <- "Cluster Name"

lin1.down.heatmap

pdf("plots/gene_dynamics/nb039_amt1_down_top50.pdf", width = 8.3, height = 8.3)
lin1.down.heatmap
dev.off()

png("plots/gene_dynamics/nb039_amt1_down_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin1.down.heatmap
dev.off()

## Genes that increase in expression with increasing pseudotime along lineage 1
lineage1.genes.up <- as.data.frame(lineage1.genes[lineage1.genes$logFC > 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage1.genes.up, 10)

saveRDS(lineage1.genes.up, "data/organoid_dynamic_genes/nb039/lineage1_up_genes.rds")
write.csv(lineage1.genes.up,"data/organoid_dynamic_genes/nb039/lineage1_up_genes.csv", row.names = FALSE)

lineage1.genes.up <- readRDS("data/organoid_dynamic_genes/nb039/lineage1_up_genes.rds")

lin1.up.gg <- plotExpression(nb039.slingshot.sce, features = head(lineage1.genes.up$gene, 10), x = "slingPseudotime_1", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 1: upregulated")

lin1.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/nb039_amt1_up_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/nb039_amt1_up_top10.png", width = 11.7, height = 8.3)

lin1.up.gg + lin1.down.gg + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Lineage 1") &
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/nb039_amt1_lineage_expression.pdf", width = 16.5, height = 8.3)
ggsave("plots/gene_dynamics/nb039_amt1_lineage_expression.png", width = 16.5, height = 8.3)

lin1.up.heatmap <- plotHeatmap(nb039.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage1.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin1.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin1.up.heatmap$gtable$grobs[[6]]$children$GRID.text.19204$label <- "Pseudotime"
lin1.up.heatmap$gtable$grobs[[6]]$children$GRID.text.19208$label <- "Cluster Name"

lin1.up.heatmap

pdf("plots/gene_dynamics/nb039_amt1_up_top50.pdf", width = 8.3, height = 8.3)
lin1.up.heatmap
dev.off()

png("plots/gene_dynamics/nb039_amt1_up_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin1.up.heatmap
dev.off()

## Linear model fitting for lineage 2
lineage2 <- testPseudotime(nb039.slingshot.sce, pseudotime = nb039.slingshot.sce$slingPseudotime_2)
lineage2.genes <- lineage1[order(lineage2$p.value), ]
head(lineage2.genes, 10)

saveRDS(lineage2.genes, "data/organoid_dynamic_genes/nb039/full_lists/lineage2_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 2
lineage2.genes.down <- as.data.frame(lineage2.genes[lineage2.genes$logFC < 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage2.genes.down, 10)

saveRDS(lineage2.genes.down, "data/organoid_dynamic_genes/nb039/lineage2_down_genes.rds")
write.csv(lineage2.genes.down,"data/organoid_dynamic_genes/nb039/lineage2_down_genes.csv", row.names = FALSE)

lineage2.genes.down <- readRDS("data/organoid_dynamic_genes/nb039/lineage2_down_genes.rds")

lin2.down.gg <- plotExpression(nb039.slingshot.sce, features = head(lineage2.genes.down$gene, 10), x = "slingPseudotime_2", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 2: downregulated")

lin2.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/nb039_amt2_down_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/nb039_amt2_down_top10.png", width = 11.7, height = 8.3)

on.lineage2 <- !is.na(nb039.slingshot.sce$slingPseudotime_2)
lin2.down.heatmap <- plotHeatmap(nb039.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage2.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin2.down.heatmap$gtable$grobs[[6]]$children$GRID.text.20726$label <- "Pseudotime"
lin2.down.heatmap$gtable$grobs[[6]]$children$GRID.text.20730$label <- "Cluster Name"

lin2.down.heatmap

pdf("plots/gene_dynamics/nb039_amt2_down_top50.pdf", width = 8.3, height = 8.3)
lin2.down.heatmap
dev.off()

png("plots/gene_dynamics/nb039_amt2_down_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin2.down.heatmap
dev.off()

## Genes that increase in expression with increasing pseudotime along lineage 2
lineage2.genes.up <- as.data.frame(lineage2.genes[lineage2.genes$logFC > 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage2.genes.up, 10)

saveRDS(lineage2.genes.up, "data/organoid_dynamic_genes/nb039/lineage2_up_genes.rds")
write.csv(lineage2.genes.up,"data/organoid_dynamic_genes/nb039/lineage2_up_genes.csv", row.names = FALSE)

lineage2.genes.up <- readRDS("data/organoid_dynamic_genes/nb039/lineage2_up_genes.rds")

lin2.up.gg <- plotExpression(nb039.slingshot.sce, features = head(lineage2.genes.up$gene, 10), x = "slingPseudotime_2", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 2: upregulated")

lin2.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/nb039_amt2_up_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/nb039_amt2_up_top10.png", width = 11.7, height = 8.3)

lin2.up.gg + lin2.down.gg + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Lineage 2") &
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/nb039_amt2_lineage_expression.pdf", width = 16.5, height = 8.3)
ggsave("plots/gene_dynamics/nb039_amt2_lineage_expression.png", width = 16.5, height = 8.3)

lin2.up.heatmap <- plotHeatmap(nb039.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage2.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin2.up.heatmap$gtable$grobs[[6]]$children$GRID.text.25299$label <- "Pseudotime"
lin2.up.heatmap$gtable$grobs[[6]]$children$GRID.text.25303$label <- "Cluster Name"

lin2.up.heatmap

pdf("plots/gene_dynamics/nb039_amt2_up_top50.pdf", width = 8.3, height = 8.3)
lin2.up.heatmap
dev.off()

png("plots/gene_dynamics/nb039_amt2_up_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin2.up.heatmap
dev.off()

## NB067
nb067.slingshot.sce <- readRDS("data/nb067_sce_slingshot_no_stretch.rds")

## Linear model fitting for lineage 1
lineage1 <- testPseudotime(nb067.slingshot.sce, pseudotime = nb067.slingshot.sce$slingPseudotime_1)
lineage1.genes <- lineage1[order(lineage1$p.value), ]
head(lineage1.genes, 10)

dir.create("data/organoid_dynamic_genes/nb067/full_lists", recursive = TRUE)

saveRDS(lineage1.genes, "data/organoid_dynamic_genes/nb067/full_lists/lineage1_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 1
lineage1.genes.down <- as.data.frame(lineage1.genes[lineage1.genes$logFC < 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage1.genes.down, 10)

saveRDS(lineage1.genes.down, "data/organoid_dynamic_genes/nb067/lineage1_down_genes.rds")
write.csv(lineage1.genes.down,"data/organoid_dynamic_genes/nb067/lineage1_down_genes.csv", row.names = FALSE)

lineage1.genes.down <- readRDS("data/organoid_dynamic_genes/nb067/lineage1_down_genes.rds")

lin1.down.gg <- plotExpression(nb067.slingshot.sce, features = head(lineage1.genes.down$gene, 10), x = "slingPseudotime_1", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 1: downregulated")

lin1.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/nb067_amt1_down_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/nb067_amt1_down_top10.png", width = 11.7, height = 8.3)

on.lineage1 <- !is.na(nb067.slingshot.sce$slingPseudotime_1)
lin1.down.heatmap <- plotHeatmap(nb067.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage1.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin1.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin1.down.heatmap$gtable$grobs[[6]]$children$GRID.text.27310$label <- "Pseudotime"
lin1.down.heatmap$gtable$grobs[[6]]$children$GRID.text.27314$label <- "Cluster Name"

lin1.down.heatmap

pdf("plots/gene_dynamics/nb067_amt1_down_top50.pdf", width = 8.3, height = 8.3)
lin1.down.heatmap
dev.off()

png("plots/gene_dynamics/nb067_amt1_down_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin1.down.heatmap
dev.off()

## Genes that increase in expression with increasing pseudotime along lineage 1
lineage1.genes.up <- as.data.frame(lineage1.genes[lineage1.genes$logFC > 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage1.genes.up, 10)

saveRDS(lineage1.genes.up, "data/organoid_dynamic_genes/nb067/lineage1_up_genes.rds")
write.csv(lineage1.genes.up,"data/organoid_dynamic_genes/nb067/lineage1_up_genes.csv", row.names = FALSE)

lineage1.genes.up <- readRDS("data/organoid_dynamic_genes/nb067/lineage1_up_genes.rds")

lin1.up.gg <- plotExpression(nb067.slingshot.sce, features = head(lineage1.genes.up$gene, 10), x = "slingPseudotime_1", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 1: upregulated")

lin1.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/nb067_amt1_up_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/nb067_amt1_up_top10.png", width = 11.7, height = 8.3)

lin1.up.gg + lin1.down.gg + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Lineage 1") &
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/nb067_amt1_lineage_expression.pdf", width = 16.5, height = 8.3)
ggsave("plots/gene_dynamics/nb067_amt1_lineage_expression.png", width = 16.5, height = 8.3)

lin1.up.heatmap <- plotHeatmap(nb067.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage1.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin1.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin1.up.heatmap$gtable$grobs[[6]]$children$GRID.text.31883$label <- "Pseudotime"
lin1.up.heatmap$gtable$grobs[[6]]$children$GRID.text.31887$label <- "Cluster Name"

lin1.up.heatmap

pdf("plots/gene_dynamics/nb067_amt1_up_top50.pdf", width = 8.3, height = 8.3)
lin1.up.heatmap
dev.off()

png("plots/gene_dynamics/nb067_amt1_up_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin1.up.heatmap
dev.off()

## [ Gene list comparison ] ----

## Plot all comparisons as heatmap
## PDX
pdx.lin1.down <- readRDS("data/pdx_dynamic_genes/lineage1_down_genes.rds")
pdx.lin1.up <- readRDS("data/pdx_dynamic_genes/lineage1_up_genes.rds")

pdx.lin2.down <- readRDS("data/pdx_dynamic_genes/lineage2_down_genes.rds")
pdx.lin2.up <- readRDS("data/pdx_dynamic_genes/lineage2_up_genes.rds")

pdx.lin3.down <- readRDS("data/pdx_dynamic_genes/lineage3_down_genes.rds")
pdx.lin3.up <- readRDS("data/pdx_dynamic_genes/lineage3_up_genes.rds")

## Combined
organoid.lin2.down <- readRDS("data/organoid_dynamic_genes/combined/lineage2_down_genes.rds")
organoid.lin2.up <- readRDS("data/organoid_dynamic_genes/combined/lineage2_up_genes.rds")

## NB039 only
nb039.lin1.down <- readRDS("data/organoid_dynamic_genes/nb039/lineage1_down_genes.rds")
nb039.lin1.up <- readRDS("data/organoid_dynamic_genes/nb039/lineage1_up_genes.rds")

nb039.lin2.down <- readRDS("data/organoid_dynamic_genes/nb039/lineage2_down_genes.rds")
nb039.lin2.up <- readRDS("data/organoid_dynamic_genes/nb039/lineage2_up_genes.rds")

## NB067 only
nb067.lin1.down <- readRDS("data/organoid_dynamic_genes/nb067/lineage1_down_genes.rds")
nb067.lin1.up <- readRDS("data/organoid_dynamic_genes/nb067/lineage1_up_genes.rds")

## Reference from SK-N-SH barcoding experiment
ref.lin1.down <- readRDS("data/ref_dynamic_genes/lineage1_down_genes.rds")
ref.lin1.up <- readRDS("data/ref_dynamic_genes/lineage1_up_genes.rds")

ref.lin2.down <- readRDS("data/ref_dynamic_genes/lineage2_down_genes.rds")
ref.lin2.up <- readRDS("data/ref_dynamic_genes/lineage2_up_genes.rds")

down.genes <- list(NB039_1 = nb039.lin1.down$gene,
                   NB067_1 = nb067.lin1.down$gene,
                   PDX_1 = pdx.lin1.down$gene,
                   Ref_1 = ref.lin1.down$gene,
                   Ref_2 = ref.lin2.down$gene)
up.genes <- list(NB039_1 = nb039.lin1.up$gene,
                 NB067_1 = nb067.lin1.up$gene,
                 PDX_1 = pdx.lin1.up$gene,
                 Ref_1 = ref.lin1.up$gene,
                 Ref_2 = ref.lin2.up$gene)

down.mat <- crossprod(table(stack(down.genes)))

require(reshape2)
down.df <- melt(as.matrix(down.mat))
colnames(down.df) <- c("var1", "var2", "value")

down.gg <- ggplot(down.df, aes(var1, var2, fill = value)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "#0077BB") +
  labs(title = "Downregulated\nplasticity genes", x = "", y = "") +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_blank(),
        aspect.ratio = 1)

down.gg

ggsave("plots/gene_dynamics/down_genes_comparison_heatmap_manuscript.pdf", width = 5.8, height = 5.8)
ggsave("plots/gene_dynamics/down_genes_comparison_heatmap_manuscript.png", width = 5.8, height = 5.8)

up.mat <- crossprod(table(stack(up.genes)))

up.df <- melt(as.matrix(up.mat))
colnames(up.df) <- c("var1", "var2", "value")

up.gg <- ggplot(up.df, aes(var1, var2, fill = value)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "#0077BB") +
  labs(title = "Upregulated\nplasticity genes", x = "", y = "") +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.title = element_blank(),
        aspect.ratio = 1)

up.gg

ggsave("plots/gene_dynamics/up_genes_comparison_heatmap_manuscript.pdf", width = 5.8, height = 5.8)
ggsave("plots/gene_dynamics/up_genes_comparison_heatmap_manuscript.png", width = 5.8, height = 5.8)
