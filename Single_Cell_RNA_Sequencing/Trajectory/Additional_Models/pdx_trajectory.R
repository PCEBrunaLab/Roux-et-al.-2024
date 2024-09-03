## NB PDX trajectory analysis

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

identity.cols <- c("MES Myc 1" = "#CCAA11",
                   "MES TGFb" = "#225555",
                   "MES Myc 2" = "#117733",
                   "MES-INT" = "#CC3311",
                   "ADRN cycling" = "#882255",
                   "ADRN" = "#BB22EE",
                   "MES cycling" = "#EE7733",
                   "MES ROS" = "#0077BB")

## [ PDX signature plots ] ----

pdx.seurat <- readRDS("data/nb_seurat_AMT_non-batch_pdx.rds")

## Cluster colours
#cluster.cols <- iwanthue(length(unique(pdx.seurat$seurat_clusters)))
#saveRDS(cluster.cols, "data/pdx_cluster_cols.rds")
cluster.cols <- readRDS("data/pdx_cluster_cols.rds")

num.gg <- DimPlot(pdx.seurat,
                  group.by = "seurat_clusters", order = TRUE) +
  umap.theme() + labs(title = "Cluster number") +
  scale_colour_manual(values = cluster.cols)

num.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/pdx_umap_cluster_number.pdf", width = 5.8, height = 5.8)
ggsave("plots/pdx_umap_cluster_number.png", width = 5.8, height = 5.8)

dir.create("plots/signatures/", recursive = TRUE)

adrn.gg <- FeaturePlot(pdx.seurat,
                       features = "ADRN.Sig",
                       order = TRUE,
                       min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.matter(100)) +
  labs(title = "vanGroningen\nADRN signature", colour = "Expression") +
  umap.theme()

mes.gg <- FeaturePlot(pdx.seurat,
                      features = "MES.Sig",
                      order = TRUE,
                      min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.haline(100)) +
  labs(title = "vanGroningen\nMES signature", colour = "Expression") +
  umap.theme()

adrn.gg + mes.gg &
  theme(text = element_text(size = 15),
        legend.key.size = unit(1, "cm"))

ggsave("plots/signatures/pdx_vG_ADRN_MES_umaps.pdf", width = 5.8, height = 8.3)
ggsave("plots/signatures/pdx_vG_ADRN_MES_umaps.png", width = 5.8, height = 8.3)

vln1.gg <- VlnPlot(pdx.seurat,
                   features = "MES.Sig", group.by = "seurat_clusters") +
  scale_fill_manual(values = cluster.cols) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln1.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln1.gg$layers[[2]] <- NULL

vln2.gg <- VlnPlot(pdx.seurat,
                   features = "ADRN.Sig", group.by = "seurat_clusters") +
  scale_fill_manual(values = cluster.cols) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln2.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln2.gg$layers[[2]] <- NULL

vln3.gg <- VlnPlot(pdx.seurat,
                   features = "AMT.score", group.by = "seurat_clusters") +
  scale_fill_manual(values = cluster.cols) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln3.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln3.gg$layers[[2]] <- NULL

(vln3.gg | (vln1.gg / vln2.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/pdx_vG_ADRN_MES_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/pdx_vG_ADRN_MES_violins.png", width = 8.3, height = 8.3)

amt.gg <- DimPlot(pdx.seurat,
                  group.by = "AMT.state", order = TRUE) +
  umap.theme() + labs(title = "AMT state") +
  scale_colour_manual(values = c("#990099", "#F37735", "lightgrey"),
                      labels = c("ADRN", "MES", "Intermediate"))

amt.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/pdx_umap_amt_state.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/pdx_umap_amt_state.png", width = 5.8, height = 5.8)

pdx.seurat$Cluster_Name <- recode(as.character(pdx.seurat$seurat_clusters),
                                  "0" = "MES Myc 1",
                                  "1" = "MES TGFb",
                                  "2" = "MES Myc 2",
                                  "3" = "MES-INT",
                                  "4" = "ADRN cycling",
                                  "5" = "ADRN",
                                  "6" = "MES cycling",
                                  "7" = "MES ROS")

cluster.gg <- DimPlot(pdx.seurat,
                      group.by = "Cluster_Name", order = TRUE) +
  umap.theme() + labs(title = "Initial Cluster Annotation") +
  scale_colour_manual(values = identity.cols)

cluster.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/pdx_umap_initial_clusters.pdf", width = 5.8, height = 5.8)
ggsave("plots/pdx_umap_initial_clusters.png", width = 5.8, height = 5.8)

## [ PDX coexpression ] ----

## Which cells coexpress van Groningen ADRN and MES signatures
coexprs = function(x){
  if (x[["MES.Sig"]] > 0.2 & x[["ADRN.Sig"]] > 0.1) {
    result <- 1
  } else {
    result <- 0
  }
  return(result)
}

metadata.df <- subset(pdx.seurat@meta.data, select = c("MES.Sig", "ADRN.Sig"))

coexprs.df <- apply(metadata.df, 1, coexprs)
pdx.seurat@meta.data <- cbind(pdx.seurat@meta.data, Coexprs.Sig = coexprs.df)

dir.create("plots/coexpression", recursive = TRUE)

coexprs.gg <- DimPlot(pdx.seurat,
                      group.by = "Coexprs.Sig", order = TRUE) +
  umap.theme() + labs(title = "ADRN and MES Coexpression\n(MES > 0.2, ADRN > 0.1)") +
  scale_colour_manual(values = c("#4477AA", "#EE6677"), labels = c("FALSE", "TRUE")) + 
  labs(colour = "Coexpression")

coexprs.gg

ggsave("plots/coexpression/pdx_umap_coexpress.pdf", width = 5.8, height = 5.8)
ggsave("plots/coexpression/pdx_umap_coexpress.png", width = 5.8, height = 5.8)

saveRDS(pdx.seurat, "data/pdx_seurat_start.rds")

## [ PDX trajectory inference ] ----

pdx.sce <- readRDS("data/pdx_sce_start.rds")

pdx.slingshot.sce <- slingshot(pdx.sce, reducedDim = "PCA",
                               clusterLabels = "seurat_clusters", start.clus = 4,
                               stretch = 0)

saveRDS(pdx.slingshot.sce, "data/pdx_sce_slingshot_no_stretch.rds")

pdx.slingshot.sce <- readRDS("data/pdx_sce_slingshot_no_stretch.rds")

## Extract matrix of pseudotime values along each lineage
slingshot.lineages <- slingPseudotime(pdx.slingshot.sce)
head(slingshot.lineages)
## 4 lineages

## Single pseudotime for all cells
shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)

dir.create("plots/slingshot/", recursive = TRUE)

umap.gg <- plotUMAP(pdx.slingshot.sce, colour_by = I(shared.pseudotime)) +
  scale_colour_viridis() +
  labs(colour = "Pseudotime")
slingshot.embedded.all <- embedCurves(pdx.slingshot.sce, "UMAP")
slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
for (lineage in slingshot.embedded.all) {
  slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
  umap.gg <- umap.gg + geom_path(data = slingshot.embedded.all,
                                 aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2)
}

saveRDS(slingshot.embedded.all, "data/pdx_slingshot_embedded_all_no_stretch.rds")

umap.gg + theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/pdx_umap_slingshot_no_stretch.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/pdx_umap_slingshot_no_stretch.png", width = 8.3, height = 5.8)

## Save as SDS
slingshot.embedded.all.sds <- embedCurves(pdx.slingshot.sce, "UMAP")
slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)
saveRDS(slingshot.embedded.all.sds, "data/pdx_slingshot_sds_embedded_all_no_stretch.rds")

## Plot pseudotime values for cells in each lineage
no.cols <- 2
pseudotime <- slingPseudotime(pdx.slingshot.sce)
names <- colnames(pseudotime)
no.rows <- ceiling(length(names)/no.cols)
pal <- viridis(100)
pdf("plots/slingshot/pdx_slingshot_pseudotime_lineages_no_stretch.pdf", width = 8.3, height = 8.3)
png("plots/slingshot/pdx_slingshot_pseudotime_lineages_no_stretch.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
par(mfrow = c(no.rows, no.cols))
for (i in names){
  cols <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(pdx.slingshot.sce, "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
  lines(slingshot.embedded.all.sds,
        lwd = 1, col = "black", type = "lineages", cex = 1)
}
dev.off()

## Embed lineages one by one
slingshot.embedded <- embedCurves(pdx.slingshot.sce, "UMAP")

lineage.names <- colnames(slingshot.embedded)
umap.clust.gg <- plotUMAP(pdx.slingshot.sce, colour_by = "Cluster_Name") +
  umap.theme() + scale_colour_manual(values = identity.cols) +
  theme(text = element_text(size = 15)) + guides(colour = guide_legend(override.aes = list(size = 5)))
umap.amt.gg <- plotUMAP(pdx.slingshot.sce, colour_by = "AMT.score") +
  umap.theme() + scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))

dir.create("plots/slingshot/by_lineage/no_stretch", recursive = TRUE)

for (i in 1:length(lineage.names)){
  embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
  embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
  saveRDS(embedded.lineage, paste0("data/pdx_slingshot_embedded_no_stretch_", lineage.names[[i]], ".rds"))
  
  plotUMAP(pdx.slingshot.sce, colour_by = paste0("slingPseudotime_", i)) +
    geom_path(data = embedded.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]]) + umap.theme() +
    theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/pdx_slingshot_no_stretch_", lineage.names[[i]], ".pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/pdx_slingshot_no_stretch_", lineage.names[[i]], ".png"), width = 8.3, height = 5.8)
  
  umap.clust.gg +
    geom_path(data = embedded.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/pdx_slingshot_no_stretch_", lineage.names[[i]], "_cluster.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/pdx_slingshot_no_stretch_", lineage.names[[i]], "_cluster.png"), width = 8.3, height = 5.8)
  
  umap.amt.gg +
    geom_path(data = embedded.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/pdx_slingshot_no_stretch_", lineage.names[[i]], "_AMT.pdf"), width = 8.3, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/no_stretch/pdx_slingshot_no_stretch_", lineage.names[[i]], "_AMT.png"), width = 8.3, height = 5.8)
}

## Lineages of interest for downstream analysis:
## lineage 1 = ADRN-MES Myc
## Lineage 2 = ADRN-MES TGFb
## Lineage 3 = ADRN-MES Myc
## Lineage 4 = ADRN-ADRN

dir.create("plots/slingshot/annotated", recursive = TRUE)

amt1.lineage <- readRDS("data/pdx_slingshot_embedded_no_stretch_Lineage1.rds")

amt1.gg <- umap.clust.gg +
  geom_path(data = amt1.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
  labs(title = "Lineage 1:\nADRN to Proliferative MES 1", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt1.gg

ggsave("plots/slingshot/annotated/pdx_umap_slingshot_amt1.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/pdx_umap_slingshot_amt1.png", width = 8.3, height = 5.8)

amt2.lineage <- readRDS("data/pdx_slingshot_embedded_no_stretch_Lineage2.rds")

amt2.gg <- umap.clust.gg +
  geom_path(data = amt2.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
  labs(title = "Lineage 2:\nADRN to EMT MES", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt2.gg

ggsave("plots/slingshot/annotated/pdx_umap_slingshot_amt2.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/pdx_umap_slingshot_amt2.png", width = 8.3, height = 5.8)

amt3.lineage <- readRDS("data/pdx_slingshot_embedded_no_stretch_Lineage3.rds")

amt3.gg <- umap.clust.gg +
  geom_path(data = amt3.lineage, aes(x = UMAP_1, y = UMAP_2), linewidth = 1.2) +
  labs(title = "Lineage 3:\nADRN to Proliferative MES 2", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt3.gg

ggsave("plots/slingshot/annotated/pdx_umap_slingshot_amt3.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/pdx_umap_slingshot_amt3.png", width = 8.3, height = 5.8)

## Correlation between AMT score and pseudotime
pdx.slingshot.sce <- readRDS("data/pdx_sce_slingshot_no_stretch.rds")

slingshot.df <- as.data.frame(colData(pdx.slingshot.sce)[, c("slingPseudotime_1", "slingPseudotime_2",
                                                             "slingPseudotime_3", "AMT.score")])
lin1.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_1, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 40, label.y = -15) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "PDX Lineage 1",
       colour = "AMT score")

lin1.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/annotated/pdx_scatter_amt1_square.pdf", width = 8.3, height = 5.8)
ggsave("plots/slingshot/annotated/pdx_scatter_amt1_square.png", width = 8.3, height = 5.8)

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

ggsave("plots/slingshot/annotated/pdx_scatter_amt2.pdf", width = 5.8, height = 5.8)
ggsave("plots/slingshot/annotated/pdx_scatter_amt2.png", width = 5.8, height = 5.8)

lin3.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_3, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 35, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "Lineage 3",
       colour = "AMT score")

lin3.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/annotated/pdx_scatter_amt3.pdf", width = 5.8, height = 5.8)
ggsave("plots/slingshot/annotated/pdx_scatter_amt3.png", width = 5.8, height = 5.8)

## [ PDX gene dynamics ] ----

pdx.slingshot.sce <- readRDS("data/pdx_sce_slingshot_no_stretch.rds")

## Linear model fitting for lineage 1
lineage1 <- testPseudotime(pdx.slingshot.sce, pseudotime = pdx.slingshot.sce$slingPseudotime_1)
lineage1.genes <- lineage1[order(lineage1$p.value), ]
head(lineage1.genes, 10)

dir.create("data/pdx_dynamic_genes/full_lists", recursive = TRUE)

saveRDS(lineage1.genes, "data/pdx_dynamic_genes/full_lists/lineage1_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 1
lineage1.genes.down <- as.data.frame(lineage1.genes[lineage1.genes$logFC < 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage1.genes.down, 10)

saveRDS(lineage1.genes.down, "data/pdx_dynamic_genes/lineage1_down_genes.rds")
write.csv(lineage1.genes.down,"data/pdx_dynamic_genes/lineage1_down_genes.csv", row.names = FALSE)

lineage1.genes.down <- readRDS("data/pdx_dynamic_genes/lineage1_down_genes.rds")

dir.create("plots/gene_dynamics", recursive = TRUE)

lin1.down.gg <- plotExpression(pdx.slingshot.sce, features = head(lineage1.genes.down$gene, 10), x = "slingPseudotime_1", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 1: downregulated")

lin1.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/pdx_amt1_down_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/pdx_amt1_down_top10.png", width = 11.7, height = 8.3)

on.lineage1 <- !is.na(pdx.slingshot.sce$slingPseudotime_1)
lin1.down.heatmap <- plotHeatmap(pdx.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage1.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin1.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin1.down.heatmap$gtable$grobs[[6]]$children$GRID.text.8134$label <- "Pseudotime"
lin1.down.heatmap$gtable$grobs[[6]]$children$GRID.text.8138$label <- "Cluster Name"

lin1.down.heatmap

pdf("plots/gene_dynamics/pdx_amt1_down_top50.pdf", width = 8.3, height = 8.3)
lin1.down.heatmap
dev.off()

png("plots/gene_dynamics/pdx_amt1_down_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin1.down.heatmap
dev.off()

## Genes that increase in expression with increasing pseudotime along lineage 1
lineage1.genes.up <- as.data.frame(lineage1.genes[lineage1.genes$logFC > 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage1.genes.up, 10)

saveRDS(lineage1.genes.up, "data/pdx_dynamic_genes/lineage1_up_genes.rds")
write.csv(lineage1.genes.up,"data/pdx_dynamic_genes/lineage1_up_genes.csv", row.names = FALSE)

lineage1.genes.up <- readRDS("data/pdx_dynamic_genes/lineage1_up_genes.rds")

lin1.up.gg <- plotExpression(pdx.slingshot.sce, features = head(lineage1.genes.up$gene, 10), x = "slingPseudotime_1", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 1: upregulated")

lin1.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/pdx_amt1_up_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/pdx_amt1_up_top10.png", width = 11.7, height = 8.3)

lin1.up.gg + lin1.down.gg + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Lineage 1") &
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/pdx_amt1_lineage_expression.pdf", width = 16.5, height = 8.3)
ggsave("plots/gene_dynamics/pdx_amt1_lineage_expression.png", width = 16.5, height = 8.3)

lin1.up.heatmap <- plotHeatmap(pdx.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage1.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin1.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin1.up.heatmap$gtable$grobs[[6]]$children$GRID.text.12491$label <- "Pseudotime"
lin1.up.heatmap$gtable$grobs[[6]]$children$GRID.text.12495$label <- "Cluster Name"

lin1.up.heatmap

pdf("plots/gene_dynamics/pdx_amt1_up_top50.pdf", width = 8.3, height = 8.3)
lin1.up.heatmap
dev.off()

png("plots/gene_dynamics/pdx_amt1_up_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin1.up.heatmap
dev.off()

## Linear model fitting for lineage 2
lineage2 <- testPseudotime(pdx.slingshot.sce, pseudotime = pdx.slingshot.sce$slingPseudotime_2)
lineage2.genes <- lineage2[order(lineage2$p.value), ]
head(lineage2.genes, 10)

saveRDS(lineage2.genes, "data/pdx_dynamic_genes/full_lists/lineage2_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 2
lineage2.genes.down <- as.data.frame(lineage2.genes[lineage2.genes$logFC < 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage2.genes.down, 10)

saveRDS(lineage2.genes.down, "data/pdx_dynamic_genes/lineage2_down_genes.rds")
write.csv(lineage2.genes.down,"data/pdx_dynamic_genes/lineage2_down_genes.csv", row.names = FALSE)

lineage2.genes.down <- readRDS("data/pdx_dynamic_genes/lineage2_down_genes.rds")

lin2.down.gg <- plotExpression(pdx.slingshot.sce, features = head(lineage2.genes.down$gene, 10), x = "slingPseudotime_2", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 2: downregulated")

lin2.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/pdx_amt2_down_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/pdx_amt2_down_top10.png", width = 11.7, height = 8.3)

on.lineage2 <- !is.na(pdx.slingshot.sce$slingPseudotime_2)
lin2.down.heatmap <- plotHeatmap(pdx.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage2.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin2.down.heatmap$gtable$grobs[[6]]$children$GRID.text.13944$label <- "Pseudotime"
lin2.down.heatmap$gtable$grobs[[6]]$children$GRID.text.13948$label <- "Cluster Name"

lin2.down.heatmap

pdf("plots/gene_dynamics/pdx_amt2_down_top50.pdf", width = 8.3, height = 8.3)
lin2.down.heatmap
dev.off()

png("plots/gene_dynamics/pdx_amt2_down_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin2.down.heatmap
dev.off()

## Genes that increase in expression with increasing pseudotime along lineage 1
lineage2.genes.up <- as.data.frame(lineage2.genes[lineage2.genes$logFC > 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage2.genes.up, 10)

saveRDS(lineage2.genes.up, "data/pdx_dynamic_genes/lineage2_up_genes.rds")
write.csv(lineage2.genes.up,"data/pdx_dynamic_genes/lineage2_up_genes.csv", row.names = FALSE)

lineage2.genes.up <- readRDS("data/pdx_dynamic_genes/lineage2_up_genes.rds")

lin2.up.gg <- plotExpression(pdx.slingshot.sce, features = head(lineage2.genes.up$gene, 10), x = "slingPseudotime_2", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 2: upregulated")

lin2.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/pdx_amt2_up_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/pdx_amt2_up_top10.png", width = 11.7, height = 8.3)

lin2.up.gg + lin2.down.gg + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Lineage 2") &
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/pdx_amt2_lineage_expression.pdf", width = 16.5, height = 8.3)
ggsave("plots/gene_dynamics/pdx_amt2_lineage_expression.png", width = 16.5, height = 8.3)

lin2.up.heatmap <- plotHeatmap(pdx.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage2.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin2.up.heatmap$gtable$grobs[[6]]$children$GRID.text.19269$label <- "Pseudotime"
lin2.up.heatmap$gtable$grobs[[6]]$children$GRID.text.19273$label <- "Cluster Name"

lin2.up.heatmap

pdf("plots/gene_dynamics/pdx_amt2_up_top50.pdf", width = 8.3, height = 8.3)
lin2.up.heatmap
dev.off()

png("plots/gene_dynamics/pdx_amt2_up_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin2.up.heatmap
dev.off()

## Linear model fitting for lineage 3
lineage3 <- testPseudotime(pdx.slingshot.sce, pseudotime = pdx.slingshot.sce$slingPseudotime_3)
lineage3.genes <- lineage3[order(lineage3$p.value), ]
head(lineage3.genes, 10)

saveRDS(lineage3.genes, "data/pdx_dynamic_genes/full_lists/lineage3_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 3
lineage3.genes.down <- as.data.frame(lineage3.genes[lineage3.genes$logFC < 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage3.genes.down, 10)

saveRDS(lineage3.genes.down, "data/pdx_dynamic_genes/lineage3_down_genes.rds")
write.csv(lineage3.genes.down,"data/pdx_dynamic_genes/lineage3_down_genes.csv", row.names = FALSE)

lineage3.genes.down <- readRDS("data/pdx_dynamic_genes/lineage3_down_genes.rds")

lin3.down.gg <- plotExpression(pdx.slingshot.sce, features = head(lineage3.genes.down$gene, 10), x = "slingPseudotime_3", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 3: downregulated")

lin3.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/pdx_amt3_down_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/pdx_amt3_down_top10.png", width = 11.7, height = 8.3)

on.lineage3 <- !is.na(pdx.slingshot.sce$slingPseudotime_3)
lin3.down.heatmap <- plotHeatmap(pdx.slingshot.sce[ ,on.lineage3], order_columns_by = "slingPseudotime_3", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage3.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin3.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin3.down.heatmap$gtable$grobs[[6]]$children$GRID.text.20722$label <- "Pseudotime"
lin3.down.heatmap$gtable$grobs[[6]]$children$GRID.text.20726$label <- "Cluster Name"

lin3.down.heatmap

pdf("plots/gene_dynamics/pdx_amt3_down_top50.pdf", width = 8.3, height = 8.3)
lin3.down.heatmap
dev.off()

png("plots/gene_dynamics/pdx_amt3_down_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin3.down.heatmap
dev.off()

## Genes that increase in expression with increasing pseudotime along lineage 1
lineage3.genes.up <- as.data.frame(lineage3.genes[lineage3.genes$logFC > 0, ]) %>%
  dplyr::filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage3.genes.up, 10)

saveRDS(lineage3.genes.up, "data/pdx_dynamic_genes/lineage3_up_genes.rds")
write.csv(lineage3.genes.up,"data/pdx_dynamic_genes/lineage3_up_genes.csv", row.names = FALSE)

lineage3.genes.up <- readRDS("data/pdx_dynamic_genes/lineage3_up_genes.rds")

lin3.up.gg <- plotExpression(pdx.slingshot.sce, features = head(lineage3.genes.up$gene, 10), x = "slingPseudotime_3", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Lineage 3: upregulated")

lin3.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/pdx_amt3_up_top10.pdf", width = 11.7, height = 8.3)
ggsave("plots/gene_dynamics/pdx_amt3_up_top10.png", width = 11.7, height = 8.3)

lin3.up.gg + lin3.down.gg + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Lineage 3") &
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/pdx_amt3_lineage_expression.pdf", width = 16.5, height = 8.3)
ggsave("plots/gene_dynamics/pdx_amt3_lineage_expression.png", width = 16.5, height = 8.3)

lin3.up.heatmap <- plotHeatmap(pdx.slingshot.sce[ ,on.lineage3], order_columns_by = "slingPseudotime_3", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage3.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin3.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin3.up.heatmap$gtable$grobs[[6]]$children$GRID.text.25079$label <- "Pseudotime"
lin3.up.heatmap$gtable$grobs[[6]]$children$GRID.text.25083$label <- "Cluster Name"

lin3.up.heatmap

pdf("plots/gene_dynamics/pdx_amt3_up_top50.pdf", width = 8.3, height = 8.3)
lin3.up.heatmap
dev.off()

png("plots/gene_dynamics/pdx_amt3_up_top50.png", width = 8.3, height = 8.3,
    units = "in", res = 200)
lin3.up.heatmap
dev.off()
