## NB Trajectory Analysis

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)
library(scater)
library(TSCAN)
library(slingshot)
library(ggplot2)
library(patchwork)
library(scattermore)
library(ggbeeswarm)
library(hues)

setwd("Z:/echen/Updated NB analysis 2")
setwd("~/my_rds/Updated NB analysis 2")

## Cluster colours
#cluster.cols.0.3 <- iwanthue(length(unique(nb.seurat$seurat_clusters.0.3)))
#saveRDS(cluster.cols.0.3, "cluster_cols_03.rds")
cluster.cols.0.3 <- readRDS("cluster_cols_03.rds")

#cluster.cols.0.6 <- iwanthue(length(unique(nb.seurat$seurat_clusters.0.6)))
#saveRDS(cluster.cols.0.6, "cluster_cols_06.rds")
cluster.cols.0.6 <- readRDS("cluster_cols_06.rds")

identity.cols <- c("Stress 1" = "#CCAA11",
                   "EMT ROS" = "#44AA99",
                   "Developmental" = "#EE3377",
                   "NOR differentiation" = "#CC3311",
                   "Cycling 1" = "#882255",
                   "Cycling 2" = "#BB22EE",
                   "Stress 2" = "#EE7733",
                   "Chromatin remodelling" = "#0077BB",
                   "Collagen metabolism" = "#117733",
                   "NFkB signalling" = "#225555",
                   "Wnt activation" = "#332288")

condition.cols <- c("POT" = "#004488",
                    "Untreated" = "#43ACFF",
                    "Cisplatin(1)_ON" = "#DDAA33",
                    "Cisplatin(2)_ON" = "#DDAA33",
                    "Cisplatin_1weeksOFF" = "#C37F26",
                    "Cisplatin_4weeksOFF" = "#9A6A20",
                    "EZH2i_ON" = "#CC6677",
                    "EZH2i_OFF" = "#80404A",
                    "JQ1_ON" = "#44AA99",
                    "JQ1_OFF" = "#2B6A60")

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

## [ AMT signature plots ] ----

nb.seurat <- readRDS("data/nb_seurat_Wbarcodes_manuscript.RDS")

require(dplyr)
nb.seurat$Cluster_Name <- recode(as.character(nb.seurat$seurat_clusters.0.3),
                                 "0" = "Stress 1",
                                 "1" = "EMT ROS",
                                 "2" = "Developmental",
                                 "3" = "NOR differentiation",
                                 "4" = "Cycling 1",
                                 "5" = "Cycling 2",
                                 "6" = "Stress 2",
                                 "7" = "Chromatin remodelling",
                                 "8" = "Collagen metabolism",
                                 "9" = "NFkB signalling",
                                 "10" = "Wnt activation")

cluster.gg <- DimPlot(nb.seurat,
                      group.by = "Cluster_Name", order = TRUE) +
  umap.theme() + labs(title = "Initial Cluster Annotation") +
  scale_colour_manual(values = identity.cols)

cluster.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/umap_initial_clusters.pdf", width = 5.8, height = 5.8)
ggsave("plots/umap_initial_clusters.png", width = 5.8, height = 5.8)

dir.create("plots/signatures/", recursive = TRUE)

adrn.gg <- FeaturePlot(nb.seurat,
                       features = "ADRN.Sig",
                       order = TRUE,
                       min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.matter(100)) +
  labs(title = "vanGroningen\nADRN signature", colour = "Expression") +
  umap.theme()

mes.gg <- FeaturePlot(nb.seurat,
                      features = "MES.Sig",
                      order = TRUE,
                      min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.haline(100)) +
  labs(title = "vanGroningen\nMES signature", colour = "Expression") +
  umap.theme()

adrn.gg + mes.gg &
  theme(text = element_text(size = 15),
        legend.key.size = unit(1, "cm"))

ggsave("plots/signatures/vG_ADRN_MES_umaps.pdf", width = 5.8, height = 8.7)
ggsave("plots/signatures/vG_ADRN_MES_umaps.png", width = 5.8, height = 8.7)

vln1.gg <- VlnPlot(nb.seurat,
                   features = "MES.Sig", group.by = "seurat_clusters.0.3") +
  scale_fill_manual(values = cluster.cols.0.3) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln1.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln1.gg$layers[[2]] <- NULL

vln2.gg <- VlnPlot(nb.seurat,
                   features = "ADRN.Sig", group.by = "seurat_clusters.0.3") +
  scale_fill_manual(values = cluster.cols.0.3) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln2.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln2.gg$layers[[2]] <- NULL

vln3.gg <- VlnPlot(nb.seurat,
                   features = "AMT.score", group.by = "seurat_clusters.0.3") +
  scale_fill_manual(values = cluster.cols.0.3) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln3.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln3.gg$layers[[2]] <- NULL

(vln3.gg | (vln1.gg / vln2.gg)) + plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/vG_ADRN_MES_violins.pdf", width = 8.7, height = 8.7)
ggsave("plots/signatures/vG_ADRN_MES_violins.png", width = 8.7, height = 8.7)

amt.gg <- DimPlot(nb.seurat,
                  group.by = "AMT.state", order = TRUE) +
  umap.theme() + labs(title = "AMT state") +
  scale_colour_manual(values = c("#990099","lightgrey","#F37735"),
                      labels = c("ADRN", "Intermediate", "MES"))

amt.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/umap_amt_state.pdf", width = 5.8, height = 5.8)
ggsave("plots/umap_amt_state.png", width = 5.8, height = 5.8)

## [ Coexpression ] ----

## Which cells coexpress van Groningen ADRN and MES signatures
coexprs = function(x){
  if (x[["MES.Sig"]] > 0.2 & x[["ADRN.Sig"]] > 0.1) {
    result <- 1
  } else {
    result <- 0
  }
  return(result)
}

metadata.df <- subset(nb.seurat@meta.data, select = c("MES.Sig", "ADRN.Sig"))

coexprs.df <- apply(metadata.df, 1, coexprs)
nb.seurat@meta.data <- cbind(nb.seurat@meta.data, Coexprs.Sig = coexprs.df)

dir.create("plots/coexpression", recursive = TRUE)

coexprs.gg <- DimPlot(nb.seurat,
                      group.by = "Coexprs.Sig", order = TRUE) +
  umap.theme() + labs(title = "ADRN and MES Coexpression\n(MES > 0.2, ADRN > 0.1)") +
  scale_colour_manual(values = c("#4477AA", "#EE6677"), labels = c("FALSE", "TRUE")) + 
  labs(colour = "Coexpression")

coexprs.gg

ggsave("plots/coexpression/umap_coexpress.pdf", width = 5.8, height = 5.8)
ggsave("plots/coexpression/umap_coexpress.png", width = 5.8, height = 5.8)

## [ Dyer signature ] ----

nb.seurat <- readRDS("data/nb_seurat.rds")

require(readxl)
dyer.df <- as.data.frame(read_xlsx(path = "data/Dyer_Sig.xlsx", col_names = FALSE))
colnames(dyer.df) <- c("Gene", "Signature")

dyer.sig <- lapply(unique(dyer.df$Signature), function(x){
  dyer.df[dyer.df$Signature==x, "Gene"]
})

names(dyer.sig) <- unique(dyer.df$Signature)
dyer.sig

nb.seurat <- AddModuleScore(nb.seurat, features = dyer.sig, assay = "SCT", seed = 12345, 
                            name = c("MES.Dyer", "SYMP.Dyer","ADRN.Dyer"))
colnames(nb.seurat@meta.data) <- gsub("\\.Dyer[1-9]$", "\\.Dyer", colnames(nb.seurat@meta.data))

vln1.gg <- VlnPlot(nb.seurat, 
                   features = c("MES.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln1.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln1.gg$layers[[2]] <- NULL

vln2.gg <- VlnPlot(nb.seurat, 
                   features = c("SYMP.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln2.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln2.gg$layers[[2]] <- NULL

vln3.gg <- VlnPlot(nb.seurat, 
                   features = c("ADRN.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln3.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln3.gg$layers[[2]] <- NULL

(vln2.gg | (vln3.gg / vln1.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/dyer_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/dyer_violins.png", width = 8.3, height = 8.3)

vln4.gg <- VlnPlot(nb.seurat, 
                   features = c("MES.Dyer"), group.by = "Condition") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln4.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln4.gg$layers[[2]] <- NULL

vln5.gg <- VlnPlot(nb.seurat, 
                   features = c("SYMP.Dyer"), group.by = "Condition") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln5.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln5.gg$layers[[2]] <- NULL

vln6.gg <- VlnPlot(nb.seurat, 
                   features = c("ADRN.Dyer"), group.by = "Condition") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln6.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln6.gg$layers[[2]] <- NULL

(vln5.gg | (vln6.gg / vln4.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/dyer_condition_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/dyer_condition_violins.png", width = 8.3, height = 8.3)

nb.seurat$Dyer.state <- ifelse(nb.seurat$MES.Dyer > nb.seurat$ADRN.Dyer &
                                 nb.seurat$MES.Dyer > nb.seurat$SYMP.Dyer, "MES",
                               ifelse(nb.seurat$ADRN.Dyer > nb.seurat$MES.Dyer &
                                        nb.seurat$ADRN.Dyer > nb.seurat$SYMP.Dyer, "ADRN",
                                      ifelse(nb.seurat$SYMP.Dyer > nb.seurat$MES.Dyer &
                                               nb.seurat$SYMP.Dyer > nb.seurat$ADRN.Dyer, "SYMP", NA)))

table(nb.seurat$Dyer.state)
## ADRN   MES  SYMP 
## 25557 13980 22858

table(nb.seurat$Condition, nb.seurat$Dyer.state)
##                     ADRN  MES SYMP
## POT                 2234 1412 3999
## Untreated           3287 3030 3153
## Cisplatin(1)_ON      460  730  134
## Cisplatin(2)_ON      135  221   55
## Cisplatin_1weeksOFF 1389  460 1180
## Cisplatin_4weeksOFF 6376 1267 2773
## JQ1_ON              1383 1108  873
## JQ1_OFF             2811 2411 3637
## EZH2i_ON            2540 2269 2096
## EZH2i_OFF           4942 1072 4958

table(nb.seurat$AMT.state)
## ADRN          MES intermediate 
## 27195        22986        12214

table(nb.seurat$Condition, nb.seurat$AMT.state)
##                     ADRN  MES intermediate
## POT                 2415 3768         1462
## Untreated           3780 4472         1218
## Cisplatin(1)_ON       56 1046          222
## Cisplatin(2)_ON       17  361           33
## Cisplatin_1weeksOFF 1126  963          940
## Cisplatin_4weeksOFF 5862 1351         3203
## JQ1_ON               775 1755          834
## JQ1_OFF             2647 5245          967
## EZH2i_ON            2861 2650         1394
## EZH2i_OFF           7656 1375         1941

saveRDS(nb.seurat, "data/nb_seurat_dyer.rds")

nb.seurat <- readRDS("data/nb_seurat_dyer.rds")

dyer.gg <- DimPlot(nb.seurat,
                   group.by = "Dyer.state", order = TRUE) +
  umap.theme() + labs(title = "Dyer state") +
  scale_colour_manual(values = c("#990099", "#F37735", "#009999"),
                      labels = c("ADRN", "MES", "SYMP"))

dyer.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/umap_dyer_state.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/umap_dyer_state.png", width = 5.8, height = 5.8)

dyer.df <- as.data.frame(nb.seurat@meta.data["Dyer.state"])

bar.gg <-
  dyer.df %>%
  ggplot(aes(Dyer.state, fill = Dyer.state)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values = c("#990099", "#F37735", "#009999")) +
  scale_x_discrete(labels = c("ADRN", "MES", "SYMP")) +
  labs(x = "",
       y = "Frequency",
       title = "Dyer states",
       fill = "Annotation") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

bar.gg + theme(text = element_text(size = 15))

ggsave("plots/signatures/dyer_states_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/dyer_states_bar.png", width = 5.8, height = 5.8)

nb.seurat$AMT.state <- factor(nb.seurat$AMT.state, levels = c("ADRN", "intermediate", "MES"))

vln7.gg <- VlnPlot(nb.seurat, 
                   features = c("SYMP.Dyer"), group.by = "AMT.state") +
  scale_fill_manual(values = c("#990099","lightgrey","#F37735"),
                      labels = c("ADRN", "Intermediate", "MES")) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln7.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln7.gg$layers[[2]] <- NULL

vln7.gg &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/vg_symp_violins.pdf", width = 8.3, height = 5.8)
ggsave("plots/signatures/vg_symp_violins.png", width = 8.3, height = 5.8)

## [ Entropy root ] ----

nb.seurat <- readRDS("data/nb_seurat.rds")
nb.sce <- as.SingleCellExperiment(nb.seurat)
reducedDim(nb.sce, "PCA") <- Embeddings(nb.seurat, "pca")
reducedDim(nb.sce, "Harmony") <- Embeddings(nb.seurat, "Harmony")
reducedDim(nb.sce, "UMAP") <- Embeddings(nb.seurat, "umap")
saveRDS(nb.sce, "data/nb_sce.rds")
nb.sce <- readRDS("data/nb_sce.rds")

## Entropy is a measure of gene expression variability
## Higher entropy = greater diversity
## It can be used to root trajectories to the more un-/differentiated states

## Computing the entropy of each cell's expression profile 
entropy <- perCellEntropy(nb.sce)
entropy.df <- data.frame(cluster_0.6 = nb.sce$seurat_clusters.0.6, cluster_0.3 = nb.sce$seurat_clusters.0.3,
                         state = nb.sce$AMT.state, entropy = entropy)

## Cell entropy per cluster
dir.create("plots/entropy/", recursive = TRUE)

entropy.vln1 <- ggplot(entropy.df, aes(x = cluster_0.6, y = entropy, fill = cluster_0.6)) +
  geom_violin(draw_quantiles = 0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cluster.cols.0.6) +
  labs(x = "Cluster number (res = 0.6)", y = "Entropy") +
  guides(fill = guide_legend(ncol = 2, bycol = TRUE))

entropy.vln2 <- ggplot(entropy.df, aes(x = cluster_0.3, y = entropy, fill = cluster_0.3)) +
  geom_violin(draw_quantiles = 0.5) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = cluster.cols.0.3) +
  labs(x = "Cluster number (res = 0.3)", y = "Entropy")

entropy.vln1 + entropy.vln2 &
  theme(text = element_text(size = 15))

ggsave("plots/entropy/cluster_violin.pdf", width = 8.7, height = 5.8)
ggsave("plots/entropy/cluster_violin.png", width = 8.7, height = 5.8)

## Cell entropy by AMT state
entropy.vln3 <- ggplot(entropy.df, aes(x = state, y = entropy, fill = state)) +
  geom_violin(draw_quantiles = 0.5) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("purple","grey","orange")) +
  labs(x = "AMT state", y = "Entropy")

entropy.vln3 +
  theme(text = element_text(size = 15))

ggsave("plots/entropy/state_violin.pdf", width = 5.8, height = 5.8)
ggsave("plots/entropy/state_violin.png", width = 5.8, height = 5.8)

## [ Trajectory inference ] ----

nb.seurat <- readRDS("data/nb_seurat.rds")

dim1.gg <- DimPlot(nb.seurat,
                   group.by = "seurat_clusters.0.6", order = TRUE, label = TRUE, repel = TRUE) +
  scale_colour_manual(values = iwanthue(length(unique(nb.seurat$seurat_clusters.0.6)))) +
  umap.theme() + labs(title = "Harmony Clusters res = 0.6") + theme(legend.position = "none")

dim1.gg$layers[[2]]$aes_params$bg.colour <- "white"
dim1.gg$layers[[2]]$aes_params$direction <- "x"
dim1.gg$layers[[2]]$geom_params$max.overlaps <- 100

dim1.gg

ggsave("plots/umap_harmony_clusters_high_res.pdf", width = 5.8, height = 5.8)
ggsave("plots/umap_harmony_clusters_high_res.png", width = 5.8, height = 5.8)

rm(nb.seurat)

nb.sce <- readRDS("data/nb_sce.rds")

nb.slingshot.sce <- slingshot(nb.sce, reducedDim = "PCA",
                              clusterLabels = "seurat_clusters.0.6", start.clus = 2)
## Takes around 5 hours to complete with 16 cores

saveRDS(nb.slingshot.sce, "data/nb_sce_slingshot.rds")

nb.slingshot.sce <- readRDS("data/nb_sce_slingshot.rds")

## Extract matrix of pseudotime values along each lineage
slingshot.lineages <- slingPseudotime(nb.slingshot.sce)
head(slingshot.lineages)
## 7 lineages

## Single pseudotime for all cells
shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)

dir.create("plots/slingshot/", recursive = TRUE)

require(viridis)
umap.gg <- plotUMAP(nb.slingshot.sce, colour_by = I(shared.pseudotime)) +
  scale_colour_viridis() +
  labs(colour = "Pseudotime")
slingshot.embedded.all <- embedCurves(nb.slingshot.sce, "UMAP")
slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
for (lineage in slingshot.embedded.all) {
  slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
  umap.gg <- umap.gg + geom_path(data = slingshot.embedded.all,
                                 aes(x = umap_1, y = umap_2), size = 1.2)
}

saveRDS(slingshot.embedded.all, "data/slingshot_embedded_all.rds")

umap.gg + theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/umap_slingshot.pdf", width = 8.7, height = 5.8)
ggsave("plots/slingshot/umap_slingshot.png", width = 8.7, height = 5.8)

## Save as SDS
slingshot.embedded.all.sds <- embedCurves(nb.slingshot.sce, "UMAP")
slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)
saveRDS(slingshot.embedded.all.sds, "data/slingshot_sds_embedded_all.rds")

## Plot pseudotime values for cells in each lineage
no.cols <- 3
pseudotime <- slingPseudotime(nb.slingshot.sce)
names <- colnames(pseudotime)
no.rows <- ceiling(length(names)/no.cols)
pal <- viridis(100)
#pdf("plots/slingshot/slingshot_pseudotime_lineages.pdf", width = 8.7, height = 8.7)
png("plots/slingshot/slingshot_pseudotime_lineages.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
par(mfrow = c(no.rows, no.cols))
for (i in names){
  cols <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(nb.slingshot.sce, "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
  lines(slingshot.embedded.all.sds,
        lwd = 1, col = "black", type = "lineages", cex = 1)
}
dev.off()

## Embed lineages one by one
slingshot.embedded <- embedCurves(nb.slingshot.sce, "UMAP")

lineage.names <- colnames(slingshot.embedded)
umap.clust.gg <- plotUMAP(nb.slingshot.sce, colour_by = "Cluster_Name") +
  umap.theme() + scale_colour_manual(values = identity.cols) +
  theme(text = element_text(size = 15)) + guides(colour = guide_legend(override.aes = list(size = 5)))
umap.amt.gg <- plotUMAP(nb.slingshot.sce, colour_by = "AMT.score") +
  umap.theme() + scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))

dir.create("plots/slingshot/by_lineage/", recursive = TRUE)

for (i in 1:length(lineage.names)){
  embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
  embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
  saveRDS(embedded.lineage, paste0("data/slingshot_embedded_", lineage.names[[i]], ".rds"))
  
  plotUMAP(nb.slingshot.sce, colour_by = paste0("slingPseudotime_", i)) +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), size = 1.2) +
    labs(title = lineage.names[[i]]) + umap.theme() +
    theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], ".pdf"), width = 8.7, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], ".png"), width = 8.7, height = 5.8)
  
  umap.clust.gg +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), size = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], "_cluster.pdf"), width = 8.7, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], "_cluster.png"), width = 8.7, height = 5.8)
  
  umap.amt.gg +
    geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), size = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], "_AMT.pdf"), width = 8.7, height = 5.8)
  ggsave(paste0("plots/slingshot/by_lineage/slingshot_", lineage.names[[i]], "_AMT.png"), width = 8.7, height = 5.8)
}

## Lineages of interest for downstream analysis:
## lineage 1 = Collagen-AMT
## Lineage 2 = ROS-AMT_1
## Lineage 4 = NFkB-AMT
## Lineage 5 = ROS-AMT_2
## Lineage 7 = ADRN-ADRN

dir.create("plots/slingshot/lois", recursive = TRUE)

amt1.lineage <- readRDS("data/slingshot_embedded_Lineage1.rds")

amt1.gg <- umap.clust.gg +
  geom_path(data = amt1.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
  labs(title = "Lineage 1:\nDevelopmental ADRN to Collagen MES", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt1.gg

ggsave("plots/slingshot/lois/umap_slingshot_amt1.pdf", width = 8.7, height = 5.8)
ggsave("plots/slingshot/lois/umap_slingshot_amt1.png", width = 8.7, height = 5.8)

amt2.lineage <- readRDS("data/slingshot_embedded_Lineage2.rds")

amt2.gg <- umap.clust.gg +
  geom_path(data = amt2.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
  labs(title = "Lineage 2:\nDevelopmental ADRN to EMT ROS MES", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt2.gg

ggsave("plots/slingshot/lois/umap_slingshot_amt2.pdf", width = 8.7, height = 5.8)
ggsave("plots/slingshot/lois/umap_slingshot_amt2.png", width = 8.7, height = 5.8)

amt4.lineage <- readRDS("data/slingshot_embedded_Lineage4.rds")

amt4.gg <- umap.clust.gg +
  geom_path(data = amt4.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
  labs(title = "Lineage 4:\nDevelopmental ADRN to NFkB MES", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt4.gg

ggsave("plots/slingshot/lois/umap_slingshot_amt4.pdf", width = 8.7, height = 5.8)
ggsave("plots/slingshot/lois/umap_slingshot_amt4.png", width = 8.7, height = 5.8)

amt5.lineage <- readRDS("data/slingshot_embedded_Lineage5.rds")

amt5.gg <- umap.clust.gg +
  geom_path(data = amt5.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
  labs(title = "Lineage 5:\nDevelopmental ADRN to EMT ROS MES", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

amt5.gg

ggsave("plots/slingshot/lois/umap_slingshot_amt5.pdf", width = 8.7, height = 5.8)
ggsave("plots/slingshot/lois/umap_slingshot_amt5.png", width = 8.7, height = 5.8)

aat.lineage <- readRDS("data/slingshot_embedded_Lineage7.rds")

aat.gg <- umap.clust.gg +
  geom_path(data = aat.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
  labs(title = "Lineage 7:\nDevelopmental ADRN to Cycling ADRN", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

aat.gg

ggsave("plots/slingshot/lois/umap_slingshot_aat.pdf", width = 8.7, height = 5.8)
ggsave("plots/slingshot/lois/umap_slingshot_aat.png", width = 8.7, height = 5.8)

## Correlation between AMT score and pseudotime
nb.slingshot.sce <- readRDS("data/nb_sce_slingshot.rds")

slingshot.df <- as.data.frame(colData(nb.slingshot.sce)[, c("slingPseudotime_1", "slingPseudotime_2", "slingPseudotime_4",
                                                            "slingPseudotime_5", "slingPseudotime_7", "AMT.score")])
require(ggpubr)
lin1.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_1, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 100, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "Lineage 1: Collagen-AMT",
       colour = "AMT score")

lin1.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/lois/scatter_amt1.pdf", width = 5.8, height = 5.8)
ggsave("plots/slingshot/lois/scatter_amt1.png", width = 5.8, height = 5.8)

lin2.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_2, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 60, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "Lineage 2: ROS-AMT 1",
       colour = "AMT score")

lin2.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/lois/scatter_amt2.pdf", width = 5.8, height = 5.8)
ggsave("plots/slingshot/lois/scatter_amt2.png", width = 5.8, height = 5.8)

lin4.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_4, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 50, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "Lineage 4: NFkB-AMT",
       colour = "AMT score")

lin4.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/lois/scatter_amt4.pdf", width = 5.8, height = 5.8)
ggsave("plots/slingshot/lois/scatter_amt4.png", width = 5.8, height = 5.8)

lin5.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_5, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 45, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "Lineage 5: ROS-AMT 2",
       colour = "AMT score")

lin5.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/lois/scatter_amt5.pdf", width = 5.8, height = 5.8)
ggsave("plots/slingshot/lois/scatter_amt5.png", width = 5.8, height = 5.8)

lin7.gg <- ggplot(slingshot.df,
                  aes(x = slingPseudotime_7, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 40, label.y = 10) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "Lineage 7: ADRN-ADRN",
       colour = "AMT score")

lin7.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

ggsave("plots/slingshot/lois/scatter_aat.pdf", width = 5.8, height = 5.8)
ggsave("plots/slingshot/lois/scatter_aat.png", width = 5.8, height = 5.8)

## [ Gene dynamics: linear model ] ----

nb.slingshot.sce <- readRDS("data/nb_sce_slingshot.rds")

## Linear model fitting for lineage 1 (Collagen-AMT)
lineage1 <- testPseudotime(nb.slingshot.sce, pseudotime = nb.slingshot.sce$slingPseudotime_1)
lineage1.genes <- lineage1[order(lineage1$p.value), ]
head(lineage1.genes, 10)

dir.create("data/dynamic_genes/full_lists", recursive = TRUE)

saveRDS(lineage1.genes, "data/dynamic_genes/full_lists/lineage1_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 1
require(dplyr)
lineage1.genes.down <- as.data.frame(lineage1.genes[lineage1.genes$logFC < 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage1.genes.down, 10)

#dir.create("data/dynamic_genes/", recursive = TRUE)

saveRDS(lineage1.genes.down, "data/dynamic_genes/lineage1_down_genes.rds")
write.csv(lineage1.genes.down,"data/dynamic_genes/lineage1_down_genes.csv", row.names = FALSE)

lineage1.genes.down <- readRDS("data/dynamic_genes/lineage1_down_genes.rds")

dir.create("plots/gene_dynamics/", recursive = TRUE)

lin1.down.gg <- plotExpression(nb.slingshot.sce, features = head(lineage1.genes.down$gene, 10), x = "slingPseudotime_1", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Collagen-AMT: downregulated")

lin1.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt1_down_top10.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/amt1_down_top10.png", width = 11.6, height = 8.7)

on.lineage1 <- !is.na(nb.slingshot.sce$slingPseudotime_1)
lin1.down.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage1.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

require(pheatmap)
lin1.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin1.down.heatmap$gtable$grobs[[6]]$children$GRID.text.5353$label <- "Pseudotime"
lin1.down.heatmap$gtable$grobs[[6]]$children$GRID.text.5357$label <- "Cluster Name"

lin1.down.heatmap

pdf("plots/gene_dynamics/amt1_down_top50.pdf", width = 8.7, height = 8.7)
lin1.down.heatmap
dev.off()

png("plots/gene_dynamics/amt1_down_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin1.down.heatmap
dev.off()

# Genes that increase in expression with increasing pseudotime along lineage 1
lineage1.genes.up <- as.data.frame(lineage1.genes[lineage1.genes$logFC > 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage1.genes.up, 10)

saveRDS(lineage1.genes.up, "data/dynamic_genes/lineage1_up_genes.rds")
write.csv(lineage1.genes.up,"data/dynamic_genes/lineage1_up_genes.csv", row.names = FALSE)

lineage1.genes.up <- readRDS("data/dynamic_genes/lineage1_up_genes.rds")

lin1.up.gg <- plotExpression(nb.slingshot.sce, features = head(lineage1.genes.up$gene, 10), x = "slingPseudotime_1", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Collagen-AMT: upregulated")

lin1.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt1_up_top10.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/amt1_up_top10.png", width = 11.6, height = 8.7)

lin1.up.gg + lin1.down.gg + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "Collagen-AMT lineage") &
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt1_lineage_expression_manuscript.pdf", width = 18.8, height = 8.7)
ggsave("plots/gene_dynamics/amt1_lineage_expression_manuscript.png", width = 18.8, height = 8.7)

lin1.up.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage1.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin1.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin1.up.heatmap$gtable$grobs[[6]]$children$GRID.text.9818$label <- "Pseudotime"
lin1.up.heatmap$gtable$grobs[[6]]$children$GRID.text.9822$label <- "Cluster Name"

lin1.up.heatmap

pdf("plots/gene_dynamics/amt1_up_top50.pdf", width = 8.7, height = 8.7)
lin1.up.heatmap
dev.off()

png("plots/gene_dynamics/amt1_up_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin1.up.heatmap
dev.off()

## Linear model fitting for lineage 2 (ROS-EMT 1)
lineage2 <- testPseudotime(nb.slingshot.sce, pseudotime = nb.slingshot.sce$slingPseudotime_2)
lineage2.genes <- lineage2[order(lineage2$p.value), ]
head(lineage2.genes, 10)

saveRDS(lineage2.genes, "data/dynamic_genes/full_lists/lineage2_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 2
lineage2.genes.down <- as.data.frame(lineage2.genes[lineage2.genes$logFC < 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage2.genes.down, 10)

saveRDS(lineage2.genes.down, "data/dynamic_genes/lineage2_down_genes.rds")
write.csv(lineage2.genes.down,"data/dynamic_genes/lineage2_down_genes.csv", row.names = FALSE)

lineage2.genes.down <- readRDS("data/dynamic_genes/lineage2_down_genes.rds")

lin2.down.gg <- plotExpression(nb.slingshot.sce, features = head(lineage2.genes.down$gene, 10), x = "slingPseudotime_2", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "ROS-AMT 1: downregulated")

lin2.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt2_down_top10.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/amt2_down_top10.png", width = 11.6, height = 8.7)

on.lineage2 <- !is.na(nb.slingshot.sce$slingPseudotime_2)
lin2.down.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage2.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin2.down.heatmap$gtable$grobs[[6]]$children$GRID.text.15750$label <- "Pseudotime"
lin2.down.heatmap$gtable$grobs[[6]]$children$GRID.text.15754$label <- "Cluster Name"

lin2.down.heatmap

pdf("plots/gene_dynamics/amt2_down_top50.pdf", width = 8.7, height = 8.7)
lin2.down.heatmap
dev.off()

png("plots/gene_dynamics/amt2_down_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin2.down.heatmap
dev.off()

# Genes that increase in expression with increasing pseudotime along lineage 2
lineage2.genes.up <- as.data.frame(lineage2.genes[lineage2.genes$logFC > 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage2.genes.up, 10)

saveRDS(lineage2.genes.up, "data/dynamic_genes/lineage2_up_genes.rds")
write.csv(lineage2.genes.up,"data/dynamic_genes/lineage2_up_genes.csv", row.names = FALSE)

lineage2.genes.up <- readRDS("data/dynamic_genes/lineage2_up_genes.rds")

lin2.up.gg <- plotExpression(nb.slingshot.sce, features = head(lineage2.genes.up$gene, 10), x = "slingPseudotime_2", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "ROS-EMT 1: upregulated")

lin2.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt2_up_top10.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/amt2_up_top10.png", width = 11.6, height = 8.7)

lin2.up.gg + lin2.down.gg + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(title = "ROS-EMT 1 lineage") &
  theme(text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt2_lineage_expression_manuscript.pdf", width = 18.8, height = 8.7)
ggsave("plots/gene_dynamics/amt2_lineage_expression_manuscript.png", width = 18.8, height = 8.7)

lin2.up.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage2.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin2.up.heatmap$gtable$grobs[[6]]$children$GRID.text.20215$label <- "Pseudotime"
lin2.up.heatmap$gtable$grobs[[6]]$children$GRID.text.20219$label <- "Cluster Name"

lin2.up.heatmap

pdf("plots/gene_dynamics/amt2_up_top50.pdf", width = 8.7, height = 8.7)
lin2.up.heatmap
dev.off()

png("plots/gene_dynamics/amt2_up_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin2.up.heatmap
dev.off()

## Linear model fitting for lineage 4 (NFkB-AMT)
lineage4 <- testPseudotime(nb.slingshot.sce, pseudotime = nb.slingshot.sce$slingPseudotime_4)
lineage4.genes <- lineage4[order(lineage4$p.value), ]
head(lineage4.genes, 10)

saveRDS(lineage4.genes, "data/dynamic_genes/full_lists/lineage4_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 4
lineage4.genes.down <- as.data.frame(lineage4.genes[lineage4.genes$logFC < 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage4.genes.down, 10)

saveRDS(lineage4.genes.down, "data/dynamic_genes/lineage4_down_genes.rds")
write.csv(lineage4.genes.down,"data/dynamic_genes/lineage4_down_genes.csv", row.names = FALSE)

lin4.down.gg <- plotExpression(nb.slingshot.sce, features = head(lineage4.genes.down$gene, 10), x = "slingPseudotime_4", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "NFkB-AMT: downregulated")

lin4.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt4_down_top10.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/amt4_down_top10.png", width = 11.6, height = 8.7)

on.lineage4 <- !is.na(nb.slingshot.sce$slingPseudotime_4)
lin4.down.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage4], order_columns_by = "slingPseudotime_4", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage4.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin4.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin4.down.heatmap$gtable$grobs[[6]]$children$GRID.text.21704$label <- "Pseudotime"
lin4.down.heatmap$gtable$grobs[[6]]$children$GRID.text.21708$label <- "Cluster Name"

lin4.down.heatmap

pdf("plots/gene_dynamics/amt4_down_top50.pdf", width = 8.7, height = 8.7)
lin4.down.heatmap
dev.off()

png("plots/gene_dynamics/amt4_down_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin4.down.heatmap
dev.off()

# Genes that increase in expression with increasing pseudotime along lineage 4
lineage4.genes.up <- as.data.frame(lineage4.genes[lineage4.genes$logFC > 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage4.genes.up, 10)

saveRDS(lineage4.genes.up, "data/dynamic_genes/lineage4_up_genes.rds")
write.csv(lineage4.genes.up,"data/dynamic_genes/lineage4_up_genes.csv", row.names = FALSE)

lin4.up.gg <- plotExpression(nb.slingshot.sce, features = head(lineage4.genes.up$gene, 10), x = "slingPseudotime_4", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "NFkB-AMT: upregulated")

lin4.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt4_up_top10.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/amt4_up_top10.png", width = 11.6, height = 8.7)

lin4.up.gg + lin4.down.gg + plot_layout(guides = "collect") +
  plot_annotation(title = "NFkB-AMT lineage") &
  theme(text = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt4_lineage_expression.pdf", width = 18.8, height = 8.7)
ggsave("plots/gene_dynamics/amt4_lineage_expression.png", width = 18.8, height = 8.7)

lin4.up.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage4], order_columns_by = "slingPseudotime_4", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage4.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin4.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin4.up.heatmap$gtable$grobs[[6]]$children$GRID.text.30612$label <- "Pseudotime"
lin4.up.heatmap$gtable$grobs[[6]]$children$GRID.text.30616$label <- "Cluster Name"

lin4.up.heatmap

pdf("plots/gene_dynamics/amt4_up_top50.pdf", width = 8.7, height = 8.7)
lin4.up.heatmap
dev.off()

png("plots/gene_dynamics/amt4_up_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin4.up.heatmap
dev.off()

## Linear model fitting for lineage 5 (ROS-EMT 2)
lineage5 <- testPseudotime(nb.slingshot.sce, pseudotime = nb.slingshot.sce$slingPseudotime_5)
lineage5.genes <- lineage5[order(lineage5$p.value), ]
head(lineage5.genes, 10)

saveRDS(lineage5.genes, "data/dynamic_genes/full_lists/lineage5_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 5
lineage5.genes.down <- as.data.frame(lineage5.genes[lineage5.genes$logFC < 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage5.genes.down, 10)

saveRDS(lineage5.genes.down, "data/dynamic_genes/lineage5_down_genes.rds")
write.csv(lineage5.genes.down,"data/dynamic_genes/lineage5_down_genes.csv", row.names = FALSE)

lin5.down.gg <- plotExpression(nb.slingshot.sce, features = head(lineage5.genes.down$gene, 10), x = "slingPseudotime_5", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "ROS-AMT 2: downregulated")

lin5.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt5_down_top10.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/amt5_down_top10.png", width = 11.6, height = 8.7)

on.lineage5 <- !is.na(nb.slingshot.sce$slingPseudotime_5)
lin5.down.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage5], order_columns_by = "slingPseudotime_5", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage5.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin5.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin5.down.heatmap$gtable$grobs[[6]]$children$GRID.text.35077$label <- "Pseudotime"
lin5.down.heatmap$gtable$grobs[[6]]$children$GRID.text.35081$label <- "Cluster Name"

lin5.down.heatmap

pdf("plots/gene_dynamics/amt5_down_top50.pdf", width = 8.7, height = 8.7)
lin5.down.heatmap
dev.off()

png("plots/gene_dynamics/amt5_down_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin5.down.heatmap
dev.off()

# Genes that increase in expression with increasing pseudotime along lineage 5
lineage5.genes.up <- as.data.frame(lineage5.genes[lineage5.genes$logFC > 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage5.genes.up, 10)

saveRDS(lineage5.genes.up, "data/dynamic_genes/lineage5_up_genes.rds")
write.csv(lineage5.genes.up,"data/dynamic_genes/lineage5_up_genes.csv", row.names = FALSE)

lin5.up.gg <- plotExpression(nb.slingshot.sce, features = head(lineage5.genes.up$gene, 10), x = "slingPseudotime_5", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "ROS-EMT 2: upregulated")

lin5.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt5_up_top10.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/amt5_up_top10.png", width = 11.6, height = 8.7)

lin5.up.gg + lin5.down.gg + plot_layout(guides = "collect") +
  plot_annotation(title = "ROS-EMT 2 lineage") &
  theme(text = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/amt5_lineage_expression.pdf", width = 18.8, height = 8.7)
ggsave("plots/gene_dynamics/amt5_lineage_expression.png", width = 18.8, height = 8.7)

lin5.up.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage5], order_columns_by = "slingPseudotime_5", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage5.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin5.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin5.up.heatmap$gtable$grobs[[6]]$children$GRID.text.39542$label <- "Pseudotime"
lin5.up.heatmap$gtable$grobs[[6]]$children$GRID.text.39546$label <- "Cluster Name"

lin5.up.heatmap

pdf("plots/gene_dynamics/amt5_up_top50.pdf", width = 8.7, height = 8.7)
lin5.up.heatmap
dev.off()

png("plots/gene_dynamics/amt5_up_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin5.up.heatmap
dev.off()

## Linear model fitting for lineage 7 (ADRN-ADRN)
lineage7 <- testPseudotime(nb.slingshot.sce, pseudotime = nb.slingshot.sce$slingPseudotime_7)
lineage7.genes <- lineage7[order(lineage7$p.value), ]
head(lineage7.genes, 10)

saveRDS(lineage7.genes, "data/dynamic_genes/full_lists/lineage7_genes_all.rds")

## Genes that decrease in expression with increasing pseudotime along lineage 7
lineage7.genes.down <- as.data.frame(lineage7.genes[lineage7.genes$logFC < 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage7.genes.down, 10)

saveRDS(lineage7.genes.down, "data/dynamic_genes/lineage7_down_genes.rds")
write.csv(lineage7.genes.down,"data/dynamic_genes/lineage7_down_genes.csv", row.names = FALSE)

lin7.down.gg <- plotExpression(nb.slingshot.sce, features = head(lineage7.genes.down$gene, 10), x = "slingPseudotime_7", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "ADRN-ADRN: downregulated")

lin7.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/aat_down_top10.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/aat_down_top10.png", width = 11.6, height = 8.7)

on.lineage7 <- !is.na(nb.slingshot.sce$slingPseudotime_7)
lin7.down.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage7], order_columns_by = "slingPseudotime_7", colour_columns_by = "Cluster_Name",
                                 column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage7.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin7.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin7.down.heatmap$gtable$grobs[[6]]$children$GRID.text.41031$label <- "Pseudotime"
lin7.down.heatmap$gtable$grobs[[6]]$children$GRID.text.41035$label <- "Cluster Name"

lin7.down.heatmap

pdf("plots/gene_dynamics/aat_down_top50.pdf", width = 8.7, height = 8.7)
lin7.down.heatmap
dev.off()

png("plots/gene_dynamics/aat_down_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin7.down.heatmap
dev.off()

# Genes that increase in expression with increasing pseudotime along lineage 7
lineage7.genes.up <- as.data.frame(lineage7.genes[lineage7.genes$logFC > 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage7.genes.up, 10)

saveRDS(lineage7.genes.up, "data/dynamic_genes/lineage7_up_genes.rds")
write.csv(lineage7.genes.up,"data/dynamic_genes/lineage7_up_genes.csv", row.names = FALSE)

lin7.up.gg <- plotExpression(nb.slingshot.sce, features = head(lineage7.genes.up$gene, 10), x = "slingPseudotime_7", colour_by = "Cluster_Name") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "ADRN-ADRN: upregulated")

lin7.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/aat_up_top7.pdf", width = 11.6, height = 8.7)
ggsave("plots/gene_dynamics/aat_up_top7.png", width = 11.6, height = 8.7)

lin7.up.gg + lin7.down.gg + plot_layout(guides = "collect") +
  plot_annotation(title = "ADRN-ADRN lineage") &
  theme(text = element_text(size = 20)) &
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

ggsave("plots/gene_dynamics/aat_lineage_expression.pdf", width = 18.8, height = 8.7)
ggsave("plots/gene_dynamics/aat_lineage_expression.png", width = 18.8, height = 8.7)

lin7.up.heatmap <- plotHeatmap(nb.slingshot.sce[ ,on.lineage7], order_columns_by = "slingPseudotime_7", colour_columns_by = "Cluster_Name",
                               column_annotation_colours = list(Cluster_Name = identity.cols), features = head(lineage7.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin7.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster Name", "Pseudotime")

lin7.up.heatmap$gtable$grobs[[6]]$children$GRID.text.45496$label <- "Pseudotime"
lin7.up.heatmap$gtable$grobs[[6]]$children$GRID.text.45500$label <- "Cluster Name"

lin7.up.heatmap

pdf("plots/gene_dynamics/aat_up_top50.pdf", width = 8.7, height = 8.7)
lin7.up.heatmap
dev.off()

png("plots/gene_dynamics/aat_up_top50.png", width = 8.7, height = 8.7,
    units = "in", res = 200)
lin7.up.heatmap
dev.off()
