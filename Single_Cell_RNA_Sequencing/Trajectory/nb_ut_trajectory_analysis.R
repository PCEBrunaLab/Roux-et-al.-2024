## NB trajectory analysis
## Untreated samples

## Feb 2025

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)
library(SeuratDisk)
library(scater)
library(scran)
library(TSCAN)
library(slingshot)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(scattermore)
library(reshape2)
library(patchwork)
library(pheatmap)
library(hues)
library(viridis)

setwd('/Users/echen/Library/CloudStorage/OneDrive-TheInstituteofCancerResearch/Documents/nb ut analysis 2025')

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

## [ Trajectory inference ] ----

## Prep data
sknsh.seurat <- readRDS("data/sknsh_pot_seurat.RDS")
shep.seurat <- readRDS("data/shep_pot_seurat.RDS")
shsy5y.seurat <- readRDS("data/shsy5y_pot_seurat.RDS")
nb039.seurat <- readRDS("data/nb039_pot_seurat.RDS")
nb067.seurat <- readRDS("data/nb067_pot_seurat.RDS")
grnb5.seurat <- readRDS("data/grnb5_pot_seurat.RDS")

ut.samples <- list(sknsh.seurat, shep.seurat, shsy5y.seurat,
                   nb039.seurat, nb067.seurat, grnb5.seurat)
names(ut.samples) <- c("SK-N-SH", "SH-EP", "SH-SY5Y",
                       "NB039", "NB067", "GRNB5")

for (object in 1:length(ut.samples)) {
  umap.gg <- DimPlot(ut.samples[[object]],
                      group.by = "seurat_clusters.0.4", order = TRUE, label = TRUE, repel = TRUE) +
    scale_colour_manual(values = iwanthue(length(unique(ut.samples[[object]]$seurat_clusters.0.4)))) +
    umap.theme() + labs(title = paste0(names(ut.samples)[object], ": Harmony Clusters (res = 0.4)")) +
    theme(legend.position = "none")
  umap.gg$layers[[2]]$aes_params$bg.colour <- "white"
  umap.gg$layers[[2]]$aes_params$direction <- "x"
  umap.gg$layers[[2]]$geom_params$max.overlaps <- 100
  umap.gg
  ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.samples)[object])), "_pot_umap_harmony_clusters.pdf"), width = 5.8, height = 5.8)
  ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.samples)[object])), "_pot_umap_harmony_clusters.png"), width = 5.8, height = 5.8)
  
  sce <- as.SingleCellExperiment(ut.samples[[object]])
  reducedDim(sce, "PCA") <- Embeddings(ut.samples[[object]], "pca")
  reducedDim(sce, "Harmony") <- Embeddings(ut.samples[[object]], "Harmony")
  reducedDim(sce, "UMAP") <- Embeddings(ut.samples[[object]], "umap")
  saveRDS(sce, paste0("data/", gsub("-", "", tolower(names(ut.samples)[object])), "pot_sce.rds"))
}

## Run Slingshot
sknsh.sce <- readRDS("data/sknsh_pot_sce.rds")
sknsh.sling2 <- slingshot(sknsh.sce, reducedDim = "PCA",
                          clusterLabels = "seurat_clusters.0.4", start.clus = 2)
saveRDS(sknsh.sling2, "data/sknsh_pot_sce_sling2.rds")

shep.sce <- readRDS("data/shep_pot_sce.rds")
shep.sling5 <- slingshot(shep.sce, reducedDim = "PCA",
                         clusterLabels = "seurat_clusters.0.4", start.clus = 5)
saveRDS(shep.sling5, "data/shep_pot_sce_sling5.rds")

shsy5y.sce <- readRDS("data/shsy5y_pot_sce.rds")
shsy5y.sling4 <- slingshot(shsy5y.sce, reducedDim = "PCA",
                         clusterLabels = "seurat_clusters.0.4", start.clus = 4)
saveRDS(shsy5y.sling4, "data/shsy5y_pot_sce_sling4.rds")

nb039.sce <- readRDS("data/nb039_pot_sce.rds")
nb039.sling3 <- slingshot(nb039.sce, reducedDim = "PCA",
                          clusterLabels = "seurat_clusters.0.4", start.clus = 3)
saveRDS(nb039.sling3, "data/nb039_pot_sce_sling3.rds")

nb067.sce <- readRDS("data/nb067_pot_sce.rds")
nb067.sling3 <- slingshot(nb067.sce, reducedDim = "PCA",
                          clusterLabels = "seurat_clusters.0.4", start.clus = 3)
saveRDS(nb067.sling3, "data/nb067_pot_sce_sling3.rds")

grnb5.sce <- readRDS("data/grnb5_pot_sce.rds")
grnb5.sling3 <- slingshot(grnb5.sce, reducedDim = "PCA",
                          clusterLabels = "seurat_clusters.0.4", start.clus = 3)
saveRDS(grnb5.sling3, "data/grnb5_pot_sce_sling3.rds")

sknsh.sling2 <- readRDS("data/sknsh_pot_sce_sling2.rds")
shep.sling5 <- readRDS("data/shep_pot_sce_sling5.rds")
shsy5y.sling4 <- readRDS("data/shsy5y_pot_sce_sling4.rds")
nb039.sling3 <- readRDS("data/nb039_pot_sce_sling3.rds")
nb067.sling3 <- readRDS("data/nb067_pot_sce_sling3.rds")
grnb5.sling3 <- readRDS("data/grnb5_pot_sce_sling3.rds")

ut.sling <- list(sknsh.sling2, shep.sling5, shsy5y.sling4,
                 nb039.sling3, nb067.sling3, grnb5.sling3)
names(ut.sling) <- c("SK-N-SH", "SH-EP", "SH-SY5Y",
                       "NB039", "NB067", "GRNB5")

for (object in 1:length(ut.sling)) {
  slingshot.lineages <- slingPseudotime(ut.sling[[object]])
  head(slingshot.lineages)
  shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)
  
  sling.gg <- plotUMAP(ut.sling[[object]], colour_by = I(shared.pseudotime)) +
    scale_colour_viridis() +
    labs(colour = "Pseudotime")
  slingshot.embedded.all <- embedCurves(ut.sling[[object]], "UMAP")
  slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
  for (lineage in slingshot.embedded.all) {
    slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
    sling.gg <- sling.gg + geom_path(data = slingshot.embedded.all,
                                     aes(x = umap_1, y = umap_2), linewidth = 1.2)
  }
  saveRDS(slingshot.embedded.all, paste0("data/", gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_embedded_all.rds"))
  
  sling.gg + labs(title = paste0(names(ut.sling)[object], " Slingshot")) +
    theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))
  ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_umap_slingshot.pdf"), width = 8.7, height = 5.8)
  ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_umap_slingshot.png"), width = 8.7, height = 5.8)
  
  slingshot.embedded.all.sds <- embedCurves(ut.sling[[object]], "UMAP")
  slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)
  saveRDS(slingshot.embedded.all.sds, paste0("data/",  gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_sds_embedded_all.rds"))
  
  no.cols <- ceiling(ncol(slingshot.lineages)/3)
  pseudotime <- slingPseudotime(ut.sling[[object]])
  names <- colnames(pseudotime)
  no.rows <- ceiling(length(names)/no.cols)
  pal <- viridis(100)
  pdf(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_pseudotime_lineages.pdf"),
      width = 8.7, height = 8.7)
  par(mfrow = c(no.rows, no.cols))
  for (i in names){
    cols <- pal[cut(pseudotime[,i], breaks = 100)]
    plot(reducedDim(ut.sling[[object]], "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
    lines(slingshot.embedded.all.sds,
          lwd = 1, col = "black", type = "lineages", cex = 1)
  }
  dev.off()
  png(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_pseudotime_lineages.png"),
      width = 8.7, height = 8.7, units = "in", res = 200)
  par(mfrow = c(no.rows, no.cols))
  for (i in names){
    cols <- pal[cut(pseudotime[,i], breaks = 100)]
    plot(reducedDim(ut.sling[[object]], "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
    lines(slingshot.embedded.all.sds,
          lwd = 1, col = "black", type = "lineages", cex = 1)
  }
  dev.off()
  
  slingshot.embedded <- embedCurves(ut.sling[[object]], "UMAP")
  
  lineage.names <- colnames(slingshot.embedded)
  umap.clust.gg <- plotUMAP(ut.sling[[object]], colour_by = "seurat_clusters.0.4") +
    umap.theme() + scale_colour_manual(values = iwanthue(length(unique(ut.sling[[object]]$seurat_clusters.0.4)))) +
    theme(text = element_text(size = 15)) + guides(colour = guide_legend(override.aes = list(size = 5)))
  umap.amt.gg <- plotUMAP(ut.sling[[object]], colour_by = "AMT.score") +
    umap.theme() + scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
    theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
  
  for (i in 1:length(lineage.names)){
    embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
    embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
    saveRDS(embedded.lineage, paste0("data/slingshot_embedded_", lineage.names[[i]], ".rds"))
    
    plotUMAP(ut.sling[[object]], colour_by = paste0("slingPseudotime_", i)) +
      geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
      labs(title = paste0(names(ut.sling)[object], lineage.names[[i]])) + umap.theme() +
      theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
    ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_", lineage.names[[i]], ".pdf"), width = 8.7, height = 5.8)
    ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_", lineage.names[[i]], ".png"), width = 8.7, height = 5.8)
    
    umap.clust.gg +
      geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
      labs(title = paste0(names(ut.sling)[object], lineage.names[[i]]))
    ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_", lineage.names[[i]], "_cluster.pdf"), width = 8.7, height = 5.8)
    ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_", lineage.names[[i]], "_cluster.png"), width = 8.7, height = 5.8)
    
    umap.amt.gg +
      geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
      labs(title = paste0(names(ut.sling)[object], lineage.names[[i]]))
    ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_", lineage.names[[i]], "_amt.pdf"), width = 8.7, height = 5.8)
    ggsave(paste0("plots/", gsub("-", "", tolower(names(ut.sling)[object])), "_slingshot_", lineage.names[[i]], "_amt.png"), width = 8.7, height = 5.8)
  }
}

df.list <- list(as.data.frame(colData(ut.sling[[1]])[, c("slingPseudotime_1", "slingPseudotime_2", "slingPseudotime_3", "AMT.score")]),
                as.data.frame(colData(ut.sling[[2]])[, c("slingPseudotime_1", "slingPseudotime_2", "slingPseudotime_3", "AMT.score")]),
                as.data.frame(colData(ut.sling[[3]])[, c("slingPseudotime_1", "AMT.score")]),
                as.data.frame(colData(ut.sling[[4]])[, c("slingPseudotime_1", "slingPseudotime_2", "AMT.score")]),
                as.data.frame(colData(ut.sling[[5]])[, c("slingPseudotime_1", "slingPseudotime_2", "AMT.score")]),
                as.data.frame(colData(ut.sling[[6]])[, c("slingPseudotime_1", "slingPseudotime_2", "AMT.score")]))
names(df.list) <- names(ut.sling)

for (df in 1:length(df.list)) {
  df.data <- df.list[[df]]
  df.name <- names(df.list[df])
  
  pseudotime.cols <- grep("^slingPseudotime_", colnames(df.data), value = TRUE)
  
  for (pseudotime in pseudotime.cols) {
    ggplot(df.data, aes_string(x = pseudotime, y = "AMT.score", colour = "AMT.score")) +
      geom_scattermore(pointsize = 2) +
      geom_smooth(method = "lm", colour = "black") +
      stat_cor(method = "pearson", label.x = 100, label.y = -15) +
      theme_bw() +
      scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
      labs(x = "Pseudotime",
           y = "AMT score",
           title = paste0(df.name, ": ", pseudotime, " vs AMT score"),
           colour = "AMT score") +
      theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

    ggsave(paste0("plots/", gsub("-", "", tolower(df.name)),"_scatter_", gsub("_", "", tolower(pseudotime)), ".pdf"), width = 8.7, height = 5.8)
    ggsave(paste0("plots/", gsub("-", "", tolower(df.name)),"_scatter_", gsub("_", "", tolower(pseudotime)), ".png"), width = 8.7, height = 5.8)
  }
}

## Adjust text on graph for correlation plots
sknsh.sling2 <- readRDS("data/sknsh_pot_sce_sling2.rds")
shep.sling5 <- readRDS("data/shep_pot_sce_sling5.rds")
shsy5y.sling4 <- readRDS("data/shsy5y_pot_sce_sling4.rds")
nb039.sling3 <- readRDS("data/nb039_pot_sce_sling3.rds")
nb067.sling3 <- readRDS("data/nb067_pot_sce_sling3.rds")
grnb5.sling3 <- readRDS("data/grnb5_pot_sce_sling3.rds")

ut.sling <- list(sknsh.sling2, shep.sling5, shsy5y.sling4,
                 nb039.sling3, nb067.sling3, grnb5.sling3)
names(ut.sling) <- c("SK-N-SH", "SH-EP", "SH-SY5Y",
                     "NB039", "NB067", "GRNB5")

df.list <- list(as.data.frame(colData(ut.sling[[1]])[, c("slingPseudotime_3", "AMT.score")]),
                as.data.frame(colData(ut.sling[[2]])[, c("slingPseudotime_1", "AMT.score")]),
                as.data.frame(colData(ut.sling[[3]])[, c("slingPseudotime_1", "AMT.score")]),
                as.data.frame(colData(ut.sling[[4]])[, c("slingPseudotime_1", "AMT.score")]),
                as.data.frame(colData(ut.sling[[5]])[, c("slingPseudotime_1", "AMT.score")]),
                as.data.frame(colData(ut.sling[[6]])[, c("slingPseudotime_2", "AMT.score")]))
names(df.list) <- names(ut.sling)

for (df in 1:length(df.list)) {
  df.data <- df.list[[df]]
  df.name <- names(df.list[df])
  
  pseudotime.cols <- grep("^slingPseudotime_", colnames(df.data), value = TRUE)
  
  for (pseudotime in pseudotime.cols) {
    x.pos <- max(df.data[pseudotime], na.rm = TRUE)*0.6
    y.pos <- max(df.data["AMT.score"], na.rm = TRUE)*0.2
    
    ggplot(df.data, aes_string(x = pseudotime, y = "AMT.score", colour = "AMT.score")) +
      geom_scattermore(pointsize = 2) +
      geom_smooth(method = "lm", colour = "black") +
      stat_cor(method = "pearson", label.x = x.pos, label.y = -y.pos) +
      theme_bw() +
      scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
      labs(x = "Pseudotime",
           y = "AMT score",
           title = paste0(df.name, ": ", pseudotime, " vs AMT score"),
           colour = "AMT score") +
      theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))
    
    ggsave(paste0("plots/", gsub("-", "", tolower(df.name)),"_scatter_", gsub("_", "", tolower(pseudotime)), "_format.pdf"), width = 8.7, height = 5.8)
    ggsave(paste0("plots/", gsub("-", "", tolower(df.name)),"_scatter_", gsub("_", "", tolower(pseudotime)), "_format.png"), width = 8.7, height = 5.8)
  }
}

## [ Dynamic genes ] ----

## /////////////////////////////////////////////////////////////////////////////
## SK-N-SH /////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

sknsh.sling2 <- readRDS("data/sknsh_pot_sce_sling2.rds")
sknsh.cols <- grep("^slingPseudotime_", colnames(colData(sknsh.sling2)), value = TRUE)
cluster.cols <- iwanthue(length(unique(sknsh.sling2$seurat_clusters.0.4)))
names(cluster.cols) <- levels(sknsh.sling2$seurat_clusters.0.4)

for (pseudotime in sknsh.cols) {
  lineage <- testPseudotime(sknsh.sling2, pseudotime = sknsh.sling2[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("data/trajectory_genes/sknsh_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("data/trajectory_genes/sknsh_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(sknsh.sling2, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(sknsh.sling2[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(sknsh.sling2[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("data/trajectory_genes/sknsh_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("data/trajectory_genes/sknsh_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(sknsh.sling2, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(sknsh.sling2[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("plots/sknsh_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## SH-EP ///////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

shep.sling5 <- readRDS("data/shep_pot_sce_sling5.rds")
shep.cols <- grep("^slingPseudotime_", colnames(colData(shep.sling5)), value = TRUE)
cluster.cols <- iwanthue(length(unique(shep.sling5$seurat_clusters.0.4)))
names(cluster.cols) <- levels(shep.sling5$seurat_clusters.0.4)

for (pseudotime in shep.cols) {
  lineage <- testPseudotime(shep.sling5, pseudotime = shep.sling5[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("data/trajectory_genes/shep_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("data/trajectory_genes/shep_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(shep.sling5, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(shep.sling5[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(shep.sling5[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("data/trajectory_genes/shep_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("data/trajectory_genes/shep_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(shep.sling5, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(shep.sling5[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("plots/shep_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## SH-SY5Y /////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

shsy5y.sling4 <- readRDS("data/shsy5y_pot_sce_sling4.rds")
shsy5y.cols <- grep("^slingPseudotime_", colnames(colData(shsy5y.sling4)), value = TRUE)
cluster.cols <- iwanthue(length(unique(shsy5y.sling4$seurat_clusters.0.4)))
names(cluster.cols) <- levels(shsy5y.sling4$seurat_clusters.0.4)

for (pseudotime in shsy5y.cols) {
  lineage <- testPseudotime(shsy5y.sling4, pseudotime = shsy5y.sling4[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("data/trajectory_genes/shsy5y_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("data/trajectory_genes/shsy5y_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(shsy5y.sling4, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(shsy5y.sling4[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(shsy5y.sling4[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("data/trajectory_genes/shsy5y_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("data/trajectory_genes/shsy5y_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(shsy5y.sling4, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(shsy5y.sling4[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("plots/shsy5y_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## NB039 ///////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

nb039.sling3 <- readRDS("data/nb039_pot_sce_sling3.rds")
nb039.cols <- grep("^slingPseudotime_", colnames(colData(nb039.sling3)), value = TRUE)
cluster.cols <- iwanthue(length(unique(nb039.sling3$seurat_clusters.0.4)))
names(cluster.cols) <- levels(nb039.sling3$seurat_clusters.0.4)

for (pseudotime in nb039.cols) {
  lineage <- testPseudotime(nb039.sling3, pseudotime = nb039.sling3[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("data/trajectory_genes/nb039_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("data/trajectory_genes/nb039_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(nb039.sling3, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(nb039.sling3[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(nb039.sling3[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("data/trajectory_genes/nb039_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("data/trajectory_genes/nb039_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(nb039.sling3, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(nb039.sling3[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("plots/nb039_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## NB067 ///////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

nb067.sling3 <- readRDS("data/nb067_pot_sce_sling3.rds")
nb067.cols <- grep("^slingPseudotime_", colnames(colData(nb067.sling3)), value = TRUE)
cluster.cols <- iwanthue(length(unique(nb067.sling3$seurat_clusters.0.4)))
names(cluster.cols) <- levels(nb067.sling3$seurat_clusters.0.4)

for (pseudotime in nb067.cols) {
  lineage <- testPseudotime(nb067.sling3, pseudotime = nb067.sling3[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("data/trajectory_genes/nb067_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("data/trajectory_genes/nb067_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(nb067.sling3, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(nb067.sling3[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(nb067.sling3[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("data/trajectory_genes/nb067_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("data/trajectory_genes/nb067_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(nb067.sling3, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(nb067.sling3[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("plots/nb067_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## GRNB5 ///////////////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

grnb5.sling3 <- readRDS("data/grnb5_pot_sce_sling3.rds")
grnb5.cols <- grep("^slingPseudotime_", colnames(colData(grnb5.sling3)), value = TRUE)
grnb5.cols <- grnb5.cols[1:2]
cluster.cols <- iwanthue(length(unique(grnb5.sling3$seurat_clusters.0.4)))
names(cluster.cols) <- levels(grnb5.sling3$seurat_clusters.0.4)

for (pseudotime in grnb5.cols) {
  lineage <- testPseudotime(grnb5.sling3, pseudotime = grnb5.sling3[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("data/trajectory_genes/grnb5_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("data/trajectory_genes/grnb5_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(grnb5.sling3, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(grnb5.sling3[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(grnb5.sling3[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("data/trajectory_genes/grnb5_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("data/trajectory_genes/grnb5_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(grnb5.sling3, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(grnb5.sling3[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("plots/grnb5_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## [ Gene list similarity ] ----

dir.path <- "data/trajectory_genes"
rds.files <- list.files(dir.path, pattern = "\\.rds$", full.names = TRUE)

up.list <- list()
down.list <- list()

for (file in rds.files) {
  file.name <- basename(file)
  prefix <- toupper(sub("_.*", "", file.name))
  pseudotime.match <- regmatches(file.name, regexpr("slingpseudotime[0-9]+", file.name))
  pseudotime.num <- sub("slingpseudotime", "", pseudotime.match)
  name <- paste0(prefix, "_", pseudotime.num)

  data <- readRDS(file)
  
  if (grepl("_up_", file.name)) {
    up.list[[name]] <- data[["gene"]]
  } else if (grepl("_down_", file.name)) {
    down.list[[name]] <- data[["gene"]]
  }
}

up.list <- up.list[c("SKNSH_3", "SHEP_1", "SKSY5Y_1", "NB039_1", "NB067_1", "GRNB5_2")]
down.list <- down.list[c("SKNSH_3", "SHEP_1", "SKSY5Y_1", "NB039_1", "NB067_1", "GRNB5_2")]

up.mat <- crossprod(table(stack(up.list)))
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
ggsave("plots/genes_comparison_up_heatmap_examples.pdf", width = 5.8, height = 5.8)
ggsave("plots/genes_comparison_up_heatmap_examples.png", width = 5.8, height = 5.8)

down.mat <- crossprod(table(stack(down.list)))
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
ggsave("plots/genes_comparison_down_heatmap_examples.pdf", width = 5.8, height = 5.8)
ggsave("plots/genes_comparison_down_heatmap_examples.png", width = 5.8, height = 5.8)

## [ Seurat to AnnData ] ----

sknsh.seurat <- readRDS("data/sknsh_pot_seurat.RDS")
shep.seurat <- readRDS("data/shep_pot_seurat.RDS")
shsy5y.seurat <- readRDS("data/shsy5y_pot_seurat.RDS")
nb039.seurat <- readRDS("data/nb039_pot_seurat.RDS")
nb067.seurat <- readRDS("data/nb067_pot_seurat.RDS")
grnb5.seurat <- readRDS("data/grnb5_pot_seurat.RDS")

ut.samples <- list(sknsh.seurat, shep.seurat, shsy5y.seurat,
                   nb039.seurat, nb067.seurat, grnb5.seurat)
names(ut.samples) <- c("SK-N-SH", "SH-EP", "SH-SY5Y",
                       "NB039", "NB067", "GRNB5")

for (object in 1:length(ut.samples)) {
  DefaultAssay(ut.samples[[object]]) <- "RNA"
  ut.samples[[object]] <- NormalizeData(ut.samples[[object]])
  ut.samples[[object]] <- FindVariableFeatures(ut.samples[[object]], nfeatures = 1000)
  ut.samples[[object]] <- ScaleData(ut.samples[[object]])
  ut.samples[[object]][["RNA"]] <- as(ut.samples[[object]][["RNA"]], "Assay")
  ut.samples[[object]][["RNA"]]$data <- NULL
  ut.samples[[object]][["RNA"]]$scale.data <- NULL
  
  SaveH5Seurat(ut.samples[[object]], filename = paste0("data/", gsub("-", "", tolower(names(ut.samples)[object])),"_seurat.h5Seurat"), overwrite = T, verbose = T)
  Convert(paste0("data/", gsub("-", "", tolower(names(ut.samples)[object])),"_seurat.h5Seurat"), dest = "h5ad", assay = "RNA", overwrite = T)
}

# # Only works on Assay v4 objects
# # This issue is because the object returned from the [[.Assay function in the SeuratObject package was changed in v5
# # Until an official update of the SeuratDisk package, this is a workaround (not tested, taken from https://github.com/mojaveazure/seurat-disk/issues/147)
# # Assigning the previous version of the `[[` function for the Assay class to the SeuratDisk package environment
# 
# "[[.Assay" <- function(x, i, ..., drop = FALSE) {
#   if (missing(x = i)) {
#     i <- colnames(x = slot(object = x, name = 'meta.features'))
#   }
#   data.return <- slot(object = x, name = 'meta.features')[, i, drop = FALSE, ...]
#   if (drop) {
#     data.return <- unlist(x = data.return, use.names = FALSE)
#     names(x = data.return) <- rep.int(x = rownames(x = x), times = length(x = i))
#   }
#   return(data.return)
# }
# environment(`[[.Assay`) <- asNamespace("SeuratObject")
# rlang::env_unlock(asNamespace("SeuratDisk"))
# assign("[[.Assay", `[[.Assay`, asNamespace("SeuratDisk"))
# lockEnvironment(asNamespace("SeuratDisk"), bindings = TRUE)
# rm(`[[.Assay`)
