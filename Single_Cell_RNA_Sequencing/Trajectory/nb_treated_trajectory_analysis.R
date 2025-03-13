## NB trajectory analysis
## Treated samples

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
sknsh.cis.seurat <- readRDS("treated_data/sknsh_cisplatin.RDS")
sknsh.brd4i.seurat <- readRDS("treated_data/sknsh_brd4i.RDS")
sknsh.ezh2i.seurat <- readRDS("treated_data/sknsh_ezh2i.RDS")
shep.seurat <- readRDS("treated_data/shep_cisplatin.RDS")
shsy5y.seurat <- readRDS("treated_data/shsy5y_cisplatin.RDS")
nb039.seurat <- readRDS("treated_data/nb039_cisplatin.RDS")
nb067.seurat <- readRDS("treated_data/nb067_cisplatin.RDS")
grnb5.seurat <- readRDS("treated_data/grnb5_cisplatin.rds")

samples <- list(sknsh.cis.seurat, sknsh.brd4i.seurat, sknsh.ezh2i.seurat,
                shep.seurat, shsy5y.seurat,
                nb039.seurat, nb067.seurat, grnb5.seurat)
names(samples) <- c("SK-N-SH Cisplatin", "SK-N-SH BRD4i", "SK-N-SH EZH2i",
                    "SH-EP Cisplatin", "SH-SY5Y Cisplatin",
                    "NB039 Cisplatin", "NB067 Cisplatin", "GRNB5 Cisplatin")
rm(list = setdiff(ls(), c("umap.theme", "samples")))

## N.B. GRNB5 object does not have seurat_clusters.0.4 but seurat_clusters is equivalent (SH preprocessed)
samples[["GRNB5 Cisplatin"]]$seurat_clusters.0.4 <- samples[["GRNB5 Cisplatin"]]$seurat_clusters

for (object in 1:length(samples)) {
  umap.gg <- DimPlot(samples[[object]],
                      group.by = "seurat_clusters.0.4", order = TRUE, label = TRUE, repel = TRUE) +
    scale_colour_manual(values = iwanthue(length(unique(samples[[object]]$seurat_clusters.0.4)))) +
    umap.theme() + labs(title = paste0(names(samples)[object], ": Harmony Clusters (res = 0.4)")) +
    theme(legend.position = "none")
  umap.gg$layers[[2]]$aes_params$bg.colour <- "white"
  umap.gg$layers[[2]]$aes_params$direction <- "x"
  umap.gg$layers[[2]]$geom_params$max.overlaps <- 100
  umap.gg
  ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(samples)[object]))), "_umap_harmony_clusters.pdf"), width = 5.8, height = 5.8)
  ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(samples)[object]))), "_umap_harmony_clusters.png"), width = 5.8, height = 5.8)
  
  sce <- as.SingleCellExperiment(samples[[object]])
  reducedDim(sce, "PCA") <- Embeddings(samples[[object]], "pca")
  reducedDim(sce, "Harmony") <- Embeddings(samples[[object]], "Harmony")
  reducedDim(sce, "UMAP") <- Embeddings(samples[[object]], "umap")
  saveRDS(sce, paste0("treated_data/", gsub(" ", "_", gsub("-", "", tolower(names(samples)[object]))), "_sce.rds"))
}
## GRNB5 also has no Harmony batch correction so had to save as SCE object separately

## Run Slingshot
sknsh.cis.sce <- readRDS("treated_data/sknsh_cisplatin_sce.rds")
sknsh.cis.sling3 <- slingshot(sknsh.cis.sce, reducedDim = "PCA",
                              clusterLabels = "seurat_clusters.0.4", start.clus = 3)
saveRDS(sknsh.cis.sling3, "treated_data/sknsh_cisplatin_sce_sling3.rds")

sknsh.brd4i.sce <- readRDS("treated_data/sknsh_brd4i_sce.rds")
sknsh.brd4i.sling1 <- slingshot(sknsh.brd4i.sce, reducedDim = "PCA",
                                clusterLabels = "seurat_clusters.0.4", start.clus = 1)
saveRDS(sknsh.brd4i.sling1, "treated_data/sknsh_brd4i_sce_sling1.rds")

sknsh.ezh2i.sce <- readRDS("treated_data/sknsh_ezh2i_sce.rds")
sknsh.ezh2i.sling3 <- slingshot(sknsh.ezh2i.sce, reducedDim = "PCA",
                                clusterLabels = "seurat_clusters.0.4", start.clus = 3)
saveRDS(sknsh.ezh2i.sling3, "treated_data/sknsh_ezh2i_sce_sling3.rds")

shep.sce <- readRDS("treated_data/shep_cisplatin_sce.rds")
shep.sling4 <- slingshot(shep.sce, reducedDim = "PCA",
                         clusterLabels = "seurat_clusters.0.4", start.clus = 4)
saveRDS(shep.sling4, "treated_data/shep_cisplatin_sce_sling4.rds")

shsy5y.sce <- readRDS("treated_data/shsy5y_cisplatin_sce.rds")
shsy5y.sling2 <- slingshot(shsy5y.sce, reducedDim = "PCA",
                           clusterLabels = "seurat_clusters.0.4", start.clus = 2)
saveRDS(shsy5y.sling2, "treated_data/shsy5y_cisplatin_sce_sling2.rds")

nb039.sce <- readRDS("treated_data/nb039_cisplatin_sce.rds")
nb039.sling2 <- slingshot(nb039.sce, reducedDim = "PCA",
                          clusterLabels = "seurat_clusters.0.4", start.clus = 2)
saveRDS(nb039.sling2, "treated_data/nb039_cisplatin_sce_sling2.rds")
nb039.sling0 <- slingshot(nb039.sce, reducedDim = "PCA",
                          clusterLabels = "seurat_clusters.0.4", start.clus = 0)
saveRDS(nb039.sling0, "treated_data/nb039_cisplatin_sce_sling0.rds")

nb067.sce <- readRDS("treated_data/nb067_cisplatin_sce.rds")
nb067.sling0 <- slingshot(nb067.sce, reducedDim = "PCA",
                          clusterLabels = "seurat_clusters.0.4", start.clus = 0)
saveRDS(nb067.sling0, "treated_data/nb067_cisplatin_sce_sling0.rds")

grnb5.sce <- readRDS("treated_data/grnb5_cisplatin_sce.rds")
grnb5.sling4 <- slingshot(grnb5.sce, reducedDim = "PCA",
                          clusterLabels = "seurat_clusters.0.4", start.clus = 4)
saveRDS(grnb5.sling4, "treated_data/grnb5_cisplatin_sce_sling4.rds")

sknsh.cis.sling3 <- readRDS("treated_data/sknsh_cisplatin_sce_sling3.rds")
sknsh.brd4i.sling1 <- readRDS("treated_data/sknsh_brd4i_sce_sling1.rds")
sknsh.ezh2i.sling3 <- readRDS("treated_data/sknsh_ezh2i_sce_sling3.rds")
shep.sling4 <- readRDS("treated_data/shep_cisplatin_sce_sling4.rds")
shsy5y.sling2 <- readRDS("treated_data/shsy5y_cisplatin_sce_sling2.rds")
nb039.sling2 <- readRDS("treated_data/nb039_cisplatin_sce_sling2.rds")
nb039.sling0 <- readRDS("treated_data/nb039_cisplatin_sce_sling0.rds")
nb067.sling0 <- readRDS("treated_data/nb067_cisplatin_sce_sling0.rds")
grnb5.sling4 <- readRDS("treated_data/grnb5_cisplatin_sce_sling4.rds")

slings <- list(sknsh.cis.sling3, sknsh.brd4i.sling1, sknsh.ezh2i.sling3,
               shep.sling4, shsy5y.sling2,
               nb039.sling2, nb039.sling0, nb067.sling0, grnb5.sling4)
names(slings) <- c("SK-N-SH Cisplatin", "SK-N-SH BRD4i", "SK-N-SH EZH2i",
                   "SH-EP Cisplatin", "SH-SY5Y Cisplatin",
                   "NB039 Cisplatin (2)", "NB039 Cisplatin (0)", "NB067 Cisplatin", "GRNB5 Cisplatin")
rm(list = setdiff(ls(), c("umap.theme", "slings")))

for (object in 1:length(slings)) {
  slingshot.lineages <- slingPseudotime(slings[[object]])
  head(slingshot.lineages)
  shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)
  
  sling.gg <- plotUMAP(slings[[object]], colour_by = I(shared.pseudotime)) +
    scale_colour_viridis() +
    labs(colour = "Pseudotime")
  slingshot.embedded.all <- embedCurves(slings[[object]], "UMAP")
  slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
  for (lineage in slingshot.embedded.all) {
    slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
    sling.gg <- sling.gg + geom_path(data = slingshot.embedded.all,
                                     aes(x = umap_1, y = umap_2), linewidth = 1.2)
    ## N.B. "UMAP_1" and "UMAP_2" in GRNB5 object - had to run separately
  }
  saveRDS(slingshot.embedded.all, paste0("treated_data/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_embedded_all.rds"))
  
  sling.gg + labs(title = paste0(names(slings)[object], " Slingshot")) +
    theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))
  ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_umap_slingshot.pdf"), width = 8.7, height = 5.8)
  ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_umap_slingshot.png"), width = 8.7, height = 5.8)
  
  slingshot.embedded.all.sds <- embedCurves(slings[[object]], "UMAP")
  slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)
  saveRDS(slingshot.embedded.all.sds, paste0("treated_data/",  gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_sds_embedded_all.rds"))
  
  no.cols <- ceiling(ncol(slingshot.lineages)/3)
  pseudotime <- slingPseudotime(slings[[object]])
  names <- colnames(pseudotime)
  no.rows <- ceiling(length(names)/no.cols)
  pal <- viridis(100)
  pdf(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_pseudotime_lineages.pdf"),
      width = 8.3, height = 8.3)
  par(mfrow = c(no.rows, no.cols))
  for (i in names){
    cols <- pal[cut(pseudotime[,i], breaks = 100)]
    plot(reducedDim(slings[[object]], "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
    lines(slingshot.embedded.all.sds,
          lwd = 1, col = "black", type = "lineages", cex = 1)
  }
  dev.off()
  png(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_pseudotime_lineages.png"),
      width = 8.3, height = 8.3, units = "in", res = 200)
  par(mfrow = c(no.rows, no.cols))
  for (i in names){
    cols <- pal[cut(pseudotime[,i], breaks = 100)]
    plot(reducedDim(slings[[object]], "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
    lines(slingshot.embedded.all.sds,
          lwd = 1, col = "black", type = "lineages", cex = 1)
  }
  dev.off()
  
  slingshot.embedded <- embedCurves(slings[[object]], "UMAP")
  
  lineage.names <- colnames(slingshot.embedded)
  umap.clust.gg <- plotUMAP(slings[[object]], colour_by = "seurat_clusters.0.4") +
    umap.theme() + scale_colour_manual(values = iwanthue(length(unique(slings[[object]]$seurat_clusters.0.4)))) +
    theme(text = element_text(size = 15)) + guides(colour = guide_legend(override.aes = list(size = 5)))
  umap.amt.gg <- plotUMAP(slings[[object]], colour_by = "AMT.score") +
    umap.theme() + scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
    theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
  
  for (i in 1:length(lineage.names)){
    embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
    embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
    saveRDS(embedded.lineage, paste0("treated_data/slingshot_embedded_", lineage.names[[i]], ".rds"))
    
    plotUMAP(slings[[object]], colour_by = paste0("slingPseudotime_", i)) +
      geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
      labs(title = paste0(names(slings)[object], lineage.names[[i]])) + umap.theme() +
      theme(text = element_text(size = 15), legend.key.size = unit(1, "cm"))
    ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_", lineage.names[[i]], ".pdf"), width = 8.7, height = 5.8)
    ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_", lineage.names[[i]], ".png"), width = 8.7, height = 5.8)
    
    umap.clust.gg +
      geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
      labs(title = paste0(names(slings)[object], lineage.names[[i]]))
    ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_", lineage.names[[i]], "_cluster.pdf"), width = 8.7, height = 5.8)
    ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_", lineage.names[[i]], "_cluster.png"), width = 8.7, height = 5.8)
    
    umap.amt.gg +
      geom_path(data = embedded.lineage, aes(x = umap_1, y = umap_2), linewidth = 1.2) +
      labs(title = paste0(names(slings)[object], lineage.names[[i]]))
    ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_", lineage.names[[i]], "_amt.pdf"), width = 8.7, height = 5.8)
    ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(names(slings)[object]))), "_slingshot_", lineage.names[[i]], "_amt.png"), width = 8.7, height = 5.8)
  }
}

sknsh.cis.sling3 <- readRDS("treated_data/sknsh_cisplatin_sce_sling3.rds")
sknsh.brd4i.sling1 <- readRDS("treated_data/sknsh_brd4i_sce_sling1.rds")
sknsh.ezh2i.sling3 <- readRDS("treated_data/sknsh_ezh2i_sce_sling3.rds")
shep.sling4 <- readRDS("treated_data/shep_cisplatin_sce_sling4.rds")
shsy5y.sling2 <- readRDS("treated_data/shsy5y_cisplatin_sce_sling2.rds")
nb039.sling2 <- readRDS("treated_data/nb039_cisplatin_sce_sling2.rds")
nb067.sling0 <- readRDS("treated_data/nb067_cisplatin_sce_sling0.rds")
grnb5.sling4 <- readRDS("treated_data/grnb5_cisplatin_sce_sling4.rds")

slings <- list(sknsh.cis.sling3, sknsh.brd4i.sling1, sknsh.ezh2i.sling3,
               shep.sling4, shsy5y.sling2,
               nb039.sling2, nb067.sling0, grnb5.sling4)
names(slings) <- c("SK-N-SH Cisplatin", "SK-N-SH BRD4i", "SK-N-SH EZH2i",
                   "SH-EP Cisplatin", "SH-SY5Y Cisplatin",
                   "NB039 Cisplatin", "NB067 Cisplatin", "GRNB5 Cisplatin")
rm(list = setdiff(ls(), c("umap.theme", "slings")))

df.list <- list(as.data.frame(colData(slings[[1]])[, c("slingPseudotime_1", "slingPseudotime_2", "slingPseudotime_3", "AMT.score")]),
                as.data.frame(colData(slings[[2]])[, c("slingPseudotime_1", "slingPseudotime_2", "slingPseudotime_4", "AMT.score")]),
                as.data.frame(colData(slings[[3]])[, c("slingPseudotime_1", "slingPseudotime_2", "slingPseudotime_4", "AMT.score")]),
                as.data.frame(colData(slings[[4]])[, c("slingPseudotime_1", "slingPseudotime_3", "slingPseudotime_4", "AMT.score")]),
                as.data.frame(colData(slings[[5]])[, c("slingPseudotime_1", "AMT.score")]),
                as.data.frame(colData(slings[[6]])[, c("slingPseudotime_1", "slingPseudotime_2", "AMT.score")]),
                as.data.frame(colData(slings[[7]])[, c("slingPseudotime_1", "AMT.score")]),
                as.data.frame(colData(slings[[8]])[, c("slingPseudotime_1", "slingPseudotime_2", "slingPseudotime_3", "AMT.score")]))
names(df.list) <- names(slings)

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
      stat_cor(method = "pearson", label.x = x.pos, label.y = y.pos) +
      theme_bw() +
      scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
      labs(x = "Pseudotime",
           y = "AMT score",
           title = paste0(df.name, ": ", pseudotime, " vs AMT score"),
           colour = "AMT score") +
      theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

    ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(df.name))),"_scatter_", gsub("_", "", tolower(pseudotime)), ".pdf"), width = 8.7, height = 5.8)
    ggsave(paste0("treated_plots/", gsub(" ", "_", gsub("-", "", tolower(df.name))),"_scatter_", gsub("_", "", tolower(pseudotime)), ".png"), width = 8.7, height = 5.8)
  }
}

## [ Dynamic genes ] ----

## /////////////////////////////////////////////////////////////////////////////
## SK-N-SH Cisplatin ///////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

sknsh.cis.sling3 <- readRDS("treated_data/sknsh_cisplatin_sce_sling3.rds")
sknsh.cols <- grep("^slingPseudotime_", colnames(colData(sknsh.cis.sling3)), value = TRUE)
cluster.cols <- iwanthue(length(unique(sknsh.cis.sling3$seurat_clusters.0.4)))
names(cluster.cols) <- levels(sknsh.cis.sling3$seurat_clusters.0.4)

for (pseudotime in sknsh.cols) {
  lineage <- testPseudotime(sknsh.cis.sling3, pseudotime = sknsh.cis.sling3[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("treated_data/trajectory_genes/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("treated_data/trajectory_genes/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(sknsh.cis.sling3, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(sknsh.cis.sling3[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(sknsh.cis.sling3[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("treated_data/trajectory_genes/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("treated_data/trajectory_genes/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(sknsh.cis.sling3, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(sknsh.cis.sling3[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("treated_plots/sknsh_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## SK-N-SH BRD4i ///////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

sknsh.brd4i.sling1 <- readRDS("treated_data/sknsh_brd4i_sce_sling1.rds")
sknsh.cols <- grep("^slingPseudotime_", colnames(colData(sknsh.brd4i.sling1)), value = TRUE)
cluster.cols <- iwanthue(length(unique(sknsh.brd4i.sling1$seurat_clusters.0.4)))
names(cluster.cols) <- levels(sknsh.brd4i.sling1$seurat_clusters.0.4)

for (pseudotime in sknsh.cols) {
  lineage <- testPseudotime(sknsh.brd4i.sling1, pseudotime = sknsh.brd4i.sling1[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("treated_data/trajectory_genes/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("treated_data/trajectory_genes/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(sknsh.brd4i.sling1, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(sknsh.brd4i.sling1[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(sknsh.brd4i.sling1[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("treated_data/trajectory_genes/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("treated_data/trajectory_genes/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(sknsh.brd4i.sling1, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(sknsh.brd4i.sling1[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("treated_plots/sknsh_brd4i_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## SK-N-SH EZH2i ///////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

sknsh.ezh2i.sling3 <- readRDS("treated_data/sknsh_ezh2i_sce_sling3.rds")
sknsh.cols <- grep("^slingPseudotime_", colnames(colData(sknsh.ezh2i.sling3)), value = TRUE)
cluster.cols <- iwanthue(length(unique(sknsh.ezh2i.sling3$seurat_clusters.0.4)))
names(cluster.cols) <- levels(sknsh.ezh2i.sling3$seurat_clusters.0.4)

for (pseudotime in sknsh.cols) {
  lineage <- testPseudotime(sknsh.ezh2i.sling3, pseudotime = sknsh.ezh2i.sling3[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("treated_data/trajectory_genes/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("treated_data/trajectory_genes/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(sknsh.ezh2i.sling3, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(sknsh.ezh2i.sling3[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(sknsh.ezh2i.sling3[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("treated_data/trajectory_genes/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("treated_data/trajectory_genes/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(sknsh.ezh2i.sling3, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(sknsh.ezh2i.sling3[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("treated_plots/sknsh_ezh2i_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## SH-EP Cisplatin /////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

shep.sling4 <- readRDS("treated_data/shep_cisplatin_sce_sling4.rds")
shep.cols <- grep("^slingPseudotime_", colnames(colData(shep.sling4)), value = TRUE)
cluster.cols <- iwanthue(length(unique(shep.sling4$seurat_clusters.0.4)))
names(cluster.cols) <- levels(shep.sling4$seurat_clusters.0.4)

for (pseudotime in shep.cols) {
  lineage <- testPseudotime(shep.sling4, pseudotime = shep.sling4[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("treated_data/trajectory_genes/shep_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("treated_data/trajectory_genes/shep_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(shep.sling4, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(shep.sling4[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(shep.sling4[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("treated_data/trajectory_genes/shep_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("treated_data/trajectory_genes/shep_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(shep.sling4, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(shep.sling4[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("treated_plots/shep_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## SH-SY5Y Cisplatin ///////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

shsy5y.sling2 <- readRDS("treated_data/shsy5y_cisplatin_sce_sling2.rds")
shsy5y.cols <- grep("^slingPseudotime_", colnames(colData(shsy5y.sling2)), value = TRUE)
cluster.cols <- iwanthue(length(unique(shsy5y.sling2$seurat_clusters.0.4)))
names(cluster.cols) <- levels(shsy5y.sling2$seurat_clusters.0.4)

for (pseudotime in shsy5y.cols) {
  lineage <- testPseudotime(shsy5y.sling2, pseudotime = shsy5y.sling2[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("treated_data/trajectory_genes/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("treated_data/trajectory_genes/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(shsy5y.sling2, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(shsy5y.sling2[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(shsy5y.sling2[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("treated_data/trajectory_genes/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("treated_data/trajectory_genes/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(shsy5y.sling2, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(shsy5y.sling2[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("treated_plots/shsy5y_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## NB039 Cisplatin /////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

nb039.sling2 <- readRDS("treated_data/nb039_cisplatin_sce_sling2.rds")
nb039.cols <- grep("^slingPseudotime_", colnames(colData(nb039.sling2)), value = TRUE)
cluster.cols <- iwanthue(length(unique(nb039.sling2$seurat_clusters.0.4)))
names(cluster.cols) <- levels(nb039.sling2$seurat_clusters.0.4)

for (pseudotime in nb039.cols) {
  lineage <- testPseudotime(nb039.sling2, pseudotime = nb039.sling2[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("treated_data/trajectory_genes/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("treated_data/trajectory_genes/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(nb039.sling2, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(nb039.sling2[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(nb039.sling2[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("treated_data/trajectory_genes/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("treated_data/trajectory_genes/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(nb039.sling2, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(nb039.sling2[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("treated_plots/nb039_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## NB067 Cisplatin /////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

nb067.sling0 <- readRDS("treated_data/nb067_cisplatin_sce_sling0.rds")
nb067.cols <- grep("^slingPseudotime_", colnames(colData(nb067.sling0)), value = TRUE)
cluster.cols <- iwanthue(length(unique(nb067.sling0$seurat_clusters.0.4)))
names(cluster.cols) <- levels(nb067.sling0$seurat_clusters.0.4)

for (pseudotime in nb067.cols) {
  lineage <- testPseudotime(nb067.sling0, pseudotime = nb067.sling0[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("treated_data/trajectory_genes/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("treated_data/trajectory_genes/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(nb067.sling0, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(nb067.sling0[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(nb067.sling0[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("treated_data/trajectory_genes/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("treated_data/trajectory_genes/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(nb067.sling0, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(nb067.sling0[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("treated_plots/nb067_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## /////////////////////////////////////////////////////////////////////////////
## GRNB5 Cisplatin /////////////////////////////////////////////////////////////
## /////////////////////////////////////////////////////////////////////////////

grnb5.sling4 <- readRDS("treated_data/grnb5_cisplatin_sce_sling4.rds")
grnb5.cols <- grep("^slingPseudotime_", colnames(colData(grnb5.sling4)), value = TRUE)
cluster.cols <- iwanthue(length(unique(grnb5.sling4$seurat_clusters.0.4)))
names(cluster.cols) <- levels(grnb5.sling4$seurat_clusters.0.4)

for (pseudotime in grnb5.cols) {
  lineage <- testPseudotime(grnb5.sling4, pseudotime = grnb5.sling4[[pseudotime]])
  lineage.genes <- lineage[order(lineage$p.value), ]
  
  lineage.genes.down <- as.data.frame(lineage.genes[lineage.genes$logFC < 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(logFC)
  saveRDS(lineage.genes.down, paste0("treated_data/trajectory_genes/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.rds"))
  write.csv(lineage.genes.down, paste0("treated_data/trajectory_genes/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)),"_down_genes.csv"), row.names = FALSE)
  
  lin.down.gg <- plotExpression(grnb5.sling4, features = head(lineage.genes.down$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": downregulated"))
  lin.down.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top10.png"), width = 11.6, height = 8.7)
  
  on.lineage <- !is.na(grnb5.sling4[[pseudotime]])
  lin.down.heatmap <- plotHeatmap(grnb5.sling4[, on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                  column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                  features = head(lineage.genes.down$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.down.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.down.heatmap)
  dev.off()
  png(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_down_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.down.heatmap)
  dev.off()
  
  lineage.genes.up <- as.data.frame(lineage.genes[lineage.genes$logFC > 0, ]) %>%
    filter(p.value < 0.05 & FDR < 0.05) %>%
    mutate(gene = rownames(.)) %>%
    arrange(desc(logFC))
  head(lineage.genes.up, 10)
  saveRDS(lineage.genes.up, paste0("treated_data/trajectory_genes/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.rds"))
  write.csv(lineage.genes.up, paste0("treated_data/trajectory_genes/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)),"_up_genes.csv"), row.names = FALSE)
  
  lin.up.gg <- plotExpression(grnb5.sling4, features = head(lineage.genes.up$gene, 10), x = pseudotime, colour_by = "seurat_clusters.0.4") +
    scale_colour_manual(values = cluster.cols) +
    guides(colour = guide_legend(title = "Cluster")) +
    labs(x = "Pseudotime",
         y = "Expression",
         title = paste0(pseudotime, ": upregulated"))
  lin.up.gg +
    theme(text = element_text(size = 20)) +
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.pdf"), width = 11.6, height = 8.7)
  ggsave(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top10.png"), width = 11.6, height = 8.7)
  
  lin.up.gg + lin.down.gg + plot_layout(guides = "collect") +
    plot_annotation(title = pseudotime) &
    theme(text = element_text(size = 20)) &
    guides(colour = guide_legend(override.aes = list(size = 10), title = ""))
  ggsave(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.pdf"), width = 18.8, height = 8.7)
  ggsave(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_expression.png"), width = 18.8, height = 8.7)
  
  lin.up.heatmap <- plotHeatmap(grnb5.sling4[ ,on.lineage], order_columns_by = pseudotime, colour_columns_by = "seurat_clusters.0.4",
                                column_annotation_colours = list(seurat_clusters.0.4 = cluster.cols),
                                features = head(lineage.genes.up$gene, 50), center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))
  lin.up.heatmap$gtable$grobs[[5]]$label <- c("Cluster", "Pseudotime")
  pdf(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.pdf"), width = 8.7, height = 8.7)
  print(lin.up.heatmap)
  dev.off()
  png(paste0("treated_plots/grnb5_cisplatin_", gsub("_", "", tolower(pseudotime)), "_up_top50.png"), width = 8.7, height = 8.7,
      units = "in", res = 200)
  print(lin.up.heatmap)
  dev.off()
}

## [ Gene list similarity ] ----

dir.path <- "treated_data/trajectory_genes"
rds.files <- list.files(dir.path, pattern = "\\.rds$", full.names = TRUE)

up.list <- list()
down.list <- list()

for (file in rds.files) {
  file.name <- basename(file)
  split <- strsplit(file.name, "_")[[1]]
  sample <- toupper(split[1])
  drug <- split[2]
  prefix <- paste0(sample, drug)
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

# up.list <- up.list[c("SKNSHcisplatin_1", "SKNSHcisplatin_2", "SKNSHcisplatin_3",
#                      "SKNSHbrd4i_1", "SKNSHbrd4i_2", "SKNSHbrd4i_4",
#                      "SKNSHezh2i_1", "SKNSHezh2i_2", "SKNSHezh2i_4",
#                      "SHEPcisplatin_1", "SHEPcisplatin_3", "SHEPcisplatin_4",
#                      "SKSY5Ycisplatin_1",
#                      "NB039cisplatin_1", "NB039cisplatin_2",
#                      "NB067cisplatin_1",
#                      "GRNB5cisplatin_1", "GRNB5cisplatin_2", "GRNB5cisplatin_3")]
# down.list <- down.list[c("SKNSHcisplatin_1", "SKNSHcisplatin_2", "SKNSHcisplatin_3",
#                          "SKNSHbrd4i_1", "SKNSHbrd4i_2", "SKNSHbrd4i_4",
#                          "SKNSHezh2i_1", "SKNSHezh2i_2", "SKNSHezh2i_4",
#                          "SHEPcisplatin_1", "SHEPcisplatin_3", "SHEPcisplatin_4",
#                          "SKSY5Ycisplatin_1",
#                          "NB039cisplatin_1", "NB039cisplatin_2",
#                          "NB067cisplatin_1",
#                          "GRNB5cisplatin_1", "GRNB5cisplatin_2", "GRNB5cisplatin_3")]

up.list <- up.list[c("SKNSHcisplatin_3", "SKNSHbrd4i_1", "SKNSHezh2i_4",
                     "SHEPcisplatin_3", "SKSY5Ycisplatin_1",
                     "NB039cisplatin_2", "NB067cisplatin_1", "GRNB5cisplatin_3")]
down.list <- down.list[c("SKNSHcisplatin_3", "SKNSHbrd4i_1", "SKNSHezh2i_4",
                         "SHEPcisplatin_3", "SKSY5Ycisplatin_1",
                         "NB039cisplatin_2", "NB067cisplatin_1", "GRNB5cisplatin_3")]

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
# ggsave("treated_plots/genes_comparison_up_heatmap.pdf", width = 8.3, height = 8.3)
# ggsave("treated_plots/genes_comparison_up_heatmap.png", width = 8.3, height = 8.3)
ggsave("treated_plots/genes_comparison_up_heatmap_examples.pdf", width = 5.8, height = 5.8)
ggsave("treated_plots/genes_comparison_up_heatmap_examples.png", width = 5.8, height = 5.8)

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
# ggsave("treated_plots/genes_comparison_down_heatmap.pdf", width = 8.3, height = 8.3)
# ggsave("treated_plots/genes_comparison_down_heatmap.png", width = 8.3, height = 8.3)
ggsave("treated_plots/genes_comparison_down_heatmap_examples.pdf", width = 5.8, height = 5.8)
ggsave("treated_plots/genes_comparison_down_heatmap_examples.png", width = 5.8, height = 5.8)

## [ Seurat to AnnData ] ----

sknsh.cis.seurat <- readRDS("treated_data/sknsh_cisplatin.RDS")
sknsh.brd4i.seurat <- readRDS("treated_data/sknsh_brd4i.RDS")
sknsh.ezh2i.seurat <- readRDS("treated_data/sknsh_ezh2i.RDS")
shep.seurat <- readRDS("treated_data/shep_cisplatin.RDS")
shsy5y.seurat <- readRDS("treated_data/shsy5y_cisplatin.RDS")
nb039.seurat <- readRDS("treated_data/nb039_cisplatin.RDS")
nb067.seurat <- readRDS("treated_data/nb067_cisplatin.RDS")
grnb5.seurat <- readRDS("treated_data/grnb5_cisplatin.rds")

samples <- list(sknsh.cis.seurat, sknsh.brd4i.seurat, sknsh.ezh2i.seurat,
                shep.seurat, shsy5y.seurat,
                nb039.seurat, nb067.seurat, grnb5.seurat)
names(samples) <- c("SK-N-SH Cisplatin", "SK-N-SH BRD4i", "SK-N-SH EZH2i",
                    "SH-EP Cisplatin", "SH-SY5Y Cisplatin",
                    "NB039 Cisplatin", "NB067 Cisplatin", "GRNB5 Cisplatin")
rm(list = setdiff(ls(), "samples"))

for (object in 1:length(samples)) {
  DefaultAssay(samples[[object]]) <- "RNA"
  samples[[object]] <- NormalizeData(samples[[object]])
  samples[[object]] <- FindVariableFeatures(samples[[object]], nfeatures = 1000)
  samples[[object]] <- ScaleData(samples[[object]])
  samples[[object]][["RNA"]] <- as(samples[[object]][["RNA"]], "Assay")
  samples[[object]][["RNA"]]$data <- NULL
  samples[[object]][["RNA"]]$scale.data <- NULL
  
  SaveH5Seurat(samples[[object]], filename = paste0("treated_data/", gsub(" ", "_", gsub("-", "", tolower(names(samples)[object]))),"_seurat.h5Seurat"), overwrite = T, verbose = T)
  Convert(paste0("treated_data/", gsub(" ", "_", gsub("-", "", tolower(names(samples)[object]))),"_seurat.h5Seurat"), dest = "h5ad", assay = "RNA", overwrite = T)
}

## Split data by condition and save as AnnData objects
split.samples <- list()
for (object in 1:length(samples)) {
  sample.name <- names(samples[object])
  split <- SplitObject(samples[[object]], split.by = "Condition")
  names(split) <- paste0(sample.name, "_", names(split))
  split.samples <- append(split.samples, split)
}

split.filt <- split.samples[!grepl("Untreated|POT", names(split.samples), ignore.case = TRUE)]
names(split.filt)
## "SK-N-SH Cisplatin_Cisplatin_1weekOFF"  "SK-N-SH Cisplatin_Cisplatin_ON"        "SK-N-SH Cisplatin_Cisplatin_4weeksOFF"
## "SK-N-SH BRD4i_JQ1_OFF"                 "SK-N-SH BRD4i_JQ1_ON"                  "SK-N-SH EZH2i_EZH2i_OFF"              
## "SK-N-SH EZH2i_EZH2i_ON"                "SH-EP Cisplatin_cisplatin recovery"    "SH-EP Cisplatin_cisplatin"            
## "SH-SY5Y Cisplatin_cisplatin recovery"  "SH-SY5Y Cisplatin_cisplatin"           "NB039 Cisplatin_cisplatin"            
## "NB039 Cisplatin_cisplatin recovery"    "NB067 Cisplatin_cisplatin"             "NB067 Cisplatin_cisplatin recovery"   
## "GRNB5 Cisplatin_cisplatin"             "GRNB5 Cisplatin_cisplatin recovery"   
names(split.filt) <- c("SK-N-SH Cisplatin 1weekOFF", "SK-N-SH Cisplatin_Cisplatin_ON", "SK-N-SH Cisplatin 4weeksOFF",
                       "SK-N-SH BRD4i OFF", "SK-N-SH BRD4i ON", "SK-N-SH EZH2i OFF", "SK-N-SH EZH2i ON",
                       "SH-EP Cisplatin OFF", "SH-EP Cisplatin ON", "SH-SY5Y Cisplatin OFF", "SH-SY5Y Cisplatin ON",
                       "NB039 Cisplatin ON", "NB039 Cisplatin OFF",
                       "NB067 Cisplatin ON", "NB067 Cisplatin OFF",
                       "GRNB5 Cisplatin ON", "GRNB5 Cisplatin OFF") 
rm(list = setdiff(ls(), "split.filt"))

for (object in 1:length(split.filt)) {
  DefaultAssay(split.filt[[object]]) <- "RNA"
  split.filt[[object]] <- NormalizeData(split.filt[[object]])
  split.filt[[object]] <- FindVariableFeatures(split.filt[[object]], nfeatures = 1000)
  split.filt[[object]] <- ScaleData(split.filt[[object]])
  split.filt[[object]][["RNA"]] <- as(split.filt[[object]][["RNA"]], "Assay")
  split.filt[[object]][["RNA"]]$data <- NULL
  split.filt[[object]][["RNA"]]$scale.data <- NULL
  
  SaveH5Seurat(split.filt[[object]], filename = paste0("treated_data_split/", gsub(" ", "_", gsub("-", "", tolower(names(split.filt)[object]))),"_seurat.h5Seurat"), overwrite = T, verbose = T)
  Convert(paste0("treated_data_split/", gsub(" ", "_", gsub("-", "", tolower(names(split.filt)[object]))),"_seurat.h5Seurat"), dest = "h5ad", assay = "RNA", overwrite = T)
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
