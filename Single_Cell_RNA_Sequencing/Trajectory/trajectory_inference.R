## Trajectory inference

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(scater)
library(scran)
library(tidyverse)
library(TSCAN)
library(ggplot2)
library(slingshot)

nb.sce <- readRDS("nb_sce_AMT.rds")

identity.cols <- c("#CC6677", "#228833", "#0077BB",
                   "#332288", "#66CCEE", "#009988",
                   "#EECC66", "#CC3311", "#CCDDAA")
names(identity.cols) <- levels(nb.sce$seurat_clusters.0.2)

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

## See: https://bioconductor.org/books/3.14/OSCA.advanced/trajectory-analysis.html

nb.slingshot.sce <- slingshot(nb.sce, reducedDim = "PCA",
                              clusterLabels = "seurat_clusters.0.8", start.clus = 4)
## 8 lineages
saveRDS(nb.slingshot.sce, "nb_slingshot_sce.rds")

slingshot.lineages <- slingPseudotime(nb.slingshot.sce)
head(slingshot.lineages)

shared.pseudotime <- rowMeans(slingshot.lineages, na.rm = TRUE)

## Extended Data Figure 7a
require(viridis)
umap.gg <- plotUMAP(nb.slingshot.sce, colour_by = I(shared.pseudotime)) +
  scale_colour_viridis() +
  labs(colour = "Pseudotime")
slingshot.embedded.all <- embedCurves(nb.slingshot.sce, "UMAP")
slingshot.embedded.all <- slingCurves(slingshot.embedded.all)
for (lineage in slingshot.embedded.all) {
  slingshot.embedded.all <- data.frame(lineage$s[lineage$ord,])
  umap.gg <- umap.gg + geom_path(data = slingshot.embedded.all,
                                   aes(x = UMAP_1, y=UMAP_2), linewidth = 1.2)
}

saveRDS(slingshot.embedded.all, "slingshot_embedded_all.rds")

umap.gg + theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))

slingshot.embedded.all.sds <- embedCurves(nb.slingshot.sce, "UMAP")
slingshot.embedded.all.sds <- SlingshotDataSet(slingshot.embedded.all.sds)

no.cols <- 3
pseudotime <- slingPseudotime(nb.slingshot.sce)
names <- colnames(pseudotime)
no.rows <- ceiling(length(names)/no.cols)
require(viridis)
pal <- viridis(100)
par(mfrow = c(no.rows, no.cols))
for (i in names){
  cols <- pal[cut(pseudotime[,i], breaks = 100)]
  plot(reducedDim(nb.slingshot.sce, "UMAP"), col = cols, pch = 16, cex = 0.5, main = i)
  lines(slingshot.embedded.all.sds,
        lwd = 1, col = "black", type = "lineages", cex = 1)
}

## Embedding lineages one by one
slingshot.embedded <- embedCurves(nb.slingshot.sce, "UMAP")

lineage.names <- colnames(slingshot.embedded)
umap.clust.gg <- plotUMAP(nb.slingshot.sce, colour_by = "seurat_clusters.0.2") +
  scale_colour_manual(values = identity.cols)
umap.amt.gg <- plotUMAP(nb.slingshot.sce, colour_by = "AMT.score") +
  scale_colour_gradientn(colours = pals::ocean.thermal(100))

for (i in 1:length(lineage.names)){
  embedded.lineage <- slingCurves(slingshot.embedded)[[i]]
  embedded.lineage <- data.frame(embedded.lineage$s[embedded.lineage$ord, ])
  saveRDS(embedded.lineage, paste0("slingshot_embedded_", lineage.names[[i]], ".rds"))
  
  plotUMAP(nb.slingshot.sce, colour_by = paste0("slingPseudotime_", i)) +
    geom_path(data = embedded.lineage, aes(x = UMAP_1, y = UMAP_2), size = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("sce_umap_pca_slingshot_", lineage.names[[i]], ".pdf"), width = 8.7, height = 8.7)
  ggsave(paste0("sce_umap_pca_slingshot_", lineage.names[[i]], ".png"), width = 8.7, height = 8.7)
  
  umap.clust.gg +
    geom_path(data = embedded.lineage, aes(x = UMAP_1, y = UMAP_2), size = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("sce_umap_pca_slingshot_", lineage.names[[i]], "_cluster.pdf"), width = 8.7, height = 8.7)
  ggsave(paste0("sce_umap_pca_slingshot_", lineage.names[[i]], "_cluster.png"), width = 8.7, height = 8.7)
  
  umap.amt.gg +
    geom_path(data = embedded.lineage, aes(x = UMAP_1, y = UMAP_2), size = 1.2) +
    labs(title = lineage.names[[i]])
  ggsave(paste0("sce_umap_pca_slingshot_", lineage.names[[i]], "_AMT.pdf"), width = 8.7, height = 8.7)
  ggsave(paste0("sce_umap_pca_slingshot_", lineage.names[[i]], "_AMT.png"), width = 8.7, height = 8.7)
}
## Lineage 1: CyclingADRN - EMT - JQ1
## Lineage 2: CyclingADRN - EMT - QC
## Lineage 3: CyclingADRN - EMT - Cisplatin

jq1.lineage <- readRDS("slingshot_embedded_Lineage1.rds")

## Figure 5a
jq1.gg <- umap.clust.gg +
  geom_path(data = jq1.lineage, aes(x = UMAP_1, y = UMAP_2), size = 1.2) +
  labs(title = "JQ1 treatment lineage", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

jq1.gg

jq1.amt.gg <- umap.amt.gg +
  geom_path(data = jq1.lineage, aes(x = UMAP_1, y = UMAP_2), size = 1.2) +
  labs(title = "JQ1 treatment lineage", colour = "AMT score") + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

jq1.amt.gg

stress.lineage <- readRDS("slingshot_embedded_Lineage2.rds")

## Figure 5a
stress.gg <- umap.clust.gg +
  geom_path(data = stress.lineage, aes(x = UMAP_1, y = UMAP_2), size = 1.2) +
  labs(title = "Cisplatin recovery lineage", colour = "") + 
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

stress.gg

stress.amt.gg <- umap.amt.gg +
  geom_path(data = stress.lineage, aes(x = UMAP_1, y = UMAP_2), size = 1.2) +
  labs(title = "Cisplatin recovery lineage", colour = "AMT score") + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

stress.amt.gg

cis.lineage <- readRDS("slingshot_embedded_Lineage3.rds")

## Figure 5a
cis.gg <- umap.clust.gg +
  geom_path(data = cis.lineage, aes(x = UMAP_1, y = UMAP_2), size = 1.2) +
  labs(title = "Cisplatin treatment lineage", colour = "") +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

cis.gg

cis.amt.gg <- umap.amt.gg +
  geom_path(data = cis.lineage, aes(x = UMAP_1, y = UMAP_2), size = 1.2) +
  labs(title = "Cisplatin treatment lineage", colour = "AMT score") + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

cis.amt.gg

## Pseudotime-AMT correlation
nb.slingshot.df <- as.data.frame(colData(nb.slingshot.sce)[, c("slingPseudotime_1", "slingPseudotime_2",
                                                               "slingPseudotime_3", "AMT.score")])

require(scattermore)
require(ggpubr)
## Extended Data Figure 7b
lin1.gg <- ggplot(nb.slingshot.df,
                  aes(x = slingPseudotime_1, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 2) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 100, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Pseudotime",
       y = "AMT score",
       title = "JQ1 treatment lineage",
       colour = "AMT score")

lin1.gg + 
  theme(text = element_text(size = 20), legend.key.size = unit(1, "cm"))

lin2.gg <- ggplot(nb.slingshot.df,
                  aes(x = slingPseudotime_2, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 1) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 130, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Slingshot pseudotime",
       y = "AMT score",
       title = "Lineage 2",
       colour = "AMT score")

lin2.gg

lin3.gg <- ggplot(nb.slingshot.df,
                  aes(x = slingPseudotime_3, y = AMT.score, colour = AMT.score)) +
  geom_scattermore(pointsize = 1) +
  geom_smooth(method = "lm", colour = "black") +
  stat_cor(method = "pearson", label.x = 115, label.y = -15) +
  theme_bw() +
  scale_colour_gradientn(colours = pals::ocean.thermal(100)) +
  labs(x = "Slingshot pseudotime",
       y = "AMT score",
       title = "Lineage 3",
       colour = "AMT score")

lin3.gg

## [ Gene dynamics: linear model ] ----

## Linear model fitting for lineage 1 (JQ1)
lineage1 <- testPseudotime(nb.slingshot.sce, pseudotime = nb.slingshot.sce$slingPseudotime_1)
lineage1.genes <- lineage1[order(lineage1$p.value), ]
head(lineage1.genes, 10)

## Genes that decrease in expression with increasing pseudotime along lineage 1
lineage1.genes.down <- as.data.frame(lineage1.genes[lineage1.genes$logFC < 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage1.genes.down, 10)

write.csv(lineage1.genes.down,"dynamic genes/lineage1_down_genes.csv", row.names = FALSE)

## Extended Data Figure 7g
lin1.down.gg <- plotExpression(nb.slingshot.sce, features = head(lineage1.genes.down$gene, 10), x = "slingPseudotime_1", colour_by = "seurat_clusters.0.2") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "JQ1 recovery-specific: downregulated")

lin1.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

on.lineage1 <- !is.na(nb.slingshot.sce$slingPseudotime_1)
lin1.heatmap.down <- plotHeatmap(nb.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "seurat_clusters.0.2",
                                 column_annotation_colours = list(seurat_clusters.0.2 = identity.cols), features = head(lineage1.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin1.heatmap.down

## Genes that increase in expression with increasing pseudotime along lineage 1
lineage1.genes.up <- as.data.frame(lineage1.genes[lineage1.genes$logFC > 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage1.genes.up, 10)

write.csv(lineage1.genes.up,"dynamic genes/lineage1_up_genes.csv", row.names = FALSE)

## Extended Data Figure 7h
lin1.up.gg <- plotExpression(nb.slingshot.sce, features = head(lineage1.genes.up$gene, 10), x = "slingPseudotime_1", colour_by = "seurat_clusters.0.2") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "JQ1 recovery-specific: upregulated")

lin1.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

lin1.heatmap.up <- plotHeatmap(nb.slingshot.sce[ ,on.lineage1], order_columns_by = "slingPseudotime_1", colour_columns_by = "seurat_clusters.0.2",
                               column_annotation_colours = list(seurat_clusters.0.2 = identity.cols), features = head(lineage1.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin1.heatmap.up

## Linear model fitting for lineage 2 (Stress)
lineage2 <- testPseudotime(nb.slingshot.sce, pseudotime = nb.slingshot.sce$slingPseudotime_2)
lineage2.genes <- lineage2[order(lineage2$p.value), ]
head(lineage2.genes, 10)

## Genes that decrease in expression with increasing pseudotime along lineage 2
lineage2.genes.down <- as.data.frame(lineage2.genes[lineage2.genes$logFC < 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage2.genes.down, 10)

write.csv(lineage2.genes.down,"dynamic genes/lineage2_down_genes.csv", row.names = FALSE)

## Extended Data Figure 7e
lin2.down.gg <- plotExpression(nb.slingshot.sce, features = head(lineage2.genes.down$gene, 10), x = "slingPseudotime_3", colour_by = "seurat_clusters.0.2") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Cisplatin recovery-specific: downregulated")

lin2.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = "Cluster"))

on.lineage2 <- !is.na(nb.slingshot.sce$slingPseudotime_2)
lin2.heatmap.down <- plotHeatmap(nb.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "seurat_clusters.0.2",
                                 column_annotation_colours = list(seurat_clusters.0.2 = identity.cols), features = head(lineage2.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.heatmap.down

## Genes that increase in expression with increasing pseudotime along lineage 2
lineage2.genes.up <- as.data.frame(lineage2.genes[lineage2.genes$logFC > 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage2.genes.up, 10)

write.csv(lineage2.genes.up,"dynamic genes/lineage2_up_genes.csv", row.names = FALSE)

## Extended Data Figure 7f
lin2.up.gg <- plotExpression(nb.slingshot.sce, features = head(lineage2.genes.up$gene, 10), x = "slingPseudotime_2", colour_by = "seurat_clusters.0.2") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Cisplatin recovery-specific: upregulated")

lin2.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = "Cluster"))

lin2.heatmap.up <- plotHeatmap(nb.slingshot.sce[ ,on.lineage2], order_columns_by = "slingPseudotime_2", colour_columns_by = "seurat_clusters.0.2",
                               column_annotation_colours = list(seurat_clusters.0.2 = identity.cols), features = head(lineage2.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin2.heatmap.up

## Linear model fitting for lineage 3 (cisplatin)
lineage3 <- testPseudotime(nb.slingshot.sce, pseudotime = nb.slingshot.sce$slingPseudotime_3)
lineage3.genes <- lineage3[order(lineage3$p.value), ]
head(lineage3.genes, 10)

## Genes that decrease in expression with increasing pseudotime along lineage 3
lineage3.genes.down <- as.data.frame(lineage3.genes[lineage3.genes$logFC < 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(logFC)
head(lineage3.genes.down, 10)

write.csv(lineage3.genes.down,"dynamic genes/lineage3_down_genes.csv", row.names = FALSE)

## Extended Data Figure 7c
lin3.down.gg <- plotExpression(nb.slingshot.sce, features = head(lineage3.genes.down$gene, 10), x = "slingPseudotime_3", colour_by = "seurat_clusters.0.2") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Cisplatin treatment-specific: downregulated")

lin3.down.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

on.lineage3 <- !is.na(nb.slingshot.sce$slingPseudotime_3)
lin3.heatmap.down <- plotHeatmap(nb.slingshot.sce[ ,on.lineage3], order_columns_by = "slingPseudotime_3", colour_columns_by = "seurat_clusters.0.2",
                                 column_annotation_colours = list(seurat_clusters.0.2 = identity.cols), features = head(lineage3.genes.down$gene, 50),
                                 center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin3.heatmap.down 

## Genes that increase in expression with increasing pseudotime along lineage 3
lineage3.genes.up <- as.data.frame(lineage3.genes[lineage3.genes$logFC > 0, ]) %>%
  filter(p.value < 0.05 & FDR < 0.05) %>%
  mutate(gene = rownames(.)) %>%
  arrange(desc(logFC))
head(lineage3.genes.up, 10)

write.csv(lineage3.genes.up,"dynamic genes/lineage3_up_genes.csv", row.names = FALSE)

## Extended Data Figure 7d
lin3.up.gg <- plotExpression(nb.slingshot.sce, features = head(lineage3.genes.up$gene, 10), x = "slingPseudotime_3", colour_by = "seurat_clusters.0.2") +
  scale_colour_manual(values = identity.cols) +
  guides(colour = guide_legend(title = "Cluster")) +
  labs(x = "Pseudotime",
       y = "Expression",
       title = "Cisplatin treatment-specific: upregulated")

lin3.up.gg +
  theme(text = element_text(size = 20)) +
  guides(colour = guide_legend(override.aes = list(size = 10), title = ""))

lin3.heatmap.up <- plotHeatmap(nb.slingshot.sce[ ,on.lineage3], order_columns_by = "slingPseudotime_3", colour_columns_by = "seurat_clusters.0.2",
                               column_annotation_colours = list(seurat_clusters.0.2 = identity.cols), features = head(lineage3.genes.up$gene, 50),
                               center = TRUE, zlim = c(-3, 3), colour = pals::ocean.balance(100))

lin3.heatmap.up

## [ Gene dynamics: random forest ] ----

## See: https://bustools.github.io/BUS_notebooks_R/slingshot.html

nb.slingshot.sce <- readRDS("nb_slingshot_sce.rds")

## Get top highly variable genes
gene.var <- modelGeneVar(nb.slingshot.sce)
top.hvg <- getTopHVGs(gene.var, n = 300)
str(top.hvg)

logcounts.hvg <- t(logcounts(nb.slingshot.sce)[top.hvg,])

## Prepare lineage 1 (JQ1) data for random forest
lin1.data <- cbind(nb.slingshot.sce$slingPseudotime_1, logcounts.hvg)
colnames(lin1.data)[1] <- "pseudotime"
lin1.df <- as.data.frame(as.matrix(lin1.data[!is.na(lin1.data[,1]), ]))

## Subset lineage 1 data into training and testing sets
## The model is fitted on the training set and evaluated on the testing set
require(tidymodels)
lin1.data.split <- initial_split(lin1.df)
lin1.data.train <- training(lin1.data.split)
lin1.data.test <- testing(lin1.data.split)

lin1.model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = lin1.data.train)

lin1.results <- lin1.data.test %>%
  mutate(estimate = predict(lin1.model, .[,-1]) %>% pull()) %>%
  select(truth = pseudotime, estimate)
metrics(data = lin1.results, truth, estimate)
## rmse  3.44 
## rsq  0.954
## mae  2.01 

summary(lin1.df$pseudotime)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.00   91.94   95.31   92.85   98.64  133.94

## plot genes suggested to be most important to predicting lineage 1 pseudotime
lin1.genes <- sort(lin1.model$fit$variable.importance, decreasing = TRUE)
lin1.top.genes <- names(lin1.genes)[1:6]

embedded <- embedCurves(nb.slingshot.sce, "UMAP")
lin1.embedded <- slingCurves(embedded)[[1]]
embedded.lin1 <- readRDS("slingshot_embedded_Lineage1.rds")

pal <- pals::ocean.dense(100)
par(mfrow = c(3,2))
for(i in seq_along(lin1.top.genes)){
  colours <- pal[cut(lin1.data[ ,lin1.top.genes[i]], breaks = 100)]
  plot(reducedDim(nb.slingshot.sce, "UMAP"), col = colours, pch = 16, cex = 0.5, main = lin1.top.genes[i])
  lines(embedded.lin1, lwd = 1, col = "black", type = "lineages")
}

## Prepare lineage 2 (QC) data for random forest
lin2.data <- cbind(nb.slingshot.sce$slingPseudotime_2, logcounts.hvg)
colnames(lin2.data)[1] <- "pseudotime"
lin2.df <- as.data.frame(as.matrix(lin2.data[!is.na(lin2.data[,1]), ]))

## Subset lineage 2 data into training and testing sets
lin2.data.split <- initial_split(lin2.df)
lin2.data.train <- training(lin2.data.split)
lin2.data.test <- testing(lin2.data.split)

lin2.model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = lin2.data.train)

lin2.results <- lin2.data.test %>%
  mutate(estimate = predict(lin2.model, .[,-1]) %>% pull()) %>%
  select(truth = pseudotime, estimate)
metrics(data = lin2.results, truth, estimate)
## rmse  3.45 
## rsq  0.959
## mae  2.20 

summary(lin2.df$pseudotime)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.00   92.58   95.97   94.27  101.83  171.48 

## Plot genes suggested to be most important to predicting lineage 2 pseudotime
lin2.genes <- sort(lin2.model$fit$variable.importance, decreasing = TRUE)
lin2.top.genes <- names(lin2.genes)[1:6]

embedded.lin2 <- readRDS("slingshot_embedded_Lineage2.rds")

par(mfrow = c(3,2))
for(i in seq_along(lin2.top.genes)){
  colours <- pal[cut(lin2.data[ ,lin2.top.genes[i]], breaks = 100)]
  plot(reducedDim(nb.slingshot.sce, "UMAP"), col = colours, pch = 16, cex = 0.5, main = lin2.top.genes[i])
  lines(embedded.lin2, lwd = 1, col = "black", type = "lineages")
}

## Prepare lineage 3 (cisplatin) data for random forest
lin3.data <- cbind(nb.slingshot.sce$slingPseudotime_3, logcounts.hvg)
colnames(lin3.data)[1] <- "pseudotime"
lin3.df <- as.data.frame(as.matrix(lin3.data[!is.na(lin3.data[,1]), ]))

## Subset lineage 3 data into training and testing sets
lin3.data.split <- initial_split(lin3.df)
lin3.data.train <- training(lin3.data.split)
lin3.data.test <- testing(lin3.data.split)

lin3.model <- rand_forest(mtry = 200, trees = 1400, min_n = 15, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = lin3.data.train)

lin3.results <- lin3.data.test %>%
  mutate(estimate = predict(lin3.model, .[,-1]) %>% pull()) %>%
  select(truth = pseudotime, estimate)
metrics(data = lin3.results, truth, estimate)
## rmse  3.22 
## rsq  0.959
## mae  1.98 

summary(lin3.df$pseudotime)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.00   92.92   96.22   94.15  101.99  151.16

## Plot genes suggested to be most important to predicting lineage 3 pseudotime
lin3.genes <- sort(lin3.model$fit$variable.importance, decreasing = TRUE)
lin3.top.genes <- names(lin3.genes)[1:6]

embedded.lin3 <- readRDS("slingshot_embedded_Lineage3.rds")

par(mfrow = c(3,2))
for(i in seq_along(lin3.top.genes)){
  colours <- pal[cut(lin3.data[ ,lin3.top.genes[i]], breaks = 100)]
  plot(reducedDim(nb.slingshot.sce, "UMAP"), col = colours, pch = 16, cex = 0.5, main = lin3.top.genes[i])
  lines(embedded.lin3, lwd = 1, col = "black", type = "lineages")
}
