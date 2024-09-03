## NB PDX Dyer signature

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)
library(tidyverse)
library(ggbeeswarm)
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

## [ PDX Dyer signature ] ----

pdx.seurat <- readRDS("data/pdx_seurat_start.rds")

dyer.df <- as.data.frame(read_xlsx(path = "data/Dyer_Sig.xlsx", col_names = FALSE))
colnames(dyer.df) <- c("Gene", "Signature")

dyer.sig <- lapply(unique(dyer.df$Signature), function(x){
  dyer.df[dyer.df$Signature==x, "Gene"]
})

names(dyer.sig) <- unique(dyer.df$Signature)
dyer.sig

pdx.seurat <- AddModuleScore(pdx.seurat, features = dyer.sig, assay = "SCT", seed = 12345, 
                             name = c("MES.Dyer", "SYMP.Dyer","ADRN.Dyer"))
colnames(pdx.seurat@meta.data) <- gsub("\\.Dyer[1-9]$", "\\.Dyer", colnames(pdx.seurat@meta.data))

vln1.gg <- VlnPlot(pdx.seurat, 
                   features = c("MES.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln1.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln1.gg$layers[[2]] <- NULL

vln2.gg <- VlnPlot(pdx.seurat, 
                   features = c("SYMP.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln2.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln2.gg$layers[[2]] <- NULL

vln3.gg <- VlnPlot(pdx.seurat, 
                   features = c("ADRN.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln3.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln3.gg$layers[[2]] <- NULL

(vln2.gg | (vln3.gg / vln1.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/pdx_dyer_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/pdx_dyer_violins.png", width = 8.3, height = 8.3)

condition.cols <- c("GRNB5_cisplatin" = "#EE8866",
                    "GRNB5_cisplatin_recovery" = "#44BB99",
                    "GRNB5_untreated" = "#77AADD")

vln4.gg <- VlnPlot(pdx.seurat, 
                   features = c("MES.Dyer"), group.by = "Sample_Type") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln4.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln4.gg$layers[[2]] <- NULL

vln5.gg <- VlnPlot(pdx.seurat, 
                   features = c("SYMP.Dyer"), group.by = "Sample_Type") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln5.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln5.gg$layers[[2]] <- NULL

vln6.gg <- VlnPlot(pdx.seurat, 
                   features = c("ADRN.Dyer"), group.by = "Sample_Type") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln6.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln6.gg$layers[[2]] <- NULL

(vln5.gg | (vln6.gg / vln4.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/pdx_dyer_condition_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/pdx_dyer_condition_violins.png", width = 8.3, height = 8.3)

pdx.seurat$Dyer.state <- ifelse(pdx.seurat$MES.Dyer > pdx.seurat$ADRN.Dyer &
                                  pdx.seurat$MES.Dyer > pdx.seurat$SYMP.Dyer, "MES",
                                ifelse(pdx.seurat$ADRN.Dyer > pdx.seurat$MES.Dyer &
                                         pdx.seurat$ADRN.Dyer > pdx.seurat$SYMP.Dyer, "ADRN",
                                       ifelse(pdx.seurat$SYMP.Dyer > pdx.seurat$MES.Dyer &
                                                pdx.seurat$SYMP.Dyer > pdx.seurat$ADRN.Dyer, "SYMP", NA)))

table(pdx.seurat$Dyer.state)
## ADRN   MES  SYMP 
## 18444 15258 11938 

table(pdx.seurat$Condition, pdx.seurat$Dyer.state)
##                    ADRN  MES SYMP
## cisplatin          6169 7279 3626
## cisplatin recovery 3577 4687 1587
## untreated          8698 3292 6725

table(pdx.seurat$AMT.state)
## ADRN          MES intermediate 
## 5311        38875         1454

table(pdx.seurat$Condition, pdx.seurat$AMT.state)
##                      ADRN  MES intermediate
## cisplatin            486 16411          177
## cisplatin recovery   157  9601           93
## untreated           4668 12863         1184

saveRDS(pdx.seurat, "data/pdx_seurat_dyer.rds")

pdx.seurat <- readRDS("data/pdx_seurat_dyer.rds")

dyer.gg <- DimPlot(pdx.seurat,
                   group.by = "Dyer.state", order = TRUE) +
  umap.theme() + labs(title = "Dyer state") +
  scale_colour_manual(values = c("#990099", "#F37735", "#009999"),
                      labels = c("ADRN", "MES", "SYMP"))

dyer.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/pdx_umap_dyer_state.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/pdx_umap_dyer_state.png", width = 5.8, height = 5.8)

dyer.df <- as.data.frame(pdx.seurat@meta.data["Dyer.state"])

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

ggsave("plots/signatures/pdx_dyer_states_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/pdx_dyer_states_bar.png", width = 5.8, height = 5.8)

pdx.seurat$AMT.state <- factor(pdx.seurat$AMT.state, levels = c("ADRN", "intermediate", "MES"))

vln7.gg <- VlnPlot(pdx.seurat, 
                   features = c("SYMP.Dyer"), group.by = "AMT.state") +
  scale_fill_manual(values = c("#990099","lightgrey","#F37735"),
                    labels = c("ADRN", "Intermediate", "MES")) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln7.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln7.gg$layers[[2]] <- NULL

vln7.gg &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/pdx_vg_symp_violins.pdf", width = 8.3, height = 5.8)
ggsave("plots/signatures/pdx_vg_symp_violins.png", width = 8.3, height = 5.8)