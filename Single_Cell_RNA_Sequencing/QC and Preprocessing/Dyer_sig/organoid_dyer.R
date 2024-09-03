## NB organoids Dyer signature

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

## [ Organoid Dyer signature ] ----

organoid.seurat <- readRDS("data/organoid_seurat_start.rds")

dyer.df <- as.data.frame(read_xlsx(path = "data/Dyer_Sig.xlsx", col_names = FALSE))
colnames(dyer.df) <- c("Gene", "Signature")

dyer.sig <- lapply(unique(dyer.df$Signature), function(x){
  dyer.df[dyer.df$Signature==x, "Gene"]
})

names(dyer.sig) <- unique(dyer.df$Signature)
dyer.sig

organoid.seurat <- AddModuleScore(organoid.seurat, features = dyer.sig, assay = "SCT", seed = 12345, 
                                  name = c("MES.Dyer", "SYMP.Dyer","ADRN.Dyer"))
colnames(organoid.seurat@meta.data) <- gsub("\\.Dyer[1-9]$", "\\.Dyer", colnames(organoid.seurat@meta.data))

vln1.gg <- VlnPlot(organoid.seurat, 
                   features = c("MES.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln1.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln1.gg$layers[[2]] <- NULL

vln2.gg <- VlnPlot(organoid.seurat, 
                   features = c("SYMP.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln2.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln2.gg$layers[[2]] <- NULL

vln3.gg <- VlnPlot(organoid.seurat, 
                   features = c("ADRN.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln3.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln3.gg$layers[[2]] <- NULL

(vln2.gg | (vln3.gg / vln1.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/organoid_dyer_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/organoid_dyer_violins.png", width = 8.3, height = 8.3)

vln4.gg <- VlnPlot(organoid.seurat,
                   features = c("MES.Dyer"), group.by = "Sample_Type") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln4.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln4.gg$layers[[2]] <- NULL

vln5.gg <- VlnPlot(organoid.seurat,
                   features = c("SYMP.Dyer"), group.by = "Sample_Type") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln5.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln5.gg$layers[[2]] <- NULL

vln6.gg <- VlnPlot(organoid.seurat,
                   features = c("ADRN.Dyer"), group.by = "Sample_Type") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln6.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln6.gg$layers[[2]] <- NULL

(vln5.gg | (vln6.gg / vln4.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/organoid_dyer_condition_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/organoid_dyer_condition_violins.png", width = 8.3, height = 8.3)

organoid.seurat$Dyer.state <- ifelse(organoid.seurat$MES.Dyer > organoid.seurat$ADRN.Dyer &
                                       organoid.seurat$MES.Dyer > organoid.seurat$SYMP.Dyer, "MES",
                                     ifelse(organoid.seurat$ADRN.Dyer > organoid.seurat$MES.Dyer &
                                              organoid.seurat$ADRN.Dyer > organoid.seurat$SYMP.Dyer, "ADRN",
                                            ifelse(organoid.seurat$SYMP.Dyer > organoid.seurat$MES.Dyer &
                                                     organoid.seurat$SYMP.Dyer > organoid.seurat$ADRN.Dyer, "SYMP", NA)))

table(organoid.seurat$Dyer.state)
## ADRN   MES  SYMP 
## 10736  2863  7654 

table(organoid.seurat$Sample_Type, organoid.seurat$Dyer.state)
##                          ADRN  MES SYMP
## NB039_cisplatin          2434  696  354
## NB039_cisplatin_recovery 1149  369  215
## NB039_untreated          2691  142 4423
## NB067_cisplatin          1509  797  176
## NB067_cisplatin_recovery 2027  775  365
## NB067_untreated           926   84 2121

table(organoid.seurat$AMT.state)
## ADRN 
## 21253

table(organoid.seurat$Sample_Type, organoid.seurat$AMT.state)
## ADRN
## NB039_cisplatin          3484
## NB039_cisplatin_recovery 1733
## NB039_untreated          7256
## NB067_cisplatin          2482
## NB067_cisplatin_recovery 3167
## NB067_untreated          3131

saveRDS(organoid.seurat, "data/organoid_seurat_dyer.rds")

organoid.seurat <- readRDS("data/organoid_seurat_dyer.rds")

dyer.gg <- DimPlot(organoid.seurat,
                   group.by = "Dyer.state", order = TRUE) +
  umap.theme() + labs(title = "Dyer state") +
  scale_colour_manual(values = c("#990099", "#F37735", "#009999"),
                      labels = c("ADRN", "MES", "SYMP"))

dyer.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/organoid_umap_dyer_state.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/organoid_umap_dyer_state.png", width = 5.8, height = 5.8)

dyer.df <- as.data.frame(organoid.seurat@meta.data[c("Dyer.state", "Model")])

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

ggsave("plots/signatures/organoid_dyer_states_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/organoid_dyer_states_bar.png", width = 5.8, height = 5.8)

organoid.seurat$AMT.state <- factor(organoid.seurat$AMT.state, levels = c("ADRN", "intermediate", "MES"))

vln7.gg <- VlnPlot(organoid.seurat, 
                   features = c("SYMP.Dyer"), group.by = "AMT.state") +
  scale_fill_manual(values = c("#990099","lightgrey","#F37735"),
                    labels = c("ADRN", "Intermediate", "MES")) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln7.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln7.gg$layers[[2]] <- NULL

vln7.gg &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/organoid_vg_symp_violins.pdf", width = 8.3, height = 5.8)
ggsave("plots/signatures/organoid_vg_symp_violins.png", width = 8.3, height = 5.8)

## Bar charts per organoid model
dyer.df <- dyer.df[complete.cases(dyer.df), ]
table(dyer.df$Model, dyer.df$Dyer.state)
##       ADRN  MES SYMP
## NB039 6274 1207 4992
## NB067 4462 1656 2662

## NB039
nb039.df <- dyer.df[dyer.df$Model == "NB039", ]

nb039.gg <-
  nb039.df %>%
  ggplot(aes(Dyer.state, fill = Dyer.state)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values = c("#990099", "#F37735", "#009999")) +
  scale_x_discrete(labels = c("ADRN", "MES", "SYMP")) +
  labs(x = "",
       y = "Frequency",
       title = "NB039 Dyer states",
       fill = "Annotation") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

nb039.gg + theme(text = element_text(size = 15))

ggsave("plots/signatures/nb039_dyer_states_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/nb039_dyer_states_bar.png", width = 5.8, height = 5.8)

## NB067
nb067.df <- dyer.df[dyer.df$Model == "NB067", ]

nb067.gg <-
  nb067.df %>%
  ggplot(aes(Dyer.state, fill = Dyer.state)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values = c("#990099", "#F37735", "#009999")) +
  scale_x_discrete(labels = c("ADRN", "MES", "SYMP")) +
  labs(x = "",
       y = "Frequency",
       title = "NB067 Dyer states",
       fill = "Annotation") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

nb067.gg + theme(text = element_text(size = 15))

ggsave("plots/signatures/nb067_dyer_states_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/nb067_dyer_states_bar.png", width = 5.8, height = 5.8)
