## NB cell lines signatures

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

identity.cols <- c("EMT wound" = "#44AA99",
                   "P53 autophagy" = "#BB5566",
                   "EMT metabolism" = "#225555",
                   "Cycling MYC 1" = "#117733",
                   "Cycling Myc 2" = "#999933",
                   "Cycling EMT 1" = "#66CCEE",
                   "EMT DNA repair" = "#332288",
                   "Cycling EMT 2" = "#0077BB",
                   "Apoptosis OXPHOS" = "#882255",
                   "Cycing DNA repair" = "#EE6677")

condition.cols <- c("SK-N-SH_cisplatin" = "#DDAA33",
                    "SH-EP_untreated" = "#004488",
                    "SH-SY5Y_cisplatin_recovery" = "#A58a55",
                    "NA_NA" = "#555555",
                    "SH-SY5Y_cisplatin" = "#EEDD88", 
                    "SK-N-SH_untreated" = "#43ACFF",
                    "SH-EP_cisplatin_recovery" = "#758020",
                    "SH-SY5Y_untreated" = "#99DDFF",
                    "SK-N-SH_cisplatin_recovery" = "#9A6A20",
                    "SH-EP_cisplatin" = "#BBCC33")

## [ Cell lines signature plots ] ----

lines.seurat <- readRDS("data/nb_seurat_AMT_non-batch_cell_lines.rds")

## Cluster colours
#cluster.cols <- iwanthue(length(unique(lines.seurat$seurat_clusters)))
#saveRDS(cluster.cols, "data/cell_lines_cluster_cols.rds")
cluster.cols <- readRDS("data/cell_lines_cluster_cols.rds")

num.gg <- DimPlot(lines.seurat,
                  group.by = "seurat_clusters", order = TRUE) +
  umap.theme() + labs(title = "Cluster number") +
  scale_colour_manual(values = cluster.cols)

num.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/cell_line_umap_cluster_number.pdf", width = 5.8, height = 5.8)
ggsave("plots/cell_line_umap_cluster_number.png", width = 5.8, height = 5.8)

dir.create("plots/signatures/", recursive = TRUE)

adrn.gg <- FeaturePlot(lines.seurat,
                       features = "ADRN.Sig",
                       order = TRUE,
                       min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.matter(100)) +
  labs(title = "vanGroningen\nADRN signature", colour = "Expression") +
  umap.theme()

mes.gg <- FeaturePlot(lines.seurat,
                      features = "MES.Sig",
                      order = TRUE,
                      min.cutoff = 0, max.cutoff = 0.5) +
  scale_colour_gradientn(colours = pals::ocean.haline(100)) +
  labs(title = "vanGroningen\nMES signature", colour = "Expression") +
  umap.theme()

adrn.gg + mes.gg &
  theme(text = element_text(size = 15),
        legend.key.size = unit(1, "cm"))

ggsave("plots/signatures/cell_line_vG_ADRN_MES_umaps.pdf", width = 5.8, height = 8.3)
ggsave("plots/signatures/cell_line_vG_ADRN_MES_umaps.png", width = 5.8, height = 8.3)

vln1.gg <- VlnPlot(lines.seurat,
                   features = "MES.Sig", group.by = "seurat_clusters") +
  scale_fill_manual(values = cluster.cols) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln1.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln1.gg$layers[[2]] <- NULL

vln2.gg <- VlnPlot(lines.seurat,
                   features = "ADRN.Sig", group.by = "seurat_clusters") +
  scale_fill_manual(values = cluster.cols) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln2.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln2.gg$layers[[2]] <- NULL

vln3.gg <- VlnPlot(lines.seurat,
                   features = "AMT.score", group.by = "seurat_clusters") +
  scale_fill_manual(values = cluster.cols) &
  umap.theme() + theme(aspect.ratio = 0.5)

vln3.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln3.gg$layers[[2]] <- NULL

(vln3.gg | (vln1.gg / vln2.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/cell_line_vG_ADRN_MES_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/cell_line_vG_ADRN_MES_violins.png", width = 8.3, height = 8.3)

amt.gg <- DimPlot(lines.seurat,
                  group.by = "AMT.state", order = TRUE) +
  umap.theme() + labs(title = "AMT state") +
  scale_colour_manual(values = c("#990099", "#F37735", "lightgrey"),
                      labels = c("ADRN", "MES", "Intermediate"))

amt.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/cell_line_umap_amt_state.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/cell_line_umap_amt_state.png", width = 5.8, height = 5.8)

lines.seurat$Cluster_Name <- recode(as.character(lines.seurat$seurat_clusters),
                                    "0" = "EMT wound",
                                    "1" = "P53 autophagy",
                                    "2" = "EMT metabolism",
                                    "3" = "Cycling MYC 1",
                                    "4" = "Cycling Myc 2",
                                    "5" = "Cycling EMT 1",
                                    "6" = "EMT DNA repair",
                                    "7" = "Cycling EMT 2",
                                    "8" = "Apoptosis OXPHOS",
                                    "9" = "Cycing DNA repair")

cluster.gg <- DimPlot(lines.seurat,
                      group.by = "Cluster_Name", order = TRUE) +
  umap.theme() + labs(title = "Initial Cluster Annotation") +
  scale_colour_manual(values = identity.cols)

cluster.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/cell_line_umap_initial_clusters.pdf", width = 5.8, height = 5.8)
ggsave("plots/cell_line_umap_initial_clusters.png", width = 5.8, height = 5.8)

## [ Cell lines coexpression ] ----

## Which cells coexpress van Groningen ADRN and MES signatures
coexprs = function(x){
  if (x[["MES.Sig"]] > 0.2 & x[["ADRN.Sig"]] > 0.1) {
    result <- 1
  } else {
    result <- 0
  }
  return(result)
}

metadata.df <- subset(lines.seurat@meta.data, select = c("MES.Sig", "ADRN.Sig"))

coexprs.df <- apply(metadata.df, 1, coexprs)
lines.seurat@meta.data <- cbind(lines.seurat@meta.data, Coexprs.Sig = coexprs.df)

dir.create("plots/coexpression", recursive = TRUE)

coexprs.gg <- DimPlot(lines.seurat,
                      group.by = "Coexprs.Sig", order = TRUE) +
  umap.theme() + labs(title = "ADRN and MES Coexpression\n(MES > 0.2, ADRN > 0.1)") +
  scale_colour_manual(values = c("#4477AA", "#EE6677"), labels = c("FALSE", "TRUE")) + 
  labs(colour = "Coexpression")

coexprs.gg

ggsave("plots/coexpression/cell_line_umap_coexpress.pdf", width = 5.8, height = 5.8)
ggsave("plots/coexpression/cell_line_umap_coexpress.png", width = 5.8, height = 5.8)

saveRDS(lines.seurat, "data/cell_lines_seurat_start.rds")

## [ Cell lines Dyer signature ] ----

lines.seurat <- readRDS("data/cell_lines_seurat_start.rds")

dyer.df <- as.data.frame(read_xlsx(path = "data/Dyer_Sig.xlsx", col_names = FALSE))
colnames(dyer.df) <- c("Gene", "Signature")

dyer.sig <- lapply(unique(dyer.df$Signature), function(x){
  dyer.df[dyer.df$Signature==x, "Gene"]
})

names(dyer.sig) <- unique(dyer.df$Signature)
dyer.sig

lines.seurat <- AddModuleScore(lines.seurat, features = dyer.sig, assay = "SCT", seed = 12345, 
                               name = c("MES.Dyer", "SYMP.Dyer","ADRN.Dyer"))
colnames(lines.seurat@meta.data) <- gsub("\\.Dyer[1-9]$", "\\.Dyer", colnames(lines.seurat@meta.data))

vln1.gg <- VlnPlot(lines.seurat, 
                   features = c("MES.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln1.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln1.gg$layers[[2]] <- NULL

vln2.gg <- VlnPlot(lines.seurat, 
                   features = c("SYMP.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln2.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln2.gg$layers[[2]] <- NULL

vln3.gg <- VlnPlot(lines.seurat, 
                   features = c("ADRN.Dyer"), group.by = "Cluster_Name") +
  scale_fill_manual(values = identity.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln3.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln3.gg$layers[[2]] <- NULL

(vln2.gg | (vln3.gg / vln1.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/cell_line_dyer_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/cell_line_dyer_violins.png", width = 8.3, height = 8.3)

lines.seurat$Sample_Type <- factor(lines.seurat$Sample_Type, c("SK-N-SH_untreated", "SK-N-SH_cisplatin", "SK-N-SH_cisplatin_recovery",
                                                               "SH-EP_untreated", "SH-EP_cisplatin", "SH-EP_cisplatin_recovery",
                                                               "SH-SY5Y_untreated", "SH-SY5Y_cisplatin", "SH-SY5Y_cisplatin_recovery",
                                                               "NA_NA"))

vln4.gg <- VlnPlot(lines.seurat, 
                   features = c("MES.Dyer"), group.by = "Sample_Type") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln4.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln4.gg$layers[[2]] <- NULL

vln5.gg <- VlnPlot(lines.seurat, 
                   features = c("SYMP.Dyer"), group.by = "Sample_Type") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln5.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln5.gg$layers[[2]] <- NULL

vln6.gg <- VlnPlot(lines.seurat, 
                   features = c("ADRN.Dyer"), group.by = "Sample_Type") +
  scale_fill_manual(values = condition.cols) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln6.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln6.gg$layers[[2]] <- NULL

(vln5.gg | (vln6.gg / vln4.gg)) + patchwork::plot_layout(guides = 'collect') &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/cell_line_dyer_condition_violins.pdf", width = 8.3, height = 8.3)
ggsave("plots/signatures/cell_line_dyer_condition_violins.png", width = 8.3, height = 8.3)

lines.seurat$Dyer.state <- ifelse(lines.seurat$MES.Dyer > lines.seurat$ADRN.Dyer &
                                    lines.seurat$MES.Dyer > lines.seurat$SYMP.Dyer, "MES",
                                  ifelse(lines.seurat$ADRN.Dyer > lines.seurat$MES.Dyer &
                                           lines.seurat$ADRN.Dyer > lines.seurat$SYMP.Dyer, "ADRN",
                                         ifelse(lines.seurat$SYMP.Dyer > lines.seurat$MES.Dyer &
                                                  lines.seurat$SYMP.Dyer > lines.seurat$ADRN.Dyer, "SYMP", NA)))

table(lines.seurat$Dyer.state)
## ADRN  MES SYMP 
## 7542 8726 7020  

table(lines.seurat$Sample_Type, lines.seurat$Dyer.state)
##                             ADRN  MES SYMP
## SK-N-SH_untreated           730  830 1131
## SK-N-SH_cisplatin          1095  984  587
## SK-N-SH_cisplatin_recovery  142  263  402
## SH-EP_untreated             232 1961  237
## SH-EP_cisplatin              65  328   35
## SH-EP_cisplatin_recovery    133 1044  293
## SH-SY5Y_untreated           620  263 1291
## SH-SY5Y_cisplatin           762  333  280
## SH-SY5Y_cisplatin_recovery 1827  174 1569
## NA_NA                      1936 2546 1195

table(lines.seurat$AMT.state)
## ADRN          MES intermediate 
## 4361        15496         3431 

table(lines.seurat$Sample_Type, lines.seurat$AMT.state)
##                            ADRN  MES intermediate
## SK-N-SH_untreated           206 2396           89
## SK-N-SH_cisplatin            97 1999          570
## SK-N-SH_cisplatin_recovery  198  569           40
## SH-EP_untreated              20 2360           50
## SH-EP_cisplatin               0  425            3
## SH-EP_cisplatin_recovery      0 1468            2
## SH-SY5Y_untreated           269 1050          855
## SH-SY5Y_cisplatin            84  769          522
## SH-SY5Y_cisplatin_recovery 3018   89          463
## NA_NA                       469 4371          837

saveRDS(lines.seurat, "data/cell_lines_seurat_dyer.rds")

lines.seurat <- readRDS("data/cell_lines_seurat_dyer.rds")

dyer.gg <- DimPlot(lines.seurat,
                   group.by = "Dyer.state", order = TRUE) +
  umap.theme() + labs(title = "Dyer state") +
  scale_colour_manual(values = c("#990099", "#F37735", "#009999"),
                      labels = c("ADRN", "MES", "SYMP"))

dyer.gg +
  guides(colour = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/cell_line_umap_dyer_state.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/cell_line_umap_dyer_state.png", width = 5.8, height = 5.8)

dyer.df <- as.data.frame(lines.seurat@meta.data[c("Dyer.state", "Model")])

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

ggsave("plots/signatures/cell_line_dyer_states_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/cell_line_dyer_states_bar.png", width = 5.8, height = 5.8)

lines.seurat$AMT.state <- factor(lines.seurat$AMT.state, levels = c("ADRN", "intermediate", "MES"))

vln7.gg <- VlnPlot(lines.seurat, 
                   features = c("SYMP.Dyer"), group.by = "AMT.state") +
  scale_fill_manual(values = c("#990099","lightgrey","#F37735"),
                    labels = c("ADRN", "Intermediate", "MES")) &
  umap.theme() + theme(aspect.ratio = 0.5,
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

vln7.gg$layers[[1]] <- geom_quasirandom(shape = 21, alpha = 1, stroke = 0.1)
vln7.gg$layers[[2]] <- NULL

vln7.gg &
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/signatures/cell_line_vg_symp_violins.pdf", width = 8.3, height = 5.8)
ggsave("plots/signatures/cell_line_vg_symp_violins.png", width = 8.3, height = 5.8)

## Bar charts per cell line model
dyer.df <- dyer.df[complete.cases(dyer.df), ]
table(dyer.df$Model, dyer.df$Dyer.state)
##         ADRN  MES SYMP
## SH-EP    430 3333  565
## SH-SY5Y 3209  770 3140
## SK-N-SH 1967 2077 2120

## SK-N-SH
sknsh.df <- dyer.df[dyer.df$Model == "SK-N-SH", ]

sknsh.gg <-
  sknsh.df %>%
  ggplot(aes(Dyer.state, fill = Dyer.state)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values = c("#990099", "#F37735", "#009999")) +
  scale_x_discrete(labels = c("ADRN", "MES", "SYMP")) +
  labs(x = "",
       y = "Frequency",
       title = "SK-N-SH Dyer states",
       fill = "Annotation") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

sknsh.gg + theme(text = element_text(size = 15))

ggsave("plots/signatures/sknsh_dyer_states_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/sknsh_dyer_states_bar.png", width = 5.8, height = 5.8)

## SH-EP
shep.df <- dyer.df[dyer.df$Model == "SH-EP", ]

shep.gg <-
  shep.df %>%
  ggplot(aes(Dyer.state, fill = Dyer.state)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values = c("#990099", "#F37735", "#009999")) +
  scale_x_discrete(labels = c("ADRN", "MES", "SYMP")) +
  labs(x = "",
       y = "Frequency",
       title = "SH-EP Dyer states",
       fill = "Annotation") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

shep.gg + theme(text = element_text(size = 15))

ggsave("plots/signatures/shep_dyer_states_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/shep_dyer_states_bar.png", width = 5.8, height = 5.8)

## SH-SY5Y
shsy5y.df <- dyer.df[dyer.df$Model == "SH-SY5Y", ]

shsy5y.gg <-
  shsy5y.df %>%
  ggplot(aes(Dyer.state, fill = Dyer.state)) +
  geom_bar() +
  theme_bw() +
  scale_fill_manual(values = c("#990099", "#F37735", "#009999")) +
  scale_x_discrete(labels = c("ADRN", "MES", "SYMP")) +
  labs(x = "",
       y = "Frequency",
       title = "SH-SY5Y Dyer states",
       fill = "Annotation") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

shsy5y.gg + theme(text = element_text(size = 15))

ggsave("plots/signatures/shsy5y_dyer_states_bar.pdf", width = 5.8, height = 5.8)
ggsave("plots/signatures/shsy5y_dyer_states_bar.png", width = 5.8, height = 5.8)
