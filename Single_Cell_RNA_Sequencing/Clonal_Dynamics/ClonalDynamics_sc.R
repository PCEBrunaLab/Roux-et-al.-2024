#ClonalDynamics_SC.R

required.packages <- c("scater", "scran", "scuttle", "Seurat", "patchwork", 
                       "biomaRt", "reshape2", "EnsDb.Hsapiens.v86", "ggsankey", 
                       "tricycle", "hues", "dplyr", "SingleCellExperiment", "ggsci",
                       "ggthemes", "cowplot", "scattermore", "viridis", "MatrixGenerics",
                       "scale", "enrichR", "umap", "rstatix", "ggpubr", "Polychrome", 
                       "ggalluvial", "forcats", "vegan", "ggVennDiagram", "purrr")

if (!all(required.packages %in% installed.packages()[,"Package"])){
  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]])
}

library(dplyr)
library(SingleCellExperiment)
library(tidyverse)
library(ggsci)
library(ggthemes)
library(cowplot)
library(reshape2)
library(scater)
library(scattermore)
library(viridis)
library(biomaRt)
library(MatrixGenerics)
library(scales)
library(enrichR)
library(scran)
library(umap)
library(rstatix)
library(ggpubr)
library(Polychrome)
library(ggalluvial)
library(forcats)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(hues)
library(patchwork)
library(vegan)
library(ggVennDiagram)
library(purrr)

setwd("~/OneDrive - The Institute of Cancer Research/Sian/Data/2021 scRNA-seq/Preprocessing and Downstream Analysis/Sian/New Analysis/Barcode Analysis/")

#Load in sce and metadata
big.sce <- readRDS("datafiles/nb_sce_AMT.rds")
big.sce <- big.sce[!grepl(rownames(big.sce), pattern="CMO[0-9]+"),]

amt.anno<-as.data.frame(colData(big.sce)[,c("AMT.score","AMT.state")])
amt.anno$CellID<-rownames(amt.anno)

meta.df <- read.table("datafiles/Neuro_noS_mergedMeta.txt",
                      sep="\t", header=TRUE, stringsAsFactors=FALSE)
meta.df <- meta.df %>%
  select(-c(UMAP1, UMAP2, Graph.Cluster, ClusterLabel, Cluster.Sample.Condition))

bcs.df <- read.table("datafiles/All_Neuro_metaWbarcodes.txt",
                     sep="\t", header=TRUE, stringsAsFactors=FALSE)
bcs.df <- bcs.df[bcs.df$Keep.BC == 1, ]
rownames(bcs.df) <- bcs.df$CellID

#Extract UMAP information
meta.extract.df <- as.data.frame(colData(big.sce))[, c( "sample_id", "description", "seurat_clusters.0.2", "Cluster_Name", "AMT.score")]
#saveRDS(meta.extract.df, "meta_df.rds")

umap.df <- as.data.frame(reducedDim(big.sce, "UMAP"))
#saveRDS(umap.df, "umap_df.rds")

#Introduce new cluster and umap information
meta.extract.df$CellID <- rownames(meta.extract.df)
umap.df$CellID <- rownames(umap.df)

# Join all data frames in list
list_df <- list(meta.extract.df, umap.df, meta.df) 
meta.merge.df <- Reduce(function(x, y) full_join(x, y), list_df)

# Add barcode information to cells
colData(big.sce)$Full.BCS <- bcs.df[colnames(big.sce), ]$Full.BCS
full.meta.df <- merge(meta.merge.df, bcs.df[, c("CellID", "BC.30", "BC.14", "Full.BCS", "Count")], by='CellID', all.x=TRUE)
n.clust <- length(unique(full.meta.df$seurat_clusters.0.2[!(is.na(full.meta.df$seurat_clusters.0.2))]))
clust.cols <- pal_igv()(n.clust)
names(clust.cols) <- levels(full.meta.df$seurat_clusters.0.2)

full.meta.df$seurat_clusters.0.2 <- factor(full.meta.df$seurat_clusters.0.2,
                                           levels=c(0:n.clust))
sub.sce <- big.sce[, !big.sce$Batch %in% c("Neuro_S")]

#Remove any cells with no Cellecta barcodes
sub.sce <- sub.sce[, !is.na(colData(sub.sce)$Full.BCS)]
full.meta.df <- full.meta.df[full.meta.df$CellID %in% colnames(sub.sce), ]

#Change name of barcodes variable
names(full.meta.df)[names(full.meta.df) == 'Full.BCS'] <- 'real_bc44'
full.meta.df <- merge(full.meta.df, amt.anno, by=c("CellID", "AMT.score"))

#Names have changed since last time so lets change the names
full.meta.df$Cluster_Name[full.meta.df$Cluster_Name == "Large"] <- "EMT"
full.meta.df$Cluster_Name[full.meta.df$Cluster_Name == "JQ1"] <- "Chromatin Remodelling"
full.meta.df$Cluster_Name[full.meta.df$Cluster_Name == "QC"] <- "Stress Response"
full.meta.df$Cluster_Name[full.meta.df$Cluster_Name == "Cisplatin"] <- "Apoptotic"
full.meta.df$Cluster_Name[full.meta.df$Cluster_Name == "GrpA"] <- "Cycling MES"
full.meta.df$Cluster_Name[full.meta.df$Cluster_Name == "GrpB"] <- "Cycling ADRN"

#Remove unnessecary files
rm(big.sce, umap.df, meta.df, meta.extract.df, amt.anno, 
   bcs.df, list_df, meta.merge.df, n.clust, sub.sce)
gc()

nb.seurat <- readRDS("datafiles/nb_seurat_AMT.rds")

#Generate colour palette for seurat_clusters_0.2
P9 <- c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA")
P9 <- as.vector(t(matrix(P9)))
names(P9) = c("0", "1", "2", "3", "4", "5", "6", "7", "8")

#Define colour palette for conditions
evalVignette <- requireNamespace("ggplot2", quietly = TRUE)
knitr::opts_chunk$set(eval = evalVignette)
P6 = createPalette(6,  c("#ff0000", "#00ff00", "#0000ff"))
P6 <- as.vector(t(matrix(P6)))
names(P6) = unique(nb.seurat$Condition)

group.cols <- c("Untreated_rec" = "#F7000D", 
                "JQ1" = "#22FF0D", 
                "Untreated" = "#1C16FE", 
                "Cisplatin" = "#FFD1D3", 
                "Cisplatin_rec" = "#FE16E1", 
                "JQ1_rec" = "#00D1FD")

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

library(harmony)
library(ggsankey)
library(tidyverse)

# scRNA-seq barcodes - Initial Processing
#Plot Manuscript Figure 3a
DimPlot(nb.seurat,
                   group.by = "seurat_clusters.0.2", order = TRUE) +
  scale_colour_manual(values = P9) +
  umap.theme() + labs(title = "10x Clusters res=0.2")

#Plot Manuscript Figure 3c
adrn.mes.gg <- DimPlot(nb.seurat,
                       group.by = "AMT.state") +
  scale_colour_manual(breaks = c("MES", "intermediate", "ADRN"),
                      values=c("#F37735", "lightgrey", "#990099")) +
  umap.theme() +
  labs(title = "AMT state")

#Cellular cluster plots
#Plot Manuscript Figure 3b
full.meta.df %>%
  mutate(Condition = factor(Condition, levels=c("Untreated", "Untreated_rec", "Cisplatin", "Cisplatin_rec", "JQ1", "JQ1_rec"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(aes(colour = seurat_clusters.0.2), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=8)) +
  scale_colour_manual(values=P9) +
  facet_wrap(~Condition, nrow = 2, ncol = 3) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))


#AMT state plots
#Plot Manuscript Figure 3d
full.meta.df %>%
  mutate(Condition = factor(Condition, levels=c("Untreated", "Untreated_rec", "Cisplatin", "Cisplatin_rec", "JQ1", "JQ1_rec"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(aes(colour = AMT.state), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(breaks = c("MES", "intermediate", "ADRN"),
                      values=c("#F37735", "lightgrey", "#990099")) +
  facet_wrap(~Condition, nrow = 2, ncol = 3) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

tmp <- full.meta.df %>%
  group_by(AMT.state, Condition, Sample) %>%
  mutate(count=n()) %>%
  ungroup() %>%
  group_by(Condition, Sample) %>%
  mutate(total = sum(count)) %>%
  ungroup() %>%
  mutate(percentage = (count/total)*100) %>%
  ungroup()

filtered.2.meta.df <- tmp %>%
  group_by(Condition, AMT.state, Sample) %>%
  summarise(total_perc = sum(percentage)) %>%
  ungroup() %>%
  mutate(Condition = fct_relevel(Condition, "Untreated", "Untreated_rec",
                                 "Cisplatin", "Cisplatin_rec", "JQ1", "JQ1_rec"))

amt.perc.change <- filtered.2.meta.df %>%
  group_by(AMT.state, Sample) %>%
  mutate(percentage_change = total_perc/total_perc[Condition == "Untreated"]) %>%
  ungroup()

library(stats)
amt.perc.change.summary <- amt.perc.change %>%
  group_by(Condition, AMT.state) %>%
  mutate(log_perc_change = log10(percentage_change),
         mean = mean(log_perc_change),
         sd = sd(log_perc_change))

test <- merge(amt.perc.change, amt.perc.change.summary)

#Plot Manuscript Figure 3e
test %>%
  filter(Condition != "Untreated" & Condition != "Untreated_rec") %>%
  group_by(Condition, Sample) %>%
  ggplot( aes(x=Condition, y=log_perc_change, fill = AMT.state)) +
  geom_bar(stat="summary", position="dodge") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  facet_wrap(~AMT.state, scales = "free") +
  ylab("Percentage Change \nfrom Untreated") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        aspect.ratio = 1)


#Conduct pairwise t test between populations
#Create df per amt state
amt.perc.change.mes <- amt.perc.change %>%
  ungroup() %>%
  filter(AMT.state == "MES") %>%
  select(-c(AMT.state, Sample, total_perc))

amt.perc.change.inter <- amt.perc.change %>%
  ungroup() %>%
  filter(AMT.state == "intermediate") %>%
  select(-c(AMT.state, Sample, total_perc))

amt.perc.change.adrn <- amt.perc.change %>%
  ungroup() %>%
  filter(AMT.state == "ADRN") %>%
  select(-c(AMT.state, Sample, total_perc))

n_sample <- unique(amt.perc.change$Condition)
first = T
for (i in 1:length(n_sample)) {
  for(j in 1:length(n_sample)){
    if(i > j){  
      
      tmp = t.test(amt.perc.change.mes$percentage_change[amt.perc.change.mes$Condition==n_sample[i]],
                   amt.perc.change.mes$percentage_change[amt.perc.change.mes$Condition==n_sample[j]])
      
      if(first){
        df_pval = data.frame(GroupA = n_sample[i],
                             GroupB = n_sample[j],
                             pvalue = tmp$p.value)
        first = F
      } else {
        df_pval = rbind.data.frame(df_pval,
                                   data.frame(GroupA = n_sample[i],
                                              GroupB = n_sample[j],
                                              pvalue = tmp$p.value))
      }
    }
  }
}

df_pval$fdr = p.adjust(df_pval$pvalue, method = "fdr", n = length(df_pval$pvalue))
df_pval <- df_pval %>%
  mutate(signif = ifelse(df_pval$fdr <= 0.05, "TRUE", "FALSE"),
         signif_value = ifelse(df_pval$fdr <= 0.001, "***", 
                               ifelse(df_pval$fdr <= 0.01, "**",
                                      ifelse(df_pval$fdr <= 0.05, "*", "NA"))))

#Statstical test for Figure 3e - MES
df_pval_mes <- df_pval

first = T
for (i in 1:length(n_sample)) {
  for(j in 1:length(n_sample)){
    if(i > j){  
      
      tmp = t.test(amt.perc.change.inter$percentage_change[amt.perc.change.inter$Condition==n_sample[i]],
                   amt.perc.change.inter$percentage_change[amt.perc.change.inter$Condition==n_sample[j]])
      
      if(first){
        df_pval = data.frame(GroupA = n_sample[i],
                             GroupB = n_sample[j],
                             pvalue = tmp$p.value)
        first = F
      } else {
        df_pval = rbind.data.frame(df_pval,
                                   data.frame(GroupA = n_sample[i],
                                              GroupB = n_sample[j],
                                              pvalue = tmp$p.value))
      }
    }
  }
}

df_pval$fdr = p.adjust(df_pval$pvalue, method = "fdr", n = length(df_pval$pvalue))
df_pval <- df_pval %>%
  mutate(signif = ifelse(df_pval$fdr <= 0.05, "TRUE", "FALSE"),
         signif_value = ifelse(df_pval$fdr <= 0.001, "***", 
                               ifelse(df_pval$fdr <= 0.01, "**",
                                      ifelse(df_pval$fdr <= 0.05, "*", "NA"))))

#Statstical test for Figure 3e - HYB
df_pval_inter <- df_pval

first = T
for (i in 1:length(n_sample)) {
  for(j in 1:length(n_sample)){
    if(i > j){  
      
      tmp = t.test(amt.perc.change.adrn$percentage_change[amt.perc.change.adrn$Condition==n_sample[i]],
                   amt.perc.change.adrn$percentage_change[amt.perc.change.adrn$Condition==n_sample[j]])
      
      if(first){
        df_pval = data.frame(GroupA = n_sample[i],
                             GroupB = n_sample[j],
                             pvalue = tmp$p.value)
        first = F
      } else {
        df_pval = rbind.data.frame(df_pval,
                                   data.frame(GroupA = n_sample[i],
                                              GroupB = n_sample[j],
                                              pvalue = tmp$p.value))
      }
    }
  }
}

df_pval$fdr = p.adjust(df_pval$pvalue, method = "fdr", n = length(df_pval$pvalue))
df_pval <- df_pval %>%
  mutate(signif = ifelse(df_pval$fdr <= 0.05, "TRUE", "FALSE"),
         signif_value = ifelse(df_pval$fdr <= 0.001, "***", 
                               ifelse(df_pval$fdr <= 0.01, "**",
                                      ifelse(df_pval$fdr <= 0.05, "*", "NA"))))

#Statstical test for Figure 3e - ADRN
df_pval_adrn <- df_pval


# scRNA-seq barcodes - Alluvial Untreated Initial
ut.cellular.dtp.df <- full.meta.df %>%
  filter(Condition == "Untreated" | Condition == "Untreated_rec") %>%
  select(c(seurat_clusters.0.2, Sample, real_bc44, Condition, Count))

ut.cellular.dtp.df <- ut.cellular.dtp.df[, c(3, 4, 1, 2, 5)]

tmp <- ut.cellular.dtp.df %>%
  group_by(Sample, Condition, real_bc44) %>%
  mutate(total_clones=sum(Count)) %>%
  ungroup() %>%
  group_by(real_bc44, Condition, Sample) %>%
  mutate(percentage = Count/total_clones * 100) %>%
  ungroup() 

tmp$timepoint[tmp$Condition == "Untreated"] <- "T1"
tmp$timepoint[tmp$Condition == "Untreated_rec"] <- "T2"

tmp2 <- tmp %>%
  group_by(real_bc44, Condition, Sample) %>%
  summarise(count = n_distinct(seurat_clusters.0.2))

ut.cellular.states <- tmp2 %>%
  group_by(count, Condition, Sample) %>%
  summarise(clones = length(unique(real_bc44)))

ut.cellular.states.summary <- ut.cellular.states %>%
  filter(Condition == "Untreated_rec")

#Calculate percentage of clones with each representation
ut.cellular.states.perc <- ut.cellular.states.summary %>%
  group_by(Sample, Condition) %>%
  mutate(total_freq = sum(clones)) %>%
  ungroup() %>%
  mutate(perc = clones/total_freq * 100) %>%
  ungroup()

set.seed(12)
P8 = createPalette(8, c("#ff0000", "#00ff00", "#0000ff"), M=10)
P8 <- as.vector(t(matrix(P8)))
names(P8) = unique(as.character(ut.cellular.states.perc$count))
ut.cellular.states.perc$count <- as.character(ut.cellular.states.perc$count)

#Plot manuscript exteneded figure 7d - UT
ut.cellular.states.perc %>%
  ggplot( ,aes(x=Sample, y=perc))+
  geom_bar(aes(x=Sample, y=perc, fill=count), colour="black", stat = "identity", position = "stack") +
  xlab("States Occupied")+
  ylab("Number of Clones")+
  ggtitle("UT Recovery Cellular Cluster Clone Summary")+
  scale_fill_manual(values = P8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))
ggsave("manuscript_plots/ut_cellular_cluster_clones_states_summary_perc.pdf", dpi=320, width=12, height=7)

ut.cellular.states.ss <- tmp2 %>%
  filter(count == 1 & Condition == "Untreated") %>% ungroup() %>%
  select(c(real_bc44, Sample))

ut.cellular.states.ss.df <- merge(ut.cellular.states.ss, tmp, by=c("real_bc44", "Sample"))
length(unique(ut.cellular.states.ss.df$real_bc44)) #151 clones

ut.cellular.states.ss.T2 <- ut.cellular.states.ss.df %>%
  group_by(real_bc44, Sample, timepoint) %>%
  summarise(barcodes = n_distinct(real_bc44))
ut.cellular.states.ss.T2$clone_state <- NA

dtp.ut.meta.clones <- ut.cellular.states.ss.T2 %>%
  mutate(clone_state = ifelse(barcodes == "1" & timepoint == "T2", "DTP", clone_state))

dtp.ut.meta.clones.filtered <- dtp.ut.meta.clones %>%
  filter(clone_state == "DTP") %>%
  select(-c(timepoint, barcodes))

ut.cellular.dtp.df <- merge(dtp.ut.meta.clones.filtered, tmp, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp2, ut.cellular.states, ut.cellular.states.perc, ut.cellular.states.ss,
   ut.cellular.states.ss.df, ut.cellular.states.ss.T2, ut.cellular.states.summary,
   dtp.ut.meta.clones, dtp.ut.meta.clones.filtered, dim1.gg, dim2.gg, dim3.gg)
gc()

#Visualise for ggalluvial
#Remove datapoints with no information at T2
ut.dtp.A.df <- ut.cellular.dtp.df %>% filter(Sample == "Neuro_A") %>% select(c(real_bc44, timepoint, Count, seurat_clusters.0.2, Sample))

tmp <- ut.dtp.A.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = seurat_clusters.0.2) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

ut.dtp.cellular.A.df <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))

#Visualise for ggalluvial
ut.dtp.B.df <- ut.cellular.dtp.df %>% filter(Sample == "Neuro_B") %>% select(c(real_bc44, timepoint, Count, seurat_clusters.0.2, Sample))

tmp <- ut.dtp.B.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = seurat_clusters.0.2) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

ut.dtp.cellular.B.df <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))

#Visualise for ggalluvial
ut.dtp.C.df <- ut.cellular.dtp.df %>% filter(Sample == "Neuro_C") %>% select(c(real_bc44, timepoint, Count, seurat_clusters.0.2, Sample))

tmp <- ut.dtp.C.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = seurat_clusters.0.2) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

ut.dtp.cellular.C.df <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp.T1, tmp.T2, ut.dtp.A, ut.dtp.B, ut.dtp.C,
   distinct.ut.A.barcodes, distinct.ut.B.barcodes, distinct.ut.C.barcodes, ut.all.barcodes,
   ut.T2.barcodes, ut.dtp.A.df, ut.dtp.B.df, ut.dtp.C.df)
gc()

#Plot data using alluvial
#Generate UMAP including all info in wave plots from DTPs
ut.dtp.cellular.A.df <- ut.dtp.cellular.A.df %>% filter(Sample == "Neuro_A")
ut.dtp.cellular.B.df <- ut.dtp.cellular.B.df %>% filter(Sample == "Neuro_B")
ut.dtp.cellular.C.df <- ut.dtp.cellular.C.df %>% filter(Sample == "Neuro_C")

#Find overlapping clones in each replicate
A <- ut.dtp.cellular.A.df$real_bc44
B <- ut.dtp.cellular.B.df$real_bc44
C <- ut.dtp.cellular.C.df$real_bc44
clone.list <- Reduce(intersect, list(A, B, C))

test <- rbind(ut.dtp.cellular.A.df, ut.dtp.cellular.B.df, ut.dtp.cellular.C.df)
all.ut.cellular.dtp.df <- left_join(test, full.meta.df, by=c("real_bc44", "Sample"))

#Clean up files
rm(ut.dtp.cellular.A.df, ut.dtp.cellular.B.df, ut.dtp.cellular.C.df)
gc()

#Generate colour palette for seurat_clusters_0.2
P9 <- c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA")
P9 <- as.vector(t(matrix(P9)))
names(P9) = c("0", "1", "2", "3", "4", "5", "6", "7", "8")


####
#UT AMT states
ut.amt.dtp.df <- full.meta.df %>%
  filter(Condition == "Untreated" | Condition == "Untreated_rec") %>%
  select(c(AMT.state, Sample, real_bc44, Condition, Count))

ut.amt.dtp.df <- ut.amt.dtp.df[, c(3, 4, 1, 2, 5)]

tmp <- ut.amt.dtp.df %>%
  group_by(Sample, Condition, real_bc44) %>%
  mutate(total_clones=sum(Count)) %>%
  ungroup() %>%
  group_by(real_bc44, Condition, Sample) %>%
  mutate(percentage = Count/total_clones * 100) %>%
  ungroup() 

tmp$timepoint[tmp$Condition == "Untreated"] <- "T1"
tmp$timepoint[tmp$Condition == "Untreated_rec"] <- "T2"

tmp2 <- tmp %>%
  group_by(real_bc44, Condition, Sample) %>%
  summarise(count = n_distinct(AMT.state))

ut.amt.states <- tmp2 %>%
  group_by(count, Condition, Sample) %>%
  summarise(clones = length(unique(real_bc44)))

ut.amt.states.summary <- ut.amt.states %>%
  filter(Condition == "Untreated_rec")

#Calculate percentage of clones with each representation
ut.amt.states.perc <- ut.amt.states.summary %>%
  group_by(Sample, Condition) %>%
  mutate(total_freq = sum(clones)) %>%
  ungroup() %>%
  mutate(perc = clones/total_freq * 100) %>%
  ungroup()

set.seed(12)
P3 = createPalette(3, c("#ff0000", "#00ff00", "#0000ff"), M=10)
P3 <- as.vector(t(matrix(P3)))
names(P3) = unique(as.character(ut.amt.states.perc$count))

ut.amt.states.perc$count <- as.character(ut.amt.states.perc$count)

#Plot manuscript exteneded figure 7c - UT
ut.amt.states.perc %>%
  ggplot( ,aes(x=Sample, y=perc))+
  geom_bar(aes(x=Sample, y=perc, fill=count), colour="black", stat = "identity", position = "stack") +
  xlab("States Occupied")+
  ylab("Number of Clones")+
  ggtitle("UT Recovery AMT states Clone Summary")+
  scale_fill_manual(values = P3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))
ggsave("manuscript_plots/ut_amt_states_clones_states_summary_perc.pdf", dpi=320, width=12, height=7)

#Define clones which are only observed in one cellular state within untreated
ut.amt.states.ss <- tmp2 %>%
  filter(count == 1 & Condition == "Untreated") %>% ungroup() %>%
  select(c(real_bc44, Sample))

ut.amt.states.ss.df <- merge(ut.amt.states.ss, tmp, by=c("real_bc44", "Sample"))
length(unique(ut.amt.states.ss.df$real_bc44)) #315 clones

ut.amt.states.ss.T2 <- ut.amt.states.ss.df %>%
  group_by(real_bc44, Sample, timepoint) %>%
  summarise(barcodes = n_distinct(real_bc44))
ut.amt.states.ss.T2$clone_state <- NA

dtp.ut.amt.clones <- ut.amt.states.ss.T2 %>%
  mutate(clone_state = ifelse(barcodes == "1" & timepoint == "T2", "DTP", clone_state))

dtp.ut.amt.clones.filtered <- dtp.ut.amt.clones %>%
  filter(clone_state == "DTP") %>%
  select(-c(timepoint, barcodes))

ut.amt.dtp.df <- merge(dtp.ut.amt.clones.filtered, tmp, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp2, ut.amt.states, ut.amt.states.perc, ut.amt.states.ss,
   ut.amt.states.ss.df, ut.amt.states.ss.T2, ut.amt.states.summary,
   dtp.ut.amt.clones, dtp.ut.amt.clones.filtered)
gc()

#Visualise for ggalluvial
ut.dtp.A.df <- ut.amt.dtp.df %>% filter(Sample == "Neuro_A") %>% select(c(real_bc44, timepoint, Count, AMT.state, Sample))

tmp <- ut.dtp.A.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = AMT.state) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

ut.dtp.amt.A.df <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))

#Visualise for ggalluvial
ut.dtp.B.df <- ut.amt.dtp.df %>% filter(Sample == "Neuro_B") %>% select(c(real_bc44, timepoint, Count, AMT.state, Sample))

tmp <- ut.dtp.B.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = AMT.state) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

ut.dtp.amt.B.df <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))

#Visualise for ggalluvial
ut.dtp.C.df <- ut.amt.dtp.df %>% filter(Sample == "Neuro_C") %>% select(c(real_bc44, timepoint, Count, AMT.state, Sample))

tmp <- ut.dtp.C.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = AMT.state) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

ut.dtp.amt.C.df <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp.T1, tmp.T2, ut.dtp.A, ut.dtp.B, ut.dtp.C,
   distinct.ut.A.barcodes, distinct.ut.B.barcodes, distinct.ut.C.barcodes, ut.all.barcodes,
   ut.T2.barcodes, ut.dtp.A.df, ut.dtp.B.df, ut.dtp.C.df)
gc()

#Plot data using alluvial
#Generate UMAP including all info in wave plots from DTPs
ut.dtp.amt.A.df <- ut.dtp.amt.A.df %>% filter(Sample == "Neuro_A")
ut.dtp.amt.B.df <- ut.dtp.amt.B.df %>% filter(Sample == "Neuro_B")
ut.dtp.amt.C.df <- ut.dtp.amt.C.df %>% filter(Sample == "Neuro_C")

#Find overlapping clones in each replicate
A <- ut.dtp.amt.A.df$real_bc44
B <- ut.dtp.amt.B.df$real_bc44
C <- ut.dtp.amt.C.df$real_bc44
clone.list <- Reduce(intersect, list(A, B, C))

test <- rbind(ut.dtp.amt.A.df, ut.dtp.amt.B.df, ut.dtp.amt.C.df)
all.ut.amt.dtp.df <- left_join(test, full.meta.df, by=c("real_bc44", "Sample"))

#Clean up files
rm(ut.dtp.amt.A.df, ut.dtp.amt.B.df, ut.dtp.amt.C.df)
gc()

#To counteract the overlap of cloens between replicates
#Make new variable which is barcode + replicate
all.ut.amt.dtp.df <- all.ut.amt.dtp.df %>% mutate(barcode_rep = paste0(real_bc44, "_", Sample))
all.ut.cellular.dtp.df <- all.ut.cellular.dtp.df %>%mutate(barcode_rep = paste0(real_bc44, "_", Sample))

# scRNA-seq barcodes - Alluvial Untreated Plots 
#Visulaise all informaito for cellular clusters and AMT states
all.ut.amt.dtp.df %>%
  filter(Condition == "Untreated") %>%
  mutate(Condition = factor(Condition, levels=c("Untreated"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "grey", pointsize = 2)+
  geom_scattermore(aes(colour = AMT.state), pointsize = 3)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(breaks = c("MES", "intermediate", "ADRN"),
                      values=c("#F37735", "lightgrey", "#990099")) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))
ggsave("manuscript_plots/ut_dtps_amt_umap.pdf", dpi = 700, height=7, width =12)

all.ut.cellular.dtp.df %>%
  filter(Condition == "Untreated") %>%
  mutate(Condition = factor(Condition, levels=c("Untreated"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "grey", pointsize = 2)+
  geom_scattermore(aes(colour = seurat_clusters.0.2), pointsize = 3)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(values=P9) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))
ggsave("manuscript_plots/ut_dtps_cellular_umap.pdf", dpi = 700, height=7, width =12)

all.ut.amt.dtp.df %>%
  filter(Condition == "Untreated_rec") %>%
  mutate(Condition = factor(Condition, levels=c("Untreated_rec"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "grey", pointsize = 2)+
  geom_scattermore(aes(colour = AMT.state), pointsize = 3)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(breaks = c("MES", "intermediate", "ADRN"),
                      values=c("#F37735", "lightgrey", "#990099")) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))
ggsave("manuscript_plots/ut_rec_dtps_amt_umap.pdf", dpi = 700, height=7, width =12)

all.ut.cellular.dtp.df %>%
  filter(Condition == "Untreated_rec") %>%
  mutate(Condition = factor(Condition, levels=c("Untreated_rec"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "grey", pointsize = 2)+
  geom_scattermore(aes(colour = seurat_clusters.0.2), pointsize = 3)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(values=P9) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))
ggsave("manuscript_plots/ut_rec_dtps_cellular_umap.pdf", dpi = 700, height=7, width =12)


#Make compiled list of clones to generate JQ1 colour palette for barcodes
amt.list <- unique(all.ut.amt.dtp.df$barcode_rep)
cellular.list <- unique(all.ut.cellular.dtp.df$barcode_rep)
ut.dtp.barcodes <- c(amt.list, cellular.list)
unique.ut.dtp.barcodes <- unique(ut.dtp.barcodes)

#Generate barcode palette
set.seed(12)
ut_palette = createPalette(234, c("#ff0000", "#00ff00", "#0000ff"))
ut_palette <- as.vector(t(matrix(ut_palette)))
names(ut_palette) = unique(as.character(ut.dtp.barcodes))

#Plot manuscript extended figure 7e
all.ut.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Untreated"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill=T1)) +
  geom_stratum( width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())
ggsave("manuscript_plots/ut_alluvial_cellular_barcode+rep_1.pdf", dpi=700, width=12, height=7)

#Plot manuscript extended figure 7e
all.ut.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Untreated_rec"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T1), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())
ggsave("manuscript_plots/ut_alluvial_cellular_barcode+rep_2.pdf", dpi=320, width=12, height=7)

#Plot manuscript extended figure 7e
all.ut.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Untreated_rec"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T2), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())
ggsave("manuscript_plots/ut_alluvial_cellular_barcode+rep_3.pdf", dpi=320, width=12, height=7)

#Plot manuscript extended figure 7e
all.ut.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Untreated_rec"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(width = 1/8, colour="black", aes(fill=barcode_rep), linewidth=0) +
  scale_fill_manual(values = ut_palette) +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank(),
        legend.position = "none")
ggsave("manuscript_plots/ut_alluvial_cellular_barcode+rep_4.pdf", dpi=320, width=12, height=7)

#AMT
#Plot manuscript extended figure 7f
all.ut.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Untreated_rec"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill=T1)) +
  geom_stratum( width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())
ggsave("manuscript_plots/ut_alluvial_amt_barcode+rep_1.pdf", dpi=700, width=12, height=7)

#Plot manuscript extended figure 7f
all.ut.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Untreated_rec"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T1), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())
ggsave("manuscript_plots/ut_alluvial_amt_barcode+rep_2.pdf", dpi=320, width=12, height=7)

#Plot manuscript extended figure 7f
all.ut.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Untreated_rec"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T2), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())
ggsave("manuscript_plots/ut_alluvial_amt_barcode+rep_3.pdf", dpi=320, width=12, height=7)

#Plot manuscript extended figure 7f
all.ut.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Untreated_rec"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(width = 1/8, colour="black", aes(fill=barcode_rep), linewidth=0) +
  scale_fill_manual(values = ut_palette) +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank(),
        legend.position = "none")
ggsave("manuscript_plots/ut_alluvial_amt_barcode+rep_5.pdf", dpi=320, width=12, height=7)

# scRNA-seq barcodes - Untreated Naive Dynamics
ut.dtp.df <- full.meta.df %>%
  filter(Condition == "Untreated") %>%
  select(c(AMT.state, Sample, real_bc44, Condition, Count))

ut.dtp.df <- ut.dtp.df[, c(3, 4, 1, 2, 5)]

tmp <- ut.dtp.df %>%
  group_by(Sample, Condition, real_bc44) %>%
  mutate(total_clones=sum(Count)) %>%
  ungroup() %>%
  group_by(real_bc44, Condition, Sample) %>%
  mutate(percentage = Count/total_clones * 100) %>%
  ungroup() 

tmp$timepoint[tmp$Condition == "Untreated"] <- "T1"

tmp2 <- tmp %>%
  group_by(real_bc44, Condition, Sample) %>%
  summarise(count = n_distinct(AMT.state))

ut.states <- tmp2 %>%
  group_by(count, Condition, Sample) %>%
  summarise(clones = length(unique(real_bc44)))

ut.states.perc <- ut.states %>%
  group_by(Sample, Condition) %>%
  mutate(total_freq = sum(clones)) %>%
  ungroup() %>%
  mutate(perc = clones/total_freq * 100) %>%
  ungroup()

ut.states.perc$count <- as.character(ut.states.perc$count)

#Plot manuscript extended figure 7a
ut.states.perc %>%
  ggplot( ,aes(x=Sample, y=perc))+
  geom_bar(aes(x=Sample, y=perc, fill=count), colour="black", stat = "identity", position = "stack") +
  xlab("States Occupied")+
  ylab("Number of Clones")+
  ggtitle("UT AMT states Clone Summary")+
  scale_fill_manual(values = P3)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))

ut.states.perc <- ut.states.perc %>%
  group_by(count) %>%
  mutate(avg_perc = mean(perc))

library(ggforce)
#Plot manuscript extended figure 7b
ut.states.perc %>%
  ggplot( ,aes(x=count, y=perc))+
  geom_boxplot(aes(x=count, y=perc, colour=count)) +
  xlab("States Occupied")+
  ylab("Number of Clones")+
  ggtitle("UT AMT states Clone Summary")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        aspect.ratio = 2)

rm(ut.cellular.dtp.df, ut.dtp.df, ut.states, ut.states.perc,
   tmp, tmp2, test, unique.ut.dtp.barcodes, ut.dtp.barcodes)
gc()


# scRNA-seq barcodes - Alluvial Cisplatin Initial
cisplatin.cellular.dtp.df <- full.meta.df %>%
  filter(Condition == "Untreated" | Condition == "Cisplatin" | Condition == "Cisplatin_rec") %>%
  select(c(seurat_clusters.0.2, Sample, real_bc44, Condition, Count))

cisplatin.cellular.dtp.df <- cisplatin.cellular.dtp.df[, c(3, 4, 1, 2, 5)]

tmp <- cisplatin.cellular.dtp.df %>%
  group_by(Sample, Condition, real_bc44) %>%
  mutate(total_clones=sum(Count)) %>%
  ungroup() %>%
  group_by(real_bc44, Condition, Sample) %>%
  mutate(percentage = Count/total_clones * 100) %>%
  ungroup() 

tmp$timepoint[tmp$Condition == "Untreated"] <- "T1"
tmp$timepoint[tmp$Condition == "Cisplatin"] <- "T2"
tmp$timepoint[tmp$Condition == "Cisplatin_rec"] <- "T3"

tmp2 <- tmp %>%
  group_by(real_bc44, Condition, Sample) %>%
  summarise(count = n_distinct(seurat_clusters.0.2))

cis.cellular.states <- tmp2 %>%
  group_by(count, Condition, Sample) %>%
  summarise(clones = length(unique(real_bc44)))

cis.cellular.states.summary <- cis.cellular.states %>%
  filter(Condition == "Cisplatin_rec")

#Calculate percentage of clones with each representation
cis.cellular.states.perc <- cis.cellular.states.summary %>%
  group_by(Sample, Condition) %>%
  mutate(total_freq = sum(clones)) %>%
  ungroup() %>%
  mutate(perc = clones/total_freq * 100) %>%
  ungroup()

#Define colour palette for states
set.seed(12)
P8 = createPalette(8, c("#ff0000", "#00ff00", "#0000ff"), M=10)
P8 <- as.vector(t(matrix(P8)))
names(P8) = unique(as.character(cis.cellular.states.perc$count))

cis.cellular.states.perc$count <- as.character(cis.cellular.states.perc$count)

#Plot manuscript extended figure 7d - Cisplatin
cis.cellular.states.perc %>%
  ggplot( ,aes(x=Sample, y=perc))+
  geom_bar(aes(x=Sample, y=perc, fill=count), colour="black", stat = "identity", position = "stack") +
  xlab("States Occupied")+
  ylab("Number of Clones")+
  ggtitle("Cisplatin Recovery Cellular Cluster Clone Summary")+
  scale_fill_manual(values = P8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))

#Define clones which are only observed in one cellular state within untreated
cis.cellular.states.ss <- tmp2 %>%
  filter(count == 1 & Condition == "Untreated") %>% ungroup() %>%
  select(c(real_bc44, Sample))

cis.cellular.states.ss.df <- merge(cis.cellular.states.ss, tmp, by=c("real_bc44", "Sample"))
length(unique(cis.cellular.states.ss.df$real_bc44)) #316 clones

cis.cellular.states.ss.T2 <- cis.cellular.states.ss.df %>%
  group_by(real_bc44, Sample, timepoint) %>%
  summarise(barcodes = n_distinct(real_bc44))
cis.cellular.states.ss.T2$clone_state <- NA

dtp.cis.meta.clones <- cis.cellular.states.ss.T2 %>%
  mutate(clone_state = ifelse(barcodes == "1" & timepoint == "T3", "DTP", clone_state))

dtp.cis.meta.clones.filtered <- dtp.cis.meta.clones %>%
  filter(clone_state == "DTP") %>%
  select(-c(timepoint, barcodes))

cisplatin.cellular.dtp.df <- merge(dtp.cis.meta.clones.filtered, tmp, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp2, cis.cellular.states, cis.cellular.states.perc, cis.cellular.states.ss,
   cis.cellular.states.ss.df, cis.cellular.states.ss.T2, cis.cellular.states.summary,
   dtp.cis.meta.clones, dtp.cis.meta.clones.filtered)
gc()

#Remove datapoints with no information at T2
cis.T2.barcodes <- cisplatin.cellular.dtp.df$real_bc44[cisplatin.cellular.dtp.df$Sample == "Neuro_A" & cisplatin.cellular.dtp.df$timepoint == "T2"]
cis.all.barcodes <- cisplatin.cellular.dtp.df$real_bc44[cisplatin.cellular.dtp.df$Sample == "Neuro_A"]
distinct.cis.A.barcodes <- c(setdiff(cis.all.barcodes, cis.T2.barcodes), setdiff(cis.T2.barcodes, cis.all.barcodes))
cisplatin.dtp.A <- cisplatin.cellular.dtp.df %>% filter(!(real_bc44 %in% distinct.cis.A.barcodes))

#Visualise for ggalluvial
cisplatin.dtp.A.df <- cisplatin.dtp.A %>% select(c(real_bc44, timepoint, Count, seurat_clusters.0.2, Sample))

tmp <- cisplatin.dtp.A.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = seurat_clusters.0.2) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
cisplatin.dtp.cellular.A.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Remove datapoints with no information at T2
cis.T2.barcodes <- cisplatin.cellular.dtp.df$real_bc44[cisplatin.cellular.dtp.df$Sample == "Neuro_B" & cisplatin.cellular.dtp.df$timepoint == "T2"]
cis.all.barcodes <- cisplatin.cellular.dtp.df$real_bc44[cisplatin.cellular.dtp.df$Sample == "Neuro_B"]
distinct.cis.B.barcodes <- c(setdiff(cis.all.barcodes, cis.T2.barcodes), setdiff(cis.T2.barcodes, cis.all.barcodes))
cisplatin.dtp.B <- cisplatin.cellular.dtp.df %>% filter(!(real_bc44 %in% distinct.cis.B.barcodes))

#Visualise for ggalluvial
cisplatin.dtp.B.df <- cisplatin.dtp.B %>% select(c(real_bc44, timepoint, Count, seurat_clusters.0.2, Sample))

tmp <- cisplatin.dtp.B.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = seurat_clusters.0.2) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
cisplatin.dtp.cellular.B.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Remove datapoints with no information at T2
cis.T2.barcodes <- cisplatin.cellular.dtp.df$real_bc44[cisplatin.cellular.dtp.df$Sample == "Neuro_C" & cisplatin.cellular.dtp.df$timepoint == "T2"]
cis.all.barcodes <- cisplatin.cellular.dtp.df$real_bc44[cisplatin.cellular.dtp.df$Sample == "Neuro_C"]
distinct.cis.C.barcodes <- c(setdiff(cis.all.barcodes, cis.T2.barcodes), setdiff(cis.T2.barcodes, cis.all.barcodes))
cisplatin.dtp.C <- cisplatin.cellular.dtp.df %>% filter(!(real_bc44 %in% distinct.cis.C.barcodes))

#Visualise for ggalluvial
cisplatin.dtp.C.df <- cisplatin.dtp.C %>% select(c(real_bc44, timepoint, Count, seurat_clusters.0.2, Sample))

tmp <- cisplatin.dtp.C.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = seurat_clusters.0.2) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
cisplatin.dtp.cellular.C.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp.T1, tmp.T2, tmp.T3, tmp3, test.merge, cisplatin.dtp.A, cisplatin.dtp.B, cisplatin.dtp.C,
   distinct.cis.A.barcodes, distinct.cis.B.barcodes, distinct.cis.C.barcodes, cis.all.barcodes,
   cis.T2.barcodes, cisplatin.dtp.A.df, cisplatin.dtp.B.df, cisplatin.dtp.C.df)
gc()

#Plot data using alluvial
#Generate UMAP including all info in wave plots from DTPs
cisplatin.dtp.cellular.A.df <- cisplatin.dtp.cellular.A.df %>% filter(Sample == "Neuro_A")
cisplatin.dtp.cellular.B.df <- cisplatin.dtp.cellular.B.df %>% filter(Sample == "Neuro_B")
cisplatin.dtp.cellular.C.df <- cisplatin.dtp.cellular.C.df %>% filter(Sample == "Neuro_C")

#Find overlapping clones in each replicate
A <- cisplatin.dtp.cellular.A.df$real_bc44
B <- cisplatin.dtp.cellular.B.df$real_bc44
C <- cisplatin.dtp.cellular.C.df$real_bc44
clone.list <- Reduce(intersect, list(A, B, C))

test <- rbind(cisplatin.dtp.cellular.A.df, cisplatin.dtp.cellular.B.df, cisplatin.dtp.cellular.C.df)
all.cisplatin.cellular.dtp.df <- left_join(test, full.meta.df, by=c("real_bc44", "Sample"))

#Clean up files
rm(cisplatin.dtp.cellular.A.df, cisplatin.dtp.cellular.B.df, cisplatin.dtp.cellular.C.df, cisplatin.dtp.cellular.df,
   cisplatin.dtp.cellular, cisplatin.dtp.filtered, dtp.cis.meta.clones.filtered)
gc()

#Generate colour palette for seurat_clusters_0.2
P9 <- c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA")
P9 <- as.vector(t(matrix(P9)))
names(P9) = c("0", "1", "2", "3", "4", "5", "6", "7", "8")

####
#Cisplatin AMT states
cisplatin.amt.dtp.df <- full.meta.df %>%
  filter(Condition == "Untreated" | Condition == "Cisplatin" | Condition == "Cisplatin_rec") %>%
  select(c(AMT.state, Sample, real_bc44, Condition, Count))

cisplatin.amt.dtp.df <- cisplatin.amt.dtp.df[, c(3, 4, 1, 2, 5)]

tmp <- cisplatin.amt.dtp.df %>%
  group_by(Sample, Condition, real_bc44) %>%
  mutate(total_clones=sum(Count)) %>%
  ungroup() %>%
  group_by(real_bc44, Condition, Sample) %>%
  mutate(percentage = Count/total_clones * 100) %>%
  ungroup() 

tmp$timepoint[tmp$Condition == "Untreated"] <- "T1"
tmp$timepoint[tmp$Condition == "Cisplatin"] <- "T2"
tmp$timepoint[tmp$Condition == "Cisplatin_rec"] <- "T3"

tmp2 <- tmp %>%
  group_by(real_bc44, Condition, Sample) %>%
  summarise(count = n_distinct(AMT.state))

cis.amt.states <- tmp2 %>%
  group_by(count, Condition, Sample) %>%
  summarise(clones = length(unique(real_bc44)))

cis.amt.states.summary <- cis.amt.states %>%
  filter(Condition == "Cisplatin_rec")

#Calculate percentage of clones with each representation
cis.amt.states.perc <- cis.amt.states.summary %>%
  group_by(Sample, Condition) %>%
  mutate(total_freq = sum(clones)) %>%
  ungroup() %>%
  mutate(perc = clones/total_freq * 100) %>%
  ungroup()

#Define colour palette for states
set.seed(12)
P3 = createPalette(3, c("#ff0000", "#00ff00", "#0000ff"), M=10)
P3 <- as.vector(t(matrix(P3)))
names(P3) = unique(as.character(cis.amt.states.perc$count))

cis.amt.states.perc$count <- as.character(cis.amt.states.perc$count)

#Plot manuscript extanded figure 7c - Cisplatin
cis.amt.states.perc %>%
  ggplot( ,aes(x=Sample, y=perc))+
  geom_bar(aes(x=Sample, y=perc, fill=count), colour="black", stat = "identity", position = "stack") +
  scale_fill_manual(values = P3)+
  xlab("States Occupied")+
  ylab("Number of Clones")+
  ggtitle("Cisplatin Recovery AMT states Clone Summary")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))

#Define clones which are only observed in one cellular state within untreated
cis.amt.states.ss <- tmp2 %>%
  filter(count == 1 & Condition == "Untreated") %>% ungroup() %>%
  select(c(real_bc44, Sample))

cis.amt.states.ss.df <- merge(cis.amt.states.ss, tmp, by=c("real_bc44", "Sample"))
length(unique(cis.amt.states.ss.df$real_bc44)) #154 clones

cis.amt.states.ss.T2 <- cis.amt.states.ss.df %>%
  group_by(real_bc44, Sample, timepoint) %>%
  summarise(barcodes = n_distinct(real_bc44))
cis.amt.states.ss.T2$clone_state <- NA

dtp.cis.amt.clones <- cis.amt.states.ss.T2 %>%
  mutate(clone_state = ifelse(barcodes == "1" & timepoint == "T3", "DTP", clone_state))

dtp.cis.amt.clones.filtered <- dtp.cis.amt.clones %>%
  filter(clone_state == "DTP") %>%
  select(-c(timepoint, barcodes))

cisplatin.amt.dtp.df <- merge(dtp.cis.amt.clones.filtered, tmp, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp2, cis.amt.states, cis.amt.states.perc, cis.amt.states.ss,
   cis.amt.states.ss.df, cis.amt.states.ss.T2, cis.amt.states.summary,
   dtp.cis.amt.clones, dtp.cis.amt.clones.filtered)
gc()

#Remove datapoints with no information at T2
cis.T2.barcodes <- cisplatin.amt.dtp.df$real_bc44[cisplatin.amt.dtp.df$Sample == "Neuro_A" & cisplatin.amt.dtp.df$timepoint == "T2"]
cis.all.barcodes <- cisplatin.amt.dtp.df$real_bc44[cisplatin.amt.dtp.df$Sample == "Neuro_A"]
distinct.cis.A.barcodes <- c(setdiff(cis.all.barcodes, cis.T2.barcodes), setdiff(cis.T2.barcodes, cis.all.barcodes))
cisplatin.dtp.A <- cisplatin.amt.dtp.df %>% filter(!(real_bc44 %in% distinct.cis.A.barcodes))

#Visualise for ggalluvial
cisplatin.dtp.A.df <- cisplatin.dtp.A %>% select(c(real_bc44, timepoint, Count, AMT.state, Sample))

tmp <- cisplatin.dtp.A.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = AMT.state) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
cisplatin.dtp.amt.A.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Remove datapoints with no information at T2
cis.T2.barcodes <- cisplatin.amt.dtp.df$real_bc44[cisplatin.amt.dtp.df$Sample == "Neuro_B" & cisplatin.amt.dtp.df$timepoint == "T2"]
cis.all.barcodes <- cisplatin.amt.dtp.df$real_bc44[cisplatin.amt.dtp.df$Sample == "Neuro_B"]
distinct.cis.B.barcodes <- c(setdiff(cis.all.barcodes, cis.T2.barcodes), setdiff(cis.T2.barcodes, cis.all.barcodes))
cisplatin.dtp.B <- cisplatin.amt.dtp.df %>% filter(!(real_bc44 %in% distinct.cis.B.barcodes))

#Visualise for ggalluvial
cisplatin.dtp.B.df <- cisplatin.dtp.B %>% select(c(real_bc44, timepoint, Count, AMT.state, Sample))

tmp <- cisplatin.dtp.B.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = AMT.state) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
cisplatin.dtp.amt.B.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Remove datapoints with no information at T2
cis.T2.barcodes <- cisplatin.amt.dtp.df$real_bc44[cisplatin.amt.dtp.df$Sample == "Neuro_C" & cisplatin.amt.dtp.df$timepoint == "T2"]
cis.all.barcodes <- cisplatin.amt.dtp.df$real_bc44[cisplatin.amt.dtp.df$Sample == "Neuro_C"]
distinct.cis.C.barcodes <- c(setdiff(cis.all.barcodes, cis.T2.barcodes), setdiff(cis.T2.barcodes, cis.all.barcodes))
cisplatin.dtp.C <- cisplatin.amt.dtp.df %>% filter(!(real_bc44 %in% distinct.cis.C.barcodes))

#Visualise for ggalluvial
cisplatin.dtp.C.df <- cisplatin.dtp.C %>% select(c(real_bc44, timepoint, Count, AMT.state, Sample))

tmp <- cisplatin.dtp.C.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = AMT.state) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
cisplatin.dtp.amt.C.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp.T1, tmp.T2, tmp.T3, test.merge, cisplatin.dtp.A, cisplatin.dtp.B, cisplatin.dtp.C,
   distinct.cis.A.barcodes, distinct.cis.B.barcodes, distinct.cis.C.barcodes, cis.all.barcodes,
   cis.T2.barcodes, cisplatin.dtp.A.df, cisplatin.dtp.B.df, cisplatin.dtp.C.df)
gc()

#Plot data using alluvial
#Generate UMAP including all info in wave plots from DTPs
cisplatin.dtp.amt.A.df <- cisplatin.dtp.amt.A.df %>% filter(Sample == "Neuro_A")
cisplatin.dtp.amt.B.df <- cisplatin.dtp.amt.B.df %>% filter(Sample == "Neuro_B")
cisplatin.dtp.amt.C.df <- cisplatin.dtp.amt.C.df %>% filter(Sample == "Neuro_C")

#Find overlapping clones in each replicate
A <- cisplatin.dtp.amt.A.df$real_bc44
B <- cisplatin.dtp.amt.B.df$real_bc44
C <- cisplatin.dtp.amt.C.df$real_bc44
clone.list <- Reduce(intersect, list(A, B, C))

test <- rbind(cisplatin.dtp.amt.A.df, cisplatin.dtp.amt.B.df, cisplatin.dtp.amt.C.df)
all.cisplatin.amt.dtp.df <- left_join(test, full.meta.df, by=c("real_bc44", "Sample"))

#Clean up files
rm(cisplatin.dtp.amt.A.df, cisplatin.dtp.amt.B.df, cisplatin.dtp.amt.C.df, cisplatin.dtp.amt.df,
   cisplatin.dtp.amt, cisplatin.dtp.filtered, dtp.cis.meta.clones.filtered)
gc()

#To counteract the overlap of cloens between replicates
#Make new variable which is barcode + replicate
all.cisplatin.amt.dtp.df <- all.cisplatin.amt.dtp.df %>%
  mutate(barcode_rep = paste0(real_bc44, "_", Sample))

all.cisplatin.cellular.dtp.df <- all.cisplatin.cellular.dtp.df %>%
  mutate(barcode_rep = paste0(real_bc44, "_", Sample))

#Visulaise all information for cellular clusters and AMT states
#Plot manuscript Figure 4c - left hand panel
all.cisplatin.amt.dtp.df %>%
  filter(Condition == "Untreated") %>%
  mutate(Condition = factor(Condition, levels=c("Untreated"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "white", pointsize = 2)+
  geom_scattermore(aes(colour = AMT.state), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(breaks = c("MES", "intermediate", "ADRN"),
                      values=c("#F37735", "lightgrey", "#990099")) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Plot manuscript Figure 4b - left hand panel
all.cisplatin.cellular.dtp.df %>%
  filter(Condition == "Untreated") %>%
  mutate(Condition = factor(Condition, levels=c("Untreated"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "white", pointsize = 2)+
  geom_scattermore(aes(colour = seurat_clusters.0.2), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(values= P9) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Plot manuscript Figure 4c - right hand panel
all.cisplatin.amt.dtp.df %>%
  filter(Condition == "Cisplatin_rec") %>%
  mutate(Condition = factor(Condition, levels=c("Cisplatin_rec"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "white", pointsize = 2)+
  geom_scattermore(aes(colour = AMT.state), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(breaks = c("MES", "intermediate", "ADRN"),
                      values=c("#F37735", "lightgrey", "#990099")) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Plot manuscript Figure 4b - right hand panel
all.cisplatin.cellular.dtp.df %>%
  filter(Condition == "Cisplatin_rec") %>%
  mutate(Condition = factor(Condition, levels=c("Cisplatin_rec"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "white", pointsize = 2)+
  geom_scattermore(aes(colour = seurat_clusters.0.2), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(values=P9) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Make compiled list of clones to generate Cisplatin colour palette for barcodes
amt.list <- unique(all.cisplatin.amt.dtp.df$barcode_rep)
cellular.list <- unique(all.cisplatin.cellular.dtp.df$barcode_rep)
cisplatin.dtp.barcodes <- c(amt.list, cellular.list)
unique.cisplatin.dtp.barcodes <- unique(cisplatin.dtp.barcodes)

#Generate barcode palette
set.seed(12)
cisplatin_palette = createPalette(600, c("#ff0000", "#00ff00", "#0000ff"))
cisplatin_palette <- cisplatin_palette[!(cisplatin_palette %in% ut_palette)]
cisplatin_palette <- sample(cisplatin_palette, 79, replace=FALSE)
cisplatin_palette <- as.vector(t(matrix(cisplatin_palette)))
names(cisplatin_palette) = unique(as.character(cisplatin.dtp.barcodes))

#Plot manusctip extended figure 7e - Cisplatin
all.cisplatin.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill=T1)) +
  geom_stratum( width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manusctip extended figure 7e - Cisplatin
all.cisplatin.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T1), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manusctip extended figure 7e - Cisplatin
all.cisplatin.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T2), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manusctip extended figure 7e - Cisplatin
all.cisplatin.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=real_bc44, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T3), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manusctip extended figure 7e - Cisplatin
all.cisplatin.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(width = 1/8, colour="black", aes(fill=barcode_rep), linewidth=0) +
  scale_fill_manual(values = cisplatin_palette) +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank(),
        legend.position = "none")


#AMT
#Plot manusctip extended figure 7f - Cisplatin
all.cisplatin.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill=T1)) +
  geom_stratum( width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manusctip extended figure 7e - Cisplatin
all.cisplatin.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T1), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manusctip extended figure 7e - Cisplatin
all.cisplatin.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T2), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manusctip extended figure 7e - Cisplatin
all.cisplatin.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=real_bc44, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T3), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manusctip extended figure 7e - Cisplatin
all.cisplatin.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "Cisplatin", "Cisplatin Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(width = 1/8, colour="black", aes(fill=barcode_rep), linewidth=0) +
  scale_fill_manual(values = cisplatin_palette) +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank(),
        legend.position = "none")

rm(cisplatin.amt.dtp.df, cisplatin.cellular.dtp.df,
   cisplatin.dtp.barcodes, unique.cisplatin.dtp.barcodes)
gc()

# scRNA-seq barcodes - Alluvial JQ1 Plots 
jq1.cellular.dtp.df <- full.meta.df %>%
  filter(Condition == "Untreated" | Condition == "JQ1" | Condition == "JQ1_rec") %>%
  select(c(seurat_clusters.0.2, Sample, real_bc44, Condition, Count))

jq1.cellular.dtp.df <- jq1.cellular.dtp.df[, c(3, 4, 1, 2, 5)]

tmp <- jq1.cellular.dtp.df %>%
  group_by(Sample, Condition, real_bc44) %>%
  mutate(total_clones=sum(Count)) %>%
  ungroup() %>%
  group_by(real_bc44, Condition, Sample) %>%
  mutate(percentage = Count/total_clones * 100) %>%
  ungroup() 

tmp$timepoint[tmp$Condition == "Untreated"] <- "T1"
tmp$timepoint[tmp$Condition == "JQ1"] <- "T2"
tmp$timepoint[tmp$Condition == "JQ1_rec"] <- "T3"

tmp2 <- tmp %>%
  group_by(real_bc44, Condition, Sample) %>%
  summarise(count = n_distinct(seurat_clusters.0.2))

jq1.cellular.states <- tmp2 %>%
  group_by(count, Condition, Sample) %>%
  summarise(clones = length(unique(real_bc44)))

jq1.cellular.states.summary <- jq1.cellular.states %>%
  filter(Condition == "JQ1_rec")

#Calculate percentage of clones with each representation
jq1.cellular.states.perc <- jq1.cellular.states.summary %>%
  group_by(Sample, Condition) %>%
  mutate(total_freq = sum(clones)) %>%
  ungroup() %>%
  mutate(perc = clones/total_freq * 100) %>%
  ungroup()

jq1.cellular.states.perc$count <- as.character(jq1.cellular.states.perc$count)

#Plot manuscript extended figure 7d - JQ1
jq1.cellular.states.perc %>%
  ggplot( ,aes(x=Sample, y=perc))+
  geom_bar(aes(x=Sample, y=perc, fill=count), colour="black", stat = "identity", position = "stack") +
  xlab("States Occupied")+
  ylab("Number of Clones")+
  ggtitle("JQ1 Recovery Cellular Cluster Clone Summary")+
  scale_fill_manual(values = P8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))

#Define clones which are only observed in one cellular state within untreated
jq1.cellular.states.ss <- tmp2 %>%
  filter(count == 1 & Condition == "Untreated") %>% ungroup() %>%
  select(c(real_bc44, Sample))

jq1.cellular.states.ss.df <- merge(jq1.cellular.states.ss, tmp, by=c("real_bc44", "Sample"))
length(unique(jq1.cellular.states.ss.df$real_bc44)) #316 clones

jq1.cellular.states.ss.T2 <- jq1.cellular.states.ss.df %>%
  group_by(real_bc44, Sample, timepoint) %>%
  summarise(barcodes = n_distinct(real_bc44))
jq1.cellular.states.ss.T2$clone_state <- NA

dtp.jq1.meta.clones <- jq1.cellular.states.ss.T2 %>%
  mutate(clone_state = ifelse(barcodes == "1" & timepoint == "T3", "DTP", clone_state))

dtp.jq1.meta.clones.filtered <- dtp.jq1.meta.clones %>%
  filter(clone_state == "DTP") %>%
  select(-c(timepoint, barcodes))

jq1.cellular.dtp.df <- merge(dtp.jq1.meta.clones.filtered, tmp, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp2, jq1.cellular.states, jq1.cellular.states.perc, jq1.cellular.states.ss,
   jq1.cellular.states.ss.df, jq1.cellular.states.ss.T2, jq1.cellular.states.summary,
   dtp.jq1.meta.clones, dtp.jq1.meta.clones.filtered)
gc()

#Remove datapoints with no information at T2
jq1.T2.barcodes <- jq1.cellular.dtp.df$real_bc44[jq1.cellular.dtp.df$Sample == "Neuro_A" & jq1.cellular.dtp.df$timepoint == "T2"]
jq1.all.barcodes <- jq1.cellular.dtp.df$real_bc44[jq1.cellular.dtp.df$Sample == "Neuro_A"]
distinct.jq1.A.barcodes <- c(setdiff(jq1.all.barcodes, jq1.T2.barcodes), setdiff(jq1.T2.barcodes, jq1.all.barcodes))
jq1.dtp.A <- jq1.cellular.dtp.df %>% filter(!(real_bc44 %in% distinct.jq1.A.barcodes))

#Visualise for ggalluvial
jq1.dtp.A.df <- jq1.dtp.A %>% select(c(real_bc44, timepoint, Count, seurat_clusters.0.2, Sample))

tmp <- jq1.dtp.A.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = seurat_clusters.0.2) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
jq1.dtp.cellular.A.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Remove datapoints with no information at T2
jq1.T2.barcodes <- jq1.cellular.dtp.df$real_bc44[jq1.cellular.dtp.df$Sample == "Neuro_B" & jq1.cellular.dtp.df$timepoint == "T2"]
jq1.all.barcodes <- jq1.cellular.dtp.df$real_bc44[jq1.cellular.dtp.df$Sample == "Neuro_B"]
distinct.jq1.B.barcodes <- c(setdiff(jq1.all.barcodes, jq1.T2.barcodes), setdiff(jq1.T2.barcodes, jq1.all.barcodes))
jq1.dtp.B <- jq1.cellular.dtp.df %>% filter(!(real_bc44 %in% distinct.jq1.B.barcodes))

#Visualise for ggalluvial
jq1.dtp.B.df <- jq1.dtp.B %>% select(c(real_bc44, timepoint, Count, seurat_clusters.0.2, Sample))

tmp <- jq1.dtp.B.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = seurat_clusters.0.2) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
jq1.dtp.cellular.B.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Remove datapoints with no information at T2
jq1.T2.barcodes <- jq1.cellular.dtp.df$real_bc44[jq1.cellular.dtp.df$Sample == "Neuro_C" & jq1.cellular.dtp.df$timepoint == "T2"]
jq1.all.barcodes <- jq1.cellular.dtp.df$real_bc44[jq1.cellular.dtp.df$Sample == "Neuro_C"]
distinct.jq1.C.barcodes <- c(setdiff(jq1.all.barcodes, jq1.T2.barcodes), setdiff(jq1.T2.barcodes, jq1.all.barcodes))
jq1.dtp.C <- jq1.cellular.dtp.df %>% filter(!(real_bc44 %in% distinct.jq1.C.barcodes))

#Visualise for ggalluvial
jq1.dtp.C.df <- jq1.dtp.C %>% select(c(real_bc44, timepoint, Count, seurat_clusters.0.2, Sample))

tmp <- jq1.dtp.C.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = seurat_clusters.0.2) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
jq1.dtp.cellular.C.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp.T1, tmp.T2, tmp.T3, tmp3, test.merge, jq1.dtp.A, jq1.dtp.B, jq1.dtp.C,
   distinct.jq1.A.barcodes, distinct.jq1.B.barcodes, distinct.jq1.C.barcodes, jq1.all.barcodes,
   jq1.T2.barcodes, jq1.dtp.A.df, jq1.dtp.B.df, jq1.dtp.C.df)
gc()

#Plot data using alluvial
#Generate UMAP including all info in wave plots from DTPs
jq1.dtp.cellular.A.df <- jq1.dtp.cellular.A.df %>% filter(Sample == "Neuro_A")
jq1.dtp.cellular.B.df <- jq1.dtp.cellular.B.df %>% filter(Sample == "Neuro_B")
jq1.dtp.cellular.C.df <- jq1.dtp.cellular.C.df %>% filter(Sample == "Neuro_C")

#Find overlapping clones in each replicate
A <- jq1.dtp.cellular.A.df$real_bc44
B <- jq1.dtp.cellular.B.df$real_bc44
C <- jq1.dtp.cellular.C.df$real_bc44
clone.list <- Reduce(intersect, list(A, B, C))

test <- rbind(jq1.dtp.cellular.A.df, jq1.dtp.cellular.B.df, jq1.dtp.cellular.C.df)
all.jq1.cellular.dtp.df <- left_join(test, full.meta.df, by=c("real_bc44", "Sample"))

#Clean up files
rm(jq1.dtp.cellular.A.df, jq1.dtp.cellular.B.df, jq1.dtp.cellular.C.df, jq1.dtp.cellular.df,
   jq1.dtp.cellular, jq1.dtp.filtered, dtp.jq1.meta.clones.filtered)
gc()

#Generate colour palette for seurat_clusters_0.2
P9 <- c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA")
P9 <- as.vector(t(matrix(P9)))
names(P9) = c("0", "1", "2", "3", "4", "5", "6", "7", "8")

####
#JQ1 AMT states
jq1.amt.dtp.df <- full.meta.df %>%
  filter(Condition == "Untreated" | Condition == "JQ1" | Condition == "JQ1_rec") %>%
  select(c(AMT.state, Sample, real_bc44, Condition, Count))

jq1.amt.dtp.df <- jq1.amt.dtp.df[, c(3, 4, 1, 2, 5)]

tmp <- jq1.amt.dtp.df %>%
  group_by(Sample, Condition, real_bc44) %>%
  mutate(total_clones=sum(Count)) %>%
  ungroup() %>%
  group_by(real_bc44, Condition, Sample) %>%
  mutate(percentage = Count/total_clones * 100) %>%
  ungroup() 

tmp$timepoint[tmp$Condition == "Untreated"] <- "T1"
tmp$timepoint[tmp$Condition == "JQ1"] <- "T2"
tmp$timepoint[tmp$Condition == "JQ1_rec"] <- "T3"

tmp2 <- tmp %>%
  group_by(real_bc44, Condition, Sample) %>%
  summarise(count = n_distinct(AMT.state))

jq1.amt.states <- tmp2 %>%
  group_by(count, Condition, Sample) %>%
  summarise(clones = length(unique(real_bc44)))

jq1.amt.states.summary <- jq1.amt.states %>%
  filter(Condition == "JQ1_rec")

#Calculate percentage of clones with each representation
jq1.amt.states.perc <- jq1.amt.states.summary %>%
  group_by(Sample, Condition) %>%
  mutate(total_freq = sum(clones)) %>%
  ungroup() %>%
  mutate(perc = clones/total_freq * 100) %>%
  ungroup()

jq1.amt.states.perc$count <- as.character(jq1.amt.states.perc$count)

#Plot manuscript extended figure 7c - JQ1
jq1.amt.states.perc %>%
  ggplot( ,aes(x=Sample, y=perc))+
  geom_bar(aes(x=Sample, y=perc, fill=count), colour="black", stat = "identity", position = "stack") +
  xlab("States Occupied")+
  scale_fill_manual(values = P3)+
  ylab("Number of Clones")+
  ggtitle("JQ1 Recovery AMT states Clone Summary")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size =10),
        axis.text.x = element_text(size=10),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))

#Define clones which are only observed in one cellular state within untreated
jq1.amt.states.ss <- tmp2 %>%
  filter(count == 1 & Condition == "Untreated") %>% ungroup() %>%
  select(c(real_bc44, Sample))

jq1.amt.states.ss.df <- merge(jq1.amt.states.ss, tmp, by=c("real_bc44", "Sample"))
length(unique(jq1.amt.states.ss.df$real_bc44)) #154 clones

jq1.amt.states.ss.T2 <- jq1.amt.states.ss.df %>%
  group_by(real_bc44, Sample, timepoint) %>%
  summarise(barcodes = n_distinct(real_bc44))
jq1.amt.states.ss.T2$clone_state <- NA

dtp.jq1.amt.clones <- jq1.amt.states.ss.T2 %>%
  mutate(clone_state = ifelse(barcodes == "1" & timepoint == "T3", "DTP", clone_state))

dtp.jq1.amt.clones.filtered <- dtp.jq1.amt.clones %>%
  filter(clone_state == "DTP") %>%
  select(-c(timepoint, barcodes))

jq1.amt.dtp.df <- merge(dtp.jq1.amt.clones.filtered, tmp, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp2, jq1.amt.states, jq1.amt.states.perc, jq1.amt.states.ss,
   jq1.amt.states.ss.df, jq1.amt.states.ss.T2, jq1.amt.states.summary,
   dtp.jq1.amt.clones, dtp.jq1.amt.clones.filtered)
gc()

#Remove datapoints with no information at T2
jq1.T2.barcodes <- jq1.amt.dtp.df$real_bc44[jq1.amt.dtp.df$Sample == "Neuro_A" & jq1.amt.dtp.df$timepoint == "T2"]
jq1.all.barcodes <- jq1.amt.dtp.df$real_bc44[jq1.amt.dtp.df$Sample == "Neuro_A"]
distinct.jq1.A.barcodes <- c(setdiff(jq1.all.barcodes, jq1.T2.barcodes), setdiff(jq1.T2.barcodes, jq1.all.barcodes))
jq1.dtp.A <- jq1.amt.dtp.df %>% filter(!(real_bc44 %in% distinct.jq1.A.barcodes))

#Visualise for ggalluvial
jq1.dtp.A.df <- jq1.dtp.A %>% select(c(real_bc44, timepoint, Count, AMT.state, Sample))

tmp <- jq1.dtp.A.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = AMT.state) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
jq1.dtp.amt.A.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Remove datapoints with no information at T2
jq1.T2.barcodes <- jq1.amt.dtp.df$real_bc44[jq1.amt.dtp.df$Sample == "Neuro_B" & jq1.amt.dtp.df$timepoint == "T2"]
jq1.all.barcodes <- jq1.amt.dtp.df$real_bc44[jq1.amt.dtp.df$Sample == "Neuro_B"]
distinct.jq1.B.barcodes <- c(setdiff(jq1.all.barcodes, jq1.T2.barcodes), setdiff(jq1.T2.barcodes, jq1.all.barcodes))
jq1.dtp.B <- jq1.amt.dtp.df %>% filter(!(real_bc44 %in% distinct.jq1.B.barcodes))

#Visualise for ggalluvial
jq1.dtp.B.df <- jq1.dtp.B %>% select(c(real_bc44, timepoint, Count, AMT.state, Sample))

tmp <- jq1.dtp.B.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = AMT.state) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
jq1.dtp.amt.B.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Remove datapoints with no information at T2
jq1.T2.barcodes <- jq1.amt.dtp.df$real_bc44[jq1.amt.dtp.df$Sample == "Neuro_C" & jq1.amt.dtp.df$timepoint == "T2"]
jq1.all.barcodes <- jq1.amt.dtp.df$real_bc44[jq1.amt.dtp.df$Sample == "Neuro_C"]
distinct.jq1.C.barcodes <- c(setdiff(jq1.all.barcodes, jq1.T2.barcodes), setdiff(jq1.T2.barcodes, jq1.all.barcodes))
jq1.dtp.C <- jq1.amt.dtp.df %>% filter(!(real_bc44 %in% distinct.jq1.C.barcodes))

#Visualise for ggalluvial
jq1.dtp.C.df <- jq1.dtp.C %>% select(c(real_bc44, timepoint, Count, AMT.state, Sample))

tmp <- jq1.dtp.C.df %>%
  group_by(real_bc44, Sample, Count) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = AMT.state) %>%
  select(-row)

tmp.T1 <- tmp %>%
  select(-c(T2,T3)) %>%
  filter(!(is.na(T1))) %>%
  group_by(real_bc44, Sample, T1) %>%
  summarise(sum_Count = sum(Count))

tmp.T2 <- tmp %>%
  select(-c(T1,T3)) %>%
  filter(!(is.na(T2))) %>%
  group_by(real_bc44, Sample, T2) %>%
  summarise(sum_Count = sum(Count))

tmp.T3 <- tmp %>%
  select(-c(T2,T1)) %>%
  filter(!(is.na(T3))) %>%
  group_by(real_bc44, Sample, T3) %>%
  summarise(sum_Count = sum(Count))

test.merge <- left_join(tmp.T1, tmp.T2, by=c("real_bc44", "Sample"))
jq1.dtp.amt.C.df <- left_join(test.merge, tmp.T3, by=c("real_bc44", "Sample"))

#Clean up files
rm(tmp, tmp.T1, tmp.T2, tmp.T3, test.merge, jq1.dtp.A, jq1.dtp.B, jq1.dtp.C,
   distinct.jq1.A.barcodes, distinct.jq1.B.barcodes, distinct.jq1.C.barcodes, jq1.all.barcodes,
   jq1.T2.barcodes, jq1.dtp.A.df, jq1.dtp.B.df, jq1.dtp.C.df)
gc()

#Plot data using alluvial
#Generate UMAP including all info in wave plots from DTPs
jq1.dtp.amt.A.df <- jq1.dtp.amt.A.df %>% filter(Sample == "Neuro_A")
jq1.dtp.amt.B.df <- jq1.dtp.amt.B.df %>% filter(Sample == "Neuro_B")
jq1.dtp.amt.C.df <- jq1.dtp.amt.C.df %>% filter(Sample == "Neuro_C")

#Find overlapping clones in each replicate
A <- jq1.dtp.amt.A.df$real_bc44
B <- jq1.dtp.amt.B.df$real_bc44
C <- jq1.dtp.amt.C.df$real_bc44
clone.list <- Reduce(intersect, list(A, B, C))

test <- rbind(jq1.dtp.amt.A.df, jq1.dtp.amt.B.df, jq1.dtp.amt.C.df)
all.jq1.amt.dtp.df <- left_join(test, full.meta.df, by=c("real_bc44", "Sample"))

#Clean up files
rm(jq1.dtp.amt.A.df, jq1.dtp.amt.B.df, jq1.dtp.amt.C.df, jq1.dtp.amt.df,
   jq1.dtp.amt, jq1.dtp.filtered, dtp.jq1.meta.clones.filtered)
gc()

#To counteract the overlap of cloens between replicates
#Make new variable which is barcode + replicate
all.jq1.amt.dtp.df <- all.jq1.amt.dtp.df %>%
  mutate(barcode_rep = paste0(real_bc44, "_", Sample))

all.jq1.cellular.dtp.df <- all.jq1.cellular.dtp.df %>%
  mutate(barcode_rep = paste0(real_bc44, "_", Sample))

#Visulaise all informaito for cellular clusters and AMT states
#Plot manuscript figure 4e - left hand panel
all.jq1.amt.dtp.df %>%
  filter(Condition == "Untreated") %>%
  mutate(Condition = factor(Condition, levels=c("Untreated"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "white", pointsize = 2)+
  geom_scattermore(aes(colour = AMT.state), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(breaks = c("MES", "intermediate", "ADRN"),
                      values=c("#F37735", "lightgrey", "#990099")) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Plot manuscript figure 4d - left hand panel
all.jq1.cellular.dtp.df %>%
  filter(Condition == "Untreated") %>%
  mutate(Condition = factor(Condition, levels=c("Untreated"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "white", pointsize = 2)+
  geom_scattermore(aes(colour = seurat_clusters.0.2), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(values=P9) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Plot manuscript figure 4e - right hand panel
all.jq1.amt.dtp.df %>%
  filter(Condition == "JQ1_rec") %>%
  mutate(Condition = factor(Condition, levels=c("JQ1_rec"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "white", pointsize = 2)+
  geom_scattermore(aes(colour = AMT.state), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(breaks = c("MES", "intermediate", "ADRN"),
                      values=c("#F37735", "lightgrey", "#990099")) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Plot manuscript figure 4d - left hand panel
all.jq1.cellular.dtp.df %>%
  filter(Condition == "JQ1_rec") %>%
  mutate(Condition = factor(Condition, levels=c("JQ1_rec"))) %>%
  ggplot(aes(x = UMAP_1, 
             y = UMAP_2))+
  geom_scattermore(data = full.meta.df[, c("UMAP_1", "UMAP_2")], colour = "white", pointsize = 2)+
  geom_scattermore(aes(colour = seurat_clusters.0.2), pointsize = 5)+
  guides(color = guide_legend(override.aes = list(size = 6)), size=guide_legend(ncol=6)) +
  scale_colour_manual(values=P9) +
  labs(x = "UMAP_1",
       y = "UMAP_2",
       color="Experimental Condition") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.title = element_text(size = 18),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size = 8, hjust=0.95, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.border = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Make compiled list of clones to generate JQ1 colour palette for barcodes
amt.list <- unique(all.jq1.amt.dtp.df$barcode_rep)
cellular.list <- unique(all.jq1.cellular.dtp.df$barcode_rep)
jq1.dtp.barcodes <- c(amt.list, cellular.list)
unique.jq1.dtp.barcodes <- unique(jq1.dtp.barcodes)

#Generate barcode palette
set.seed(12)
jq1_palette = createPalette(600, c("#ff0000", "#00ff00", "#0000ff"))
jq1_palette <- jq1_palette[!(jq1_palette %in% cisplatin_palette)]
jq1_palette <- jq1_palette[!(jq1_palette %in% ut_palette)]
jq1_palette <- sample(jq1_palette, 149, replace=FALSE)
jq1_palette <- as.vector(t(matrix(jq1_palette)))
names(jq1_palette) = unique(as.character(jq1.dtp.barcodes))

#Plot manuscript extended figure 7e - JQ1
all.jq1.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill=T1)) +
  geom_stratum( width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manuscript extended figure 7e - JQ1
all.jq1.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T1), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manuscript extended figure 7e - JQ1
all.jq1.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T2), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manuscript extended figure 7e - JQ1
all.jq1.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=real_bc44, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T3), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
                    values=c("#CC6677", "#228833", "#0077BB", "#332288", "#66CCEE", "#009988", "#EECC66", "#CC3311", "#CCDDAA"),
                    na.value="black",
                    name = "MetaCluster at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manuscript extended figure 7e - JQ1
all.jq1.cellular.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("0", "1", "2", "3", "4", "5", "6", "7", "8", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(width = 1/8, colour="black", aes(fill=barcode_rep), linewidth=0) +
  scale_fill_manual(values = jq1_palette) +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank(),
        legend.position = "none")

#AMT
#Plot manuscript extended figure 7f - JQ1
all.jq1.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill=T1)) +
  geom_stratum( width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manuscript extended figure 7f - JQ1
all.jq1.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T1), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manuscript extended figure 7f - JQ1
all.jq1.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T2), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manuscript extended figure 7f - JQ1
all.jq1.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=real_bc44, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=T3), width = 1/8, alpha=0.8, colour="black") +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "lightgrey", "#990099"),
                    na.value="black",
                    name = "AMT State at T1") +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank())

#Plot manuscript extended figure 7f - JQ1
all.jq1.amt.dtp.df %>%
  mutate(T1 = factor(T1, levels=c("MES", "intermediate", "ADRN", NA))) %>%
  ggplot( aes(axis1=barcode_rep, axis2 = T1, axis3=T2, axis4=T3, y = log10(sum_Count.x))) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"Untreated", "JQ1", "JQ1 Recovery"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(width = 1/8, colour="black", aes(fill=barcode_rep), linewidth=0) +
  scale_fill_manual(values = jq1_palette) +
  theme_cowplot() +
  ylab("log10(Clone Frequency)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        axis.line = element_blank(),
        legend.position = "none")


#Load in UT data sets for AMT and cellular clusters
all.ut.amt.dtp.df <- read.csv("datafiles/all_ut_amt_dtp.csv")
all.cisplatin.amt.dtp.df <- read.csv("datafiles/all_cis_amt_dtp.csv")
all.jq1.amt.dtp.df <- read.csv("datafiles/all_jq1_amt_dtp.csv")

#Define levels of plotting for AMT states
all.ut.amt.dtp.df$T1 <- factor(all.ut.amt.dtp.df$T1, levels = c("ADRN", "intermediate", "MES"))
all.ut.amt.dtp.df$T2 <- factor(all.ut.amt.dtp.df$T2, levels = c("ADRN", "intermediate", "MES"))

#Untreated transition plots
ut.transitions <- all.ut.amt.dtp.df %>%
  select(c(barcode_rep, T1, T2)) %>%
  ungroup() %>%
  mutate(transition = ifelse(T1 == "MES" & T2 == "MES", "MES fixed",
                             ifelse(T1 == "ADRN" & T2 == "ADRN", "ADRN fixed",
                                    ifelse(T1 == "intermediate" & T2 == "intermediate", "Intermediate fixed",
                                           ifelse(T1 == "MES" & T2 == "ADRN", "MAT",
                                                  ifelse(T1 == "ADRN" & T2 == "MES", "AMT",
                                                         ifelse(T1 == "intermediate" & T2 == "ADRN", "IAT",
                                                                ifelse(T1 == "intermediate" & T2 == "MES", "IMT",
                                                                       ifelse(T1 == "ADRN" & T2 == "intermediate", "AIT",
                                                                              ifelse(T1 == "MES" & T2 == "intermediate", "MIT", NA))))))))))

#Check transitions present
unique(ut.transitions$transition)
#We are seeing everything we should so lets summarise how many unique clones we see in untreated conditions
#Now lets visualise in an alluvial these transitions
transitions <- c("MES fixed", "Intermediate fixed", "ADRN fixed", "AIT", "AMT", "IAT", "IMT", "MAT", "MIT")
colours <- c("#f9bb9a", "#D3D3D3", "#c266c2", "#730073", "#3d003d", "#898989", "#0b0b0b", "#cf652d", "#6d3618")
names(colours) <- transitions

#Visualise via ggplot
ut.transitions.summary <- ut.transitions %>%
  group_by(transition) %>%
  summarise(clones = n_distinct(barcode_rep))

#Plot manuscript extended figure 7g - untreated
ut.transitions.summary %>%
  mutate(transition = fct_reorder(transition, clones)) %>%
  ggplot(aes(x=transition, y=clones, fill=transition)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=seq(0,125,25), limits=c(0,125)) +
  ylab("Number of Clones per Transition Group") +
  xlab("Phenotypic Transition \n (Ordered by Increasing Clone Number)") +
  scale_x_discrete(labels=c("HYB to MES", "HYB to HYB", "HYB to ADRN", "ADRN to MES", "MES to HYB", "MES to ADRN", "ADRN to HYB", "ADRN to ADRN", "MES to MES")) +
  scale_fill_manual(values=colours,
                    breaks=c("IMT", "Intermediate fixed", "IAT", "AMT", "MIT", "MAT", "AIT", "ADRN fixed", "MES fixed"),
                    labels=c("HYB to MES", "HYB to HYB", "HYB to ADRN", "ADRN to MES", "MES to HYB", "MES to ADRN", "ADRN to HYB", "ADRN to ADRN", "MES to MES"),
                    name = "Phenotypic Transition") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1))


#Generate palette for transitions
ut.sum <- ut.transitions %>%
  distinct()

#Plot manuscript figure 4f - untreated
ut.sum %>%
  ggplot( aes(axis1=T1, axis2=T2, fill=transition), alpha=1) +
  geom_alluvium(width=1/8) +
  scale_x_discrete(limits = c("Untreated", "Untreated Recovery"), expand = c(.2, .05)) +
  geom_stratum(width = 1/8, colour="black", aes(fill=transition), linewidth=0) +
  scale_fill_manual(values=colours,
                    labels=c("ADRN to ADRN", "ADRN to HYB", "ADRN to MES", "HYB to ADRN", "HYB to HYB", "HYB to MES", "MES to ADRN", "MES to HYB", "MES to MES"),
                    breaks=c("ADRN fixed", "AIT", "AMT", "IAT", "Intermediate fixed", "IMT", "MAT", "MIT", "MES fixed"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,250,50), limits=c(0,250))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4f - untreated
ut.sum %>%
  ggplot( aes(axis1=T1, axis2=T2, fill=T1), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("Untreated", "Untreated Recovery"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,250,50), limits=c(0,250))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4f - untreated
ut.sum %>%
  ggplot( aes(axis1=T1, axis2=T2, fill=T2), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("Untreated", "Untreated Recovery"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,250,50), limits=c(0,250))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Define levels of plotting for AMT states
all.cisplatin.amt.dtp.df$T1 <- factor(all.cisplatin.amt.dtp.df$T1, levels = c("ADRN", "intermediate", "MES"))
all.cisplatin.amt.dtp.df$T2 <- factor(all.cisplatin.amt.dtp.df$T2, levels = c("ADRN", "intermediate", "MES"))

#Cisplatin ON transition plots
cis.on.transitions <- all.cisplatin.amt.dtp.df %>%
  select(c(barcode_rep, T1, T2)) %>%
  ungroup() %>%
  mutate(transition = ifelse(T1 == "MES" & T2 == "MES", "MES fixed",
                             ifelse(T1 == "ADRN" & T2 == "ADRN", "ADRN fixed",
                                    ifelse(T1 == "intermediate" & T2 == "intermediate", "Intermediate fixed",
                                           ifelse(T1 == "MES" & T2 == "ADRN", "MAT",
                                                  ifelse(T1 == "ADRN" & T2 == "MES", "AMT",
                                                         ifelse(T1 == "intermediate" & T2 == "ADRN", "IAT",
                                                                ifelse(T1 == "intermediate" & T2 == "MES", "IMT",
                                                                       ifelse(T1 == "ADRN" & T2 == "intermediate", "AIT",
                                                                              ifelse(T1 == "MES" & T2 == "intermediate", "MIT", NA))))))))))

#Check transitions present
unique(cis.on.transitions$transition)
#We are seeing everything we should so lets summarise how many unique clones we see in untreated conditions

cis.on.transitions.summary <- cis.on.transitions %>%
  group_by(transition) %>%
  summarise(clones = n_distinct(barcode_rep))

#Plot manuscript extended figure 7g - cisplatin on
cis.on.transitions.summary %>%
  mutate(transition = fct_reorder(transition, clones)) %>%
  ggplot(aes(x=transition, y=clones, fill=transition)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=seq(0,125,25), limits=c(0,125)) +
  ylab("Number of Clones per Transition Group") +
  xlab("Phenotypic Transition \n (Ordered by Increasing Clone Number)") +
  scale_x_discrete(labels=c("MES to ADRN", "ADRN to MES", "ADRN to ADRN", "MES to HYB", "MES to MES")) +
  scale_fill_manual(values=colours,
                    breaks=c("MAT", "AMT", "ADRN fixed", "MIT", "MES fixed"),
                    labels=c("MES to ADRN", "ADRN to MES", "ADRN to ADRN", "MES to HYB", "MES to MES"),
                    name = "Phenotypic Transition") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1))

#Now lets visualise in an alluvial these transitions
#Generate palette for transitions
cis.on.sum <- cis.on.transitions %>%
  distinct()

#Plot manuscript figure 4g - cisplatin on
cis.on.sum %>%
  ggplot( aes(axis1=T1, axis2=T2, fill=transition)) +
  geom_alluvium(width=1/8) +
  scale_x_discrete(limits = c("Untreated", "Cisplatin"), expand = c(.2, .05)) +
  geom_stratum(width = 1/8, colour="black", aes(fill=transition), linewidth=0) +
  scale_fill_manual(values=colours,
                    labels=c("ADRN to ADRN", "ADRN to HYB", "ADRN to MES", "HYB to ADRN", "HYB to HYB", "HYB to MES", "MES to ADRN", "MES to HYB", "MES to MES"),
                    breaks=c("ADRN fixed", "AIT", "AMT", "IAT", "Intermediate fixed", "IMT", "MAT", "MIT", "MES fixed"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,90,10), limits=c(0,90))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4g - cisplatin on
cis.on.sum %>%
  ggplot( aes(axis1=T1, axis2=T2, fill=T1), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("Untreated", "Cisplatin"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,90,10), limits=c(0,90))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4g - cisplatin on
cis.on.sum %>%
  ggplot( aes(axis1=T1, axis2=T2, fill=T2), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("Untreated", "Cisplatin"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,90,10), limits=c(0,90))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Define levels of plotting for AMT states
all.cisplatin.amt.dtp.df$T3 <- factor(all.cisplatin.amt.dtp.df$T3, levels = c("ADRN", "intermediate", "MES"))

#Cisplatin OFF transition plots
cis.off.transitions <- all.cisplatin.amt.dtp.df %>%
  select(c(barcode_rep, T2, T3)) %>%
  ungroup() %>%
  mutate(transition = ifelse(T2 == "MES" & T3 == "MES", "MES fixed",
                             ifelse(T2 == "ADRN" & T3 == "ADRN", "ADRN fixed",
                                    ifelse(T2 == "intermediate" & T3 == "intermediate", "Intermediate fixed",
                                           ifelse(T2 == "MES" & T3 == "ADRN", "MAT",
                                                  ifelse(T2 == "ADRN" & T3 == "MES", "AMT",
                                                         ifelse(T2 == "intermediate" & T3 == "ADRN", "IAT",
                                                                ifelse(T2 == "intermediate" & T3 == "MES", "IMT",
                                                                       ifelse(T2 == "ADRN" & T3 == "intermediate", "AIT",
                                                                              ifelse(T2 == "MES" & T3 == "intermediate", "MIT", NA))))))))))

#Check transitions present
unique(cis.off.transitions$transition)
#We are seeing everything we should so lets summarise how many unique clones we see in untreated conditions

cis.off.transitions.summary <- cis.off.transitions %>%
  group_by(transition) %>%
  summarise(clones = n_distinct(barcode_rep))

#Plot manuscript extended figure 7g - cisplatin off
cis.off.transitions.summary %>%
  mutate(transition = fct_reorder(transition, clones)) %>%
  ggplot(aes(x=transition, y=clones, fill=transition)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=seq(0,125,25), limits=c(0,125)) +
  ylab("Number of Clones per Transition Group") +
  xlab("Phenotypic Transition \n (Ordered by Increasing Clone Number)") +
  scale_x_discrete(labels=c("ADRN to ADRN", "ADRN to HYB", "MES to ADRN", "MES to HYB", "ADRN to MES", "HYB to MES", "MES to MES")) +
  scale_fill_manual(values=colours,
                    breaks=c("ADRN fixed", "AIT", "MAT", "MIT", "AMT", "IMT", "MES fixed"),
                    labels=c("ADRN to ADRN", "ADRN to HYB", "MES to ADRN", "MES to HYB", "ADRN to MES", "HYB to MES", "MES to MES"),
                    name = "Phenotypic Transition") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1))

#Now lets visualise in an alluvial these transitions
#Generate palette for transitions
cis.off.sum <- cis.off.transitions %>%
  distinct()

#Plot manuscript figure 4g - cisplatin off
cis.off.sum %>%
  ggplot( aes(axis1=T2, axis2=T3, fill=transition)) +
  geom_alluvium(width=1/8) +
  scale_x_discrete(limits = c("Cisplatin", "Cisplatin_rec"), expand = c(.2, .05)) +
  geom_stratum(width = 1/8, colour="black", aes(fill=transition), linewidth=0) +
  scale_fill_manual(values=colours,
                    labels=c("ADRN to ADRN", "ADRN to HYB", "ADRN to MES", "HYB to ADRN", "HYB to HYB", "HYB to MES", "MES to ADRN", "MES to HYB", "MES to MES"),
                    breaks=c("ADRN fixed", "AIT", "AMT", "IAT", "Intermediate fixed", "IMT", "MAT", "MIT", "MES fixed"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,90,10), limits=c(0,90))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4g - cisplatin off
cis.off.sum %>%
  ggplot( aes(axis1=T2, axis2=T3, fill=T2), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("Cisplatin", "Cisplatin_rec"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,90,10), limits=c(0,90))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4g - cisplatin off
cis.off.sum %>%
  ggplot( aes(axis1=T2, axis2=T3, fill=T3), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("Cisplatin", "Cisplatin_rec"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,90,10), limits=c(0,90))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Define levels of plotting for AMT states
all.jq1.amt.dtp.df$T1 <- factor(all.jq1.amt.dtp.df$T1, levels = c("ADRN", "intermediate", "MES"))
all.jq1.amt.dtp.df$T2 <- factor(all.jq1.amt.dtp.df$T2, levels = c("ADRN", "intermediate", "MES"))

#JQ1 ON transition plots
jq1.on.transitions <- all.jq1.amt.dtp.df %>%
  select(c(barcode_rep, T1, T2)) %>%
  ungroup() %>%
  mutate(transition = ifelse(T1 == "MES" & T2 == "MES", "MES fixed",
                             ifelse(T1 == "ADRN" & T2 == "ADRN", "ADRN fixed",
                                    ifelse(T1 == "intermediate" & T2 == "intermediate", "Intermediate fixed",
                                           ifelse(T1 == "MES" & T2 == "ADRN", "MAT",
                                                  ifelse(T1 == "ADRN" & T2 == "MES", "AMT",
                                                         ifelse(T1 == "intermediate" & T2 == "ADRN", "IAT",
                                                                ifelse(T1 == "intermediate" & T2 == "MES", "IMT",
                                                                       ifelse(T1 == "ADRN" & T2 == "intermediate", "AIT",
                                                                              ifelse(T1 == "MES" & T2 == "intermediate", "MIT", NA))))))))))

#Check transitions present
unique(jq1.on.transitions$transition)
#We are seeing everything we should so lets summarise how many unique clones we see in untreated conditions

jq1.on.transitions.summary <- jq1.on.transitions %>%
  group_by(transition) %>%
  summarise(clones = n_distinct(barcode_rep))

#Plot manuscript extended figure 7g - jq1 on
jq1.on.transitions.summary %>%
  mutate(transition = fct_reorder(transition, clones)) %>%
  ggplot(aes(x=transition, y=clones, fill=transition)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=seq(0,125,25), limits=c(0,125)) +
  ylab("Number of Clones per Transition Group") +
  xlab("Phenotypic Transition \n (Ordered by Increasing Clone Number)") +
  scale_x_discrete(labels=c("HYB to ADRN", "HYB to HYB", "MES to ADRN", "MES to HYB", "ADRN to MES", "ADRN to HYB",  "ADRN to ADRN", "MES to MES")) +
  scale_fill_manual(values=colours,
                    breaks=c("IAT", "Intermediate fixed", "MAT", "MIT", "AMT", "AIT", "ADRN fixed", "MES fixed"),
                    labels=c("HYB to ADRN", "HYB to HYB", "MES to ADRN", "MES to HYB", "ADRN to MES", "ADRN to HYB",  "ADRN to ADRN", "MES to MES"),
                    name = "Phenotypic Transition") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1))

#Now lets visualise in an alluvial these transitions
#Generate palette for transitions
jq1.on.sum <- jq1.on.transitions %>%
  distinct()

#Plot manuscript figure 4h - jq1 on
jq1.on.sum %>%
  ggplot( aes(axis1=T1, axis2=T2, fill=transition)) +
  geom_alluvium() +
  scale_x_discrete(limits = c("Untreated", "JQ1"), expand = c(.2, .05)) +
  geom_stratum(width = 1/8, colour="black", aes(fill=transition), linewidth=0) +
  scale_fill_manual(values=colours,
                    labels=c("ADRN to ADRN", "ADRN to HYB", "ADRN to MES", "HYB to ADRN", "HYB to HYB", "HYB to MES", "MES to ADRN", "MES to HYB", "MES to MES"),
                    breaks=c("ADRN fixed", "AIT", "AMT", "IAT", "Intermediate fixed", "IMT", "MAT", "MIT", "MES fixed"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,150,25), limits=c(0,150))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4h - jq1 on
jq1.on.sum %>%
  ggplot( aes(axis1=T1, axis2=T2, fill=T1), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("Untreated", "JQ1"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,150,25), limits=c(0,150))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4h - jq1 on
jq1.on.sum %>%
  ggplot( aes(axis1=T1, axis2=T2, fill=T2), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("Untreated", "JQ1"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,150,25), limits=c(0,150))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Define levels of plotting for AMT states
all.jq1.amt.dtp.df$T3 <- factor(all.jq1.amt.dtp.df$T3, levels = c("ADRN", "intermediate", "MES"))

#Cisplatin OFF transition plots
jq1.off.transitions <- all.jq1.amt.dtp.df %>%
  select(c(barcode_rep, T2, T3)) %>%
  ungroup() %>%
  mutate(transition = ifelse(T2 == "MES" & T3 == "MES", "MES fixed",
                             ifelse(T2 == "ADRN" & T3 == "ADRN", "ADRN fixed",
                                    ifelse(T2 == "intermediate" & T3 == "intermediate", "Intermediate fixed",
                                           ifelse(T2 == "MES" & T3 == "ADRN", "MAT",
                                                  ifelse(T2 == "ADRN" & T3 == "MES", "AMT",
                                                         ifelse(T2 == "intermediate" & T3 == "ADRN", "IAT",
                                                                ifelse(T2 == "intermediate" & T3 == "MES", "IMT",
                                                                       ifelse(T2 == "ADRN" & T3 == "intermediate", "AIT",
                                                                              ifelse(T2 == "MES" & T3 == "intermediate", "MIT", NA))))))))))

#Check transitions present
unique(jq1.off.transitions$transition)
#We are seeing everything we should so lets summarise how many unique clones we see in untreated conditions

jq1.off.transitions.summary <- jq1.off.transitions %>%
  group_by(transition) %>%
  summarise(clones = n_distinct(barcode_rep))

#Plot manuscript extended figure 7g - jq1 off
jq1.off.transitions.summary %>%
  mutate(transition = fct_reorder(transition, clones)) %>%
  ggplot(aes(x=transition, y=clones, fill=transition)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=seq(0,125,25), limits=c(0,125)) +
  ylab("Number of Clones per Transition Group") +
  xlab("Phenotypic Transition \n (Ordered by Increasing Clone Number)") +
  scale_x_discrete(labels=c("HYB to HYB", "ADRN to HYB", "MES to HYB", "HYB to MES", "ADRN to MES", "HYB to ADRN",  "ADRN to ADRN", "MES to ADRN", "MES to MES")) +
  scale_fill_manual(values=colours,
                    breaks=c("Intermediate fixed", "AIT", "MIT", "IMT", "AMT", "IAT", "ADRN fixed", "MAT", "MES fixed"),
                    labels=c("HYB to HYB", "ADRN to HYB", "MES to HYB", "HYB to MES", "ADRN to MES", "HYB to ADRN",  "ADRN to ADRN", "MES to ADRN", "MES to MES"),
                    name = "Phenotypic Transition") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1))

#Now lets visualise in an alluvial these transitions
#Generate palette for transitions
jq1.off.sum <- jq1.off.transitions %>%
  distinct()

#Plot manuscript figure 4h - jq1 off
jq1.off.sum %>%
  ggplot( aes(axis1=T2, axis2=T3, fill=transition)) +
  geom_alluvium(width=1/8) +
  scale_x_discrete(limits = c("JQ1", "JQ1_rec"), expand = c(.2, .05)) +
  geom_stratum(width = 1/8, colour="black", aes(fill=transition), linewidth=0) +
  scale_fill_manual(values=colours,
                    labels=c("ADRN to ADRN", "ADRN to HYB", "ADRN to MES", "HYB to ADRN", "HYB to HYB", "HYB to MES", "MES to ADRN", "MES to HYB", "MES to MES"),
                    breaks=c("ADRN fixed", "AIT", "AMT", "IAT", "Intermediate fixed", "IMT", "MAT", "MIT", "MES fixed"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,200,25), limits=c(0,200))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4h - jq1 off
jq1.off.sum %>%
  ggplot( aes(axis1=T2, axis2=T3, fill=T2), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("JQ1", "JQ1_rec"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,200,25), limits=c(0,200))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")

#Plot manuscript figure 4h - jq1 off
jq1.off.sum %>%
  ggplot( aes(axis1=T2, axis2=T3, fill=T3), alpha=1) +
  geom_alluvium(width = 1/8, colour="black", linewidth=0) +
  scale_x_discrete(limits = c("JQ1", "JQ1_rec"), expand = c(.2, .05)) +
  geom_stratum(width=1/8) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic Transition") +
  scale_y_continuous(breaks=seq(0,200,25), limits=c(0,200))+
  theme_cowplot() +
  xlab("") +
  ylab("Unique Clones per Transition Group")


#####
#Now we need to generate bar plots summarising the number of DTPs in each cell state to accompany the UMAPS
#Load in UT data sets for AMT and cellular clusters
all.cisplatin.amt.dtp.df <- read.csv("datafiles/all_cis_amt_dtp.csv")
all.jq1.amt.dtp.df <- read.csv("datafiles/all_jq1_amt_dtp.csv")

#Cisplatin DTPS
cisplatin.rec <- all.cisplatin.amt.dtp.df %>%
  select(c(barcode_rep, T1, T2, T3, Condition)) %>%
  filter(Condition == "Cisplatin_rec" | Condition == "Untreated")

cisplatin.rec.sum <- cisplatin.rec %>%
  filter(Condition == "Cisplatin_rec") %>%
  group_by(T3, Condition) %>%
  summarise(count = n_distinct(barcode_rep))
colnames(cisplatin.rec.sum)[which(names(cisplatin.rec.sum) == "T3")] <- "state"

untreated.cis.sum <- cisplatin.rec %>%
  filter(Condition == "Untreated") %>%
  group_by(T1, Condition) %>%
  summarise(count = n_distinct(barcode_rep))
colnames(untreated.cis.sum)[which(names(untreated.cis.sum) == "T1")] <- "state"

ut.cis.rec <- rbind(cisplatin.rec.sum, untreated.cis.sum)

ut.cis.rec$Condition <- factor(ut.cis.rec$Condition, levels = c("Untreated", "Cisplatin_rec"))

ut.cis.rec.sum <- ut.cis.rec %>%
  group_by(Condition) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  group_by(Condition, state) %>%
  mutate(proportion = (count/total_count))

#Plot manuscript figure 4c - cisplatin
ut.cis.rec.sum %>%
  ggplot(aes(x=Condition, y=proportion, fill=state))+
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(labels=c("Untreated", "Cisplatin OFF")) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic State") +
  xlab("Proportion of AMT state") +
  ylab("") +
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#JQ1 DTPS
jq1.rec <- all.jq1.amt.dtp.df %>%
  select(c(barcode_rep, T1, T2, T3, Condition)) %>%
  filter(Condition == "JQ1_rec" | Condition == "Untreated")

jq1.rec.sum <- jq1.rec %>%
  filter(Condition == "JQ1_rec") %>%
  group_by(T3, Condition) %>%
  summarise(count = n_distinct(barcode_rep))
colnames(jq1.rec.sum)[which(names(jq1.rec.sum) == "T3")] <- "state"

untreated.jq1.sum <- jq1.rec %>%
  filter(Condition == "Untreated") %>%
  group_by(T1, Condition) %>%
  summarise(count = n_distinct(barcode_rep))
colnames(untreated.jq1.sum)[which(names(untreated.jq1.sum) == "T1")] <- "state"

ut.jq1.rec <- rbind(jq1.rec.sum, untreated.jq1.sum)

ut.jq1.rec$Condition <- factor(ut.jq1.rec$Condition, levels = c("Untreated", "JQ1_rec"))

ut.jq1.rec.sum <- ut.jq1.rec %>%
  group_by(Condition) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  group_by(Condition, state) %>%
  mutate(proportion = (count/total_count))

#Plot manuscript figure 4e - JQ1
ut.jq1.rec.sum %>%
  ggplot(aes(x=Condition, y=proportion, fill=state))+
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(labels=c("Untreated", "JQ1 OFF")) +
  scale_fill_manual(breaks = c("MES", "intermediate", "ADRN"),
                    values=c("#F37735", "#D3D3D3", "#990099"),
                    labels=c("MES", "HYB", "ADRN"),
                    name = "Phenotypic State") +
  xlab("Proportion of AMT state") +
  ylab("") +
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

