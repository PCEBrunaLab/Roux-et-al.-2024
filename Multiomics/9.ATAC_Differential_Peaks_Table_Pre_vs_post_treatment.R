#Differential Peaks Table Generation

#Initial R Env with correct installation of packages - ArchR_Renv
renv::activate()
renv::load()

#Set Working Directory and Arch R parameters
setwd("~/OneDrive/PhD/Thesis/Thesis Plots/Results Chapter 3/Sub Chapter 1/")

addArchRThreads(threads = 1) 
addArchRGenome("hg38")

# define thresholds for differential vs marker peaks study 
threshold_loose_logFC = log2(2)
threshold_loose_FDR = 0.1

#Define plot variables for condition
#BiocManager::install("paletteer")
library(paletteer)
group.cols <- paletteer_d("rcartocolor::Prism")
names(group.cols) <- c("Clone2", "Clone3", "Clone4", "Clone5", "Clone7",
                       "Clone9", "Clone11", "Clone13")

## Load in Data Files ----
project <- readRDS("ArchR_projects/project_clones_working/Save-ArchR-Project.rds")

#Save ArchR project version for Pre Vs Post Treatment Analysis
saveArchRProject(ArchRProj = project, outputDirectory = "project_clones_working/", load = TRUE)
atac_metadata <- getCellColData(project)

#Load RNA subtype annotations
#Load RNA subtype annotations
seurat_sample_all <- readRDS("datafiles/nb_seurat_AMT.rds")

unique(seurat_sample_all$Sample)
seurat_sample_all <- RenameCells(seurat_sample_all, new.names = gsub("Clone-2_", "Clone_2#", colnames(seurat_sample_all)))
seurat_sample_all <- RenameCells(seurat_sample_all, new.names = gsub("Clone-3_", "Clone_3#", colnames(seurat_sample_all)))
seurat_sample_all <- RenameCells(seurat_sample_all, new.names = gsub("Clone-4_", "Clone_4#", colnames(seurat_sample_all)))
seurat_sample_all <- RenameCells(seurat_sample_all, new.names = gsub("Clone-5_", "Clone_5#", colnames(seurat_sample_all)))
seurat_sample_all <- RenameCells(seurat_sample_all, new.names = gsub("Clone-7_", "Clone_7#", colnames(seurat_sample_all)))
seurat_sample_all <- RenameCells(seurat_sample_all, new.names = gsub("Clone-9_", "Clone_9#", colnames(seurat_sample_all)))
seurat_sample_all <- RenameCells(seurat_sample_all, new.names = gsub("Clone-11_", "Clone_11#", colnames(seurat_sample_all)))
seurat_sample_all <- RenameCells(seurat_sample_all, new.names = gsub("Clone-13_", "Clone_13#", colnames(seurat_sample_all)))

#Assign plasticity values
seurat_sample_all$Plasticity <- NA
seurat_sample_all$Plasticity[seurat_sample_all$Sample == "Clone 2" | seurat_sample_all$Sample == "Clone 3" | seurat_sample_all$Sample == "Clone 4"] <- "Fixed"
seurat_sample_all$Plasticity[seurat_sample_all$Sample == "Clone 5" | seurat_sample_all$Sample == "Clone 7" | seurat_sample_all$Sample == "Clone 11" | seurat_sample_all$Sample == "Clone 13"] <- "Plastic"
seurat_sample_all$Plasticity[seurat_sample_all$Sample == "Clone 9"] <- "Parental"

saveRDS(seurat_sample_all, "datafiles/nb_seurat_AMT_updated.rds")

#Load in seurat object
seurat_sample_all <- readRDS("datafiles/nb_seurat_AMT_updated.rds")

rna_metadata <- seurat_sample_all@meta.data

#We need to match ArchR project rownames/CellIDs and the Seurat object rownames/CellIDS
#Examine each first 
cells_archr <- rownames(atac_metadata) #57805 cells
cells_seurat <- rownames(rna_metadata) #59261 cells

#Check intersect between two sets of cells
cells_archr <- rownames(atac_metadata)
cells_seurat <- rownames(rna_metadata)
intersect_cells <- intersect(cells_seurat, cells_archr) #47851 cells overlap between the datasets 
#i.e. 82.8% cells with ATAC have RNA
#i.e. 80.7% cells with RNA have ATAC

#Add RNA subtypes info to ATAC ArchR project
atac_metadata$AMT.state <- NA
atac_metadata[intersect(rownames(atac_metadata), rownames(rna_metadata)), "AMT.state"] <- as.character(rna_metadata[intersect(rownames(atac_metadata), rownames(rna_metadata)),"AMT.state"]) # 10939 common cells between ArchR (all cells) and Seurat (tumor cells)
project$AMT.state <- atac_metadata$AMT.state # add to ArchR project

# Visualise AMT score across UMAP Harmony ATAC
pdf("plots/multiomics_AMT_state_on_ATAC_UMAP.pdf")
plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "AMT.state", embedding = "UMAP_LSI")
dev.off()

## Compute Differential Analysis ----
#We will only take cells which have intersect of ATAC and RNA information not just ATAC alone
#First lets generate a variable with AMT state and condition to isolate properly
project$Condition_AMT <- paste0(project$Sample, "_", project$AMT.state)
atac_metadata <- getCellColData(project)
unique(atac_metadata$Condition_AMT)

#Find number of cells in each group
table(atac_metadata$Condition_AMT)

#Add plasticity subtypes info to ATAC ArchR project
atac_metadata$Plasticity <- NA
atac_metadata[intersect(rownames(atac_metadata), rownames(rna_metadata)), "Plasticity"] <- as.character(rna_metadata[intersect(rownames(atac_metadata), rownames(rna_metadata)),"Plasticity"]) # 10939 common cells between ArchR (all cells) and Seurat (tumor cells)
project$Plasticity <- atac_metadata$Plasticity # add to ArchR project

#Lets also generated a Plasticity and state variable, taking whether the clone is plastic rather than the clone itself
project$Plasticity_AMT <- paste0(project$Plasticity, "_", project$AMT.state)
atac_metadata <- getCellColData(project)
unique(atac_metadata$Plasticity_AMT)

#Find number of cells in each group
table(atac_metadata$Plasticity_AMT)

#Fixed_ADRN    Fixed_intermediate             Fixed_MES                 NA_NA         Parental_ADRN Parental_intermediate 
#458                  7715                 11561                  9234                  2301                  2605 
#Parental_MES          Plastic_ADRN  Plastic_intermediate           Plastic_MES 
#2048                  5858                 11708                  3597                   

#Conduct MES Plastic Vs MES Fixed
diffPeaks_MES_plastic_fixed <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "Plasticity_AMT",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "Plastic_MES",
  bgdGroups = "Fixed_MES",
  maxCells = 12000 # choose high number of max cells to make sure all cells in each condition are included
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

saveRDS(diffPeaks_MES_plastic_fixed, "datafiles/diffPeaks_MES_plastic_vs_fixed_summarized_exp.RDS")

#Conduct MES Plastic Vs MES Parental
diffPeaks_MES_plastic_parental <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "Plasticity_AMT",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "Plastic_MES",
  bgdGroups = "Parental_MES",
  maxCells=4000 # choose high number of max cells to make sure all cells in each condition are included
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

saveRDS(diffPeaks_MES_plastic_parental, "datafiles/diffPeaks_MES_plastic_vs_parental_summarized_exp.RDS")

#Conduct MES Fixed Vs MES Parental
diffPeaks_MES_fixed_parental <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "Plasticity_AMT",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "Fixed_MES",
  bgdGroups = "Parental_MES",
  maxCells=10500 # choose high number of max cells to make sure all cells in each condition are included
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

saveRDS(diffPeaks_MES_fixed_parental, "datafiles/diffPeaks_MES_fixed_vs_parental_summarized_exp.RDS")





#Conduct ADRN Plastic Vs ADRN Fixed
diffPeaks_ADRN_plastic_fixed <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "Plasticity_AMT",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "Plastic_ADRN",
  bgdGroups = "Fixed_ADRN",
  maxCells=10500 # choose high number of max cells to make sure all cells in each condition are included
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

saveRDS(diffPeaks_ADRN_plastic_fixed, "datafiles/diffPeaks_ADRN_plastic_vs_fixed_summarized_exp.RDS")

#Conduct ADRN Plastic Vs ADRN Parental
diffPeaks_ADRN_plastic_parental <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "Plasticity_AMT",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "Plastic_ADRN",
  bgdGroups = "Parental_ADRN",
  maxCells=4000 # choose high number of max cells to make sure all cells in each condition are included
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

saveRDS(diffPeaks_ADRN_plastic_parental, "datafiles/diffPeaks_ADRN_plastic_vs_parental_summarized_exp.RDS")

#Conduct ADRN Fixed Vs ADRN Parental
diffPeaks_ADRN_fixed_parental <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix", 
  groupBy = "Plasticity_AMT",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  normBy="ReadsInTSS",
  useGroups = "Fixed_ADRN",
  bgdGroups = "Parental_ADRN",
  maxCells=10500 # choose high number of max cells to make sure all cells in each condition are included
) # the result is a summarized experiment, useGroups will have >0 FC while bgdGroups will have <0 FC

saveRDS(diffPeaks_ADRN_fixed_parental, "datafiles/diffPeaks_ADRN_fixed_vs_parental_summarized_exp.RDS")





#Convert Differential Sets to Dataframes
#Remove NA values from dataset for FDR
diffPeaks_MES_plastic_fixed_markers <- getMarkers(diffPeaks_MES_plastic_fixed, cutOff = "FDR < 10 ") #Threshold set high to keep all peaks
diffPeaks_MES_plastic_fixed.df <- as.data.frame(diffPeaks_MES_plastic_fixed_markers[[1]]) 
# is sorted by FDR
saveRDS(diffPeaks_MES_plastic_fixed_markers, file="datafiles/diffPeaks_MES_plastic_fixed.RDS")
write.csv(diffPeaks_MES_plastic_fixed.df, "datafiles/diffPeaks_MES_plastic_fixed.csv")

diffPeaks_MES_plastic_parental_markers <- getMarkers(diffPeaks_MES_plastic_parental, cutOff = "FDR < 10 ") #Threshold set high to keep all peaks
diffPeaks_MES_plastic_parental.df <- as.data.frame(diffPeaks_MES_plastic_parental_markers[[1]]) 
# is sorted by FDR
saveRDS(diffPeaks_MES_plastic_parental_markers, file="datafiles/diffPeaks_MES_plastic_parental.RDS")
write.csv(diffPeaks_MES_plastic_parental.df, "datafiles/diffPeaks_MES_plastic_parental.csv")

diffPeaks_MES_fixed_parental_markers <- getMarkers(diffPeaks_MES_fixed_parental, cutOff = "FDR < 10 ") #Threshold set high to keep all peaks
diffPeaks_MES_fixed_parental.df <- as.data.frame(diffPeaks_MES_fixed_parental_markers[[1]]) 
# is sorted by FDR
saveRDS(diffPeaks_MES_fixed_parental_markers, file="datafiles/diffPeaks_MES_fixed_parental.RDS")
write.csv(diffPeaks_MES_fixed_parental.df, "datafiles/diffPeaks_MES_fixed_parental.csv")


diffPeaks_ADRN_plastic_fixed_markers <- getMarkers(diffPeaks_ADRN_plastic_fixed, cutOff = "FDR < 10 ") #Threshold set high to keep all peaks
diffPeaks_ADRN_plastic_fixed.df <- as.data.frame(diffPeaks_ADRN_plastic_fixed_markers[[1]]) 
# is sorted by FDR
saveRDS(diffPeaks_ADRN_plastic_fixed_markers, file="datafiles/diffPeaks_ADRN_plastic_fixed.RDS")
write.csv(diffPeaks_ADRN_plastic_fixed.df, "datafiles/diffPeaks_ADRN_plastic_fixed.csv")

diffPeaks_ADRN_plastic_parental_markers <- getMarkers(diffPeaks_ADRN_plastic_parental, cutOff = "FDR < 10 ") #Threshold set high to keep all peaks
diffPeaks_ADRN_plastic_parental.df <- as.data.frame(diffPeaks_ADRN_plastic_parental_markers[[1]]) 
# is sorted by FDR
saveRDS(diffPeaks_ADRN_plastic_parental_markers, file="datafiles/diffPeaks_ADRN_plastic_parental.RDS")
write.csv(diffPeaks_ADRN_plastic_parental.df, "datafiles/diffPeaks_ADRN_plastic_parental.csv")

diffPeaks_ADRN_fixed_parental_markers <- getMarkers(diffPeaks_ADRN_fixed_parental, cutOff = "FDR < 10 ") #Threshold set high to keep all peaks
diffPeaks_ADRN_fixed_parental.df <- as.data.frame(diffPeaks_ADRN_fixed_parental_markers[[1]]) 
# is sorted by FDR
saveRDS(diffPeaks_ADRN_fixed_parental_markers, file="datafiles/diffPeaks_ADRN_fixed_parental.RDS")
write.csv(diffPeaks_ADRN_fixed_parental.df, "datafiles/diffPeaks_ADRN_fixed_parental.csv")



## Identify Differential Peaks ----
#Each plot is divided the same
#Plastic MES Vs Fixed MES
#Plastic MES Vs Parental MES
#Fixed MES Vs Parental MES
#Plastic ADRN Vs Fixed ADRN
#Plastic ADRN Vs Parental ADRN
#Fixed ADRN Vs Parental ADRN
#threshold cut offs log2foldchange = 1, FDR 0.1 as per ArchR documentation

## Loose thresholds
pdf(paste0("plots/Volcano_plots_", threshold_loose_FDR, "_log2FC_", threshold_loose_logFC, ".pdf"))

### Mean/Log2fc plot
plotMarkers(seMarker = diffPeaks_MES_plastic_fixed, name = "Plastic_MES", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "MA") +
### Volcano plot
plotMarkers(seMarker = diffPeaks_MES_plastic_fixed, name = "Plastic_MES", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
plotMarkers(seMarker = diffPeaks_MES_plastic_parental, name = "Plastic_MES", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "MA") +
### Volcano plot
plotMarkers(seMarker = diffPeaks_MES_plastic_parental, name = "Plastic_MES", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
plotMarkers(seMarker = diffPeaks_MES_fixed_parental, name = "Fixed_MES", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "MA") +
### Volcano plot
plotMarkers(seMarker = diffPeaks_MES_fixed_parental, name = "Fixed_MES", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
plotMarkers(seMarker = diffPeaks_ADRN_plastic_fixed, name = "Plastic_ADRN", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "MA") +
### Volcano plot
plotMarkers(seMarker = diffPeaks_ADRN_plastic_fixed, name = "Plastic_ADRN", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
plotMarkers(seMarker = diffPeaks_ADRN_plastic_parental, name = "Plastic_ADRN", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "MA") +
### Volcano plot
plotMarkers(seMarker = diffPeaks_ADRN_plastic_parental, name = "Plastic_ADRN", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "Volcano")
### Mean/Log2fc plot
plotMarkers(seMarker = diffPeaks_ADRN_fixed_parental, name = "Fixed_ADRN", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "MA") +
### Volcano plot
plotMarkers(seMarker = diffPeaks_ADRN_fixed_parental, name = "Fixed_ADRN", cutOff = paste0("FDR <= ",threshold_loose_FDR, " & abs(Log2FC) >= ", threshold_loose_logFC), plotAs = "Volcano")

dev.off()


#Save ArchR project
saveArchRProject(ArchRProj = project, outputDirectory = "ArchR_projects/project_clones_working/", load = TRUE)


#Plot number of differential genes between MES POT and MES EZH2i OFF
diffPeaks_MES_plastic_fixed.df <- read.csv("datafiles/diffPeaks_MES_plastic_fixed.csv")
diffPeaks_MES_plastic_parental.df <- read.csv("datafiles/diffPeaks_MES_plastic_parental.csv")
diffPeaks_MES_fixed_parental.df <- read.csv("datafiles/diffPeaks_MES_fixed_parental.csv")
diffPeaks_ADRN_plastic_fixed.df <- read.csv("datafiles/diffPeaks_ADRN_plastic_fixed.csv")
diffPeaks_ADRN_plastic_parental.df <- read.csv("datafiles/diffPeaks_ADRN_plastic_parental.csv")
diffPeaks_ADRN_fixed_parental.df <- read.csv("datafiles/diffPeaks_ADRN_fixed_parental.csv")

#Isolate only significant foldchanges using loose cut-off fold change = 1.5, FDR = 0.05
diffPeaks_MES_plastic_fixed.filt.df <- diffPeaks_MES_plastic_fixed.df %>% dplyr::filter((Log2FC >= 1 & FDR < 0.1) | (Log2FC <= -1 & FDR < 0.1))
diffPeaks_MES_plastic_parental.filt.df <- diffPeaks_MES_plastic_parental.df %>% dplyr::filter((Log2FC >= 1 & FDR < 0.1) | (Log2FC <= -1 & FDR < 0.1))
diffPeaks_MES_fixed_parental.filt.df <- diffPeaks_MES_fixed_parental.df %>% dplyr::filter((Log2FC >= 1 & FDR < 0.1) | (Log2FC <= -1 & FDR < 0.1))
diffPeaks_ADRN_plastic_fixed.filt.df <- diffPeaks_ADRN_plastic_fixed.df %>% dplyr::filter((Log2FC >= 1 & FDR < 0.1) | (Log2FC <= -1 & FDR < 0.1))
diffPeaks_ADRN_plastic_parental.filt.df <- diffPeaks_ADRN_plastic_parental.df %>% dplyr::filter((Log2FC >= 1 & FDR < 0.1) | (Log2FC <= -1 & FDR < 0.1))
diffPeaks_ADRN_fixed_parental.filt.df <- diffPeaks_ADRN_fixed_parental.df %>% dplyr::filter((Log2FC >= 1 & FDR < 0.1) | (Log2FC <= -1 & FDR < 0.1))

#Summarise how many differentially accessible regions there are in each comparison
#Add variable to show positive or negative enrichment
diffPeaks_MES_plastic_fixed.filt.df$enrichment <- ifelse(diffPeaks_MES_plastic_fixed.filt.df$Log2FC >= 1, "positive", "negative")
diffPeaks_MES_plastic_parental.filt.df$enrichment <- ifelse(diffPeaks_MES_plastic_parental.filt.df$Log2FC >= 1, "positive", "negative")
diffPeaks_MES_fixed_parental.filt.df$enrichment <- ifelse(diffPeaks_MES_fixed_parental.filt.df$Log2FC >= 1, "positive", "negative")
diffPeaks_ADRN_plastic_fixed.filt.df$enrichment <- ifelse(diffPeaks_ADRN_plastic_fixed.filt.df$Log2FC >= 1, "positive", "negative")
diffPeaks_ADRN_plastic_parental.filt.df$enrichment <- ifelse(diffPeaks_ADRN_plastic_parental.filt.df$Log2FC >= 1, "positive", "negative")
diffPeaks_ADRN_fixed_parental.filt.df$enrichment <- ifelse(diffPeaks_ADRN_fixed_parental.filt.df$Log2FC >= 1, "positive", "negative")

#Summarise number of regions per positive and negative comparison
summary_MES_plastic_fixed <- diffPeaks_MES_plastic_fixed.filt.df %>%
  group_by(enrichment) %>%
  dplyr::summarise(region_number = n()) %>%
  mutate(comparison = "Plastic MES Vs. Fixed MES")

summary_MES_plastic_parental <- diffPeaks_MES_plastic_parental.filt.df %>%
  group_by(enrichment) %>%
  dplyr::summarise(region_number = n())%>%
  mutate(comparison = "Plastic MES Vs. Parental MES")

summary_MES_fixed_parental <- diffPeaks_MES_fixed_parental.filt.df %>%
  group_by(enrichment) %>%
  dplyr::summarise(region_number = n())%>%
  mutate(comparison = "Fixed MES Vs. Parental MES")

summary_ADRN_plastic_fixed <- diffPeaks_ADRN_plastic_fixed.filt.df %>%
  group_by(enrichment) %>%
  dplyr::summarise(region_number = n())%>%
  mutate(comparison = "Plastic ADRN Vs. Fixed ADRN")

summary_ADRN_plastic_parental <- diffPeaks_ADRN_plastic_parental.filt.df %>%
  group_by(enrichment) %>%
  dplyr::summarise(region_number = n())%>%
  mutate(comparison = "Plastic ADRN Vs. Parental ADRN")

summary_ADRN_fixed_parental <- diffPeaks_ADRN_fixed_parental.filt.df %>%
  group_by(enrichment) %>%
  dplyr::summarise(region_number = n())%>%
  mutate(comparison = "Fixed ADRN Vs. Parental ADRN")

#Combine all dataframes to visualise information of differential accessibility
summary.df <- do.call("rbind", list(summary_MES_plastic_fixed, summary_MES_plastic_parental, summary_MES_fixed_parental,
                                    summary_ADRN_plastic_fixed, summary_ADRN_plastic_parental, summary_ADRN_fixed_parental))

#Make negative enrichment value, negative to plot nicer
summary.df$region_number <- ifelse(summary.df$enrichment == "negative", -summary.df$region_number, summary.df$region_number)

#Define colours for positive and negative values
library(RColorBrewer)
display.brewer.pal(n = 8, name = 'RdBu')
brewer.pal(n = 8, name = 'RdBu')
#"#B2182B" "#D6604D" "#F4A582" "#FDDBC7" "#D1E5F0" "#92C5DE" "#4393C3" "#2166AC"

ggplot(summary.df, aes(x=comparison, y=region_number, fill = enrichment)) +
  geom_bar(stat = "identity", position = "stack", colour = "black") +
  theme_bw() +
  scale_fill_manual(values = c("#4393C3", "#D6604D"),
                    name = "Enrichment") +
  labs(y = "Number of Differentiall Accessible Regions", x = "Comparison") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))
ggsave("plots/differential_accessibility_regions_summary.pdf", width = 5, height = 6)

