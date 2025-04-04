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

#Load in sce and metadata
#Load in RDS object with barcodes for processing
nb.seurat <- readRDS("datafiles/nb_seurat_Wbarcodes.RDS")

meta.df <- as.data.frame(nb.seurat@meta.data)

meta.df$replicate[meta.df$replicate == " A"] <- "A"
meta.df$replicate[meta.df$replicate == " B"] <- "B"
meta.df$replicate[meta.df$replicate == " C"] <- "C"

#Add in POT replicates
meta.df$description[meta.df$batch == "multiplex5" & meta.df$Condition == "POT"] <- "POT A"
meta.df$description[meta.df$batch == "multiplex11" & meta.df$Condition == "POT"] <- "POT B"
meta.df$description[meta.df$batch == "multiplex12" & meta.df$Condition == "POT"] <- "POT C"

#Set palette so that 
set.seed(12)
ut_palette = createPalette(382, c("#ff0000", "#00ff00", "#0000ff"))
ut_palette <- as.vector(t(matrix(ut_palette)))
names(ut_palette) = unique(as.character(meta.df$real_bc44))

#Edit names now
unique(meta.df$Condition)
meta.df$Condition <- as.character(meta.df$Condition)
meta.df$Condition[meta.df$Condition == "Cisplatin_1weeksOFF"] <- "Cisplatin_1weekOFF"

meta.df$description[meta.df$description == "Cisplatin_1weeksOFF C"] <- "Cisplatin_1weekOFF C"
meta.df$description[meta.df$description == "Cisplatin_1weeksOFF B"] <- "Cisplatin_1weekOFF B"
meta.df$description[meta.df$description == "Cisplatin_4weekOFF B"] <- "Cisplatin_4weeksOFF B"
meta.df$description[meta.df$description == "Cisplatin_4weekOFF C"] <- "Cisplatin_4weeksOFF C"

meta.df$Condition <- as.factor(meta.df$Condition)
unique(meta.df$Condition)

## Visualise phenotypic plasticity of lineages across conditions ----
#Remove unnecessary variables
subset.meta.df <- meta.df %>% dplyr::select(c(real_bc44, CellID, AMT.state, Condition, description, replicate))

length(unique(subset.meta.df$CellID)) #74891 cells total
length(unique(subset.meta.df$CellID[!is.na(subset.meta.df$real_bc44)])) #29285 cells total  with barcode information

#Remove information without barcodes
subset.meta.df <- subset.meta.df %>% filter(!is.na(real_bc44))
length(unique(subset.meta.df$real_bc44)) #381 lineages

#Remove clones which do not appear more than 3 times
subset.meta.df <- subset.meta.df %>%
  group_by(real_bc44, description) %>%
  mutate(count = n()) %>%
  filter(count >= 3) %>%
  select(-c(count))

length(unique(subset.meta.df$CellID)) #28258 cells total
length(unique(subset.meta.df$real_bc44)) #200 cells total - therefore loosing 52% of clones

# Count unique real_bc44 per description
result <- subset.meta.df %>%
  group_by(description) %>%
  summarise(unique_real_bc44_count = n_distinct(real_bc44))

#Visualise results
print(result)

#Compute how many cell states are present for each clone in each condition
summary.df <- subset.meta.df %>%
  group_by(real_bc44, description, Condition) %>%
  summarize(num_cell_states = n_distinct(AMT.state)) %>%
  ungroup()

prop.df <- summary.df %>%
  group_by(description, Condition) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(num_cell_states, description, Condition) %>%
  mutate(count = n(),
         proportion = count / total) %>%
  ungroup()

#Summarise the number of lineages 
sum <- prop.df %>%
  group_by(description, num_cell_states) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(description) %>%
  mutate(total = sum(count))

write.csv(sum, "datafiles/number_of_barcodes_per_cell_state_group_new.csv")

#Remove duplicated information
plot.prop.df <- prop.df %>%
  dplyr::select(-c(real_bc44, total, count)) %>%
  distinct()

plot.prop.df$num_cell_states <- as.character(plot.prop.df$num_cell_states)

#Set plotting variables
conditions <- c("POT A", "POT B", "POT C", 
                "UT A", "UT B", "UT C", 
                "JQ1_ON A", "JQ1_ON B", "JQ1_ON C", 
                "JQ1_OFF A", "JQ1_OFF B", "JQ1_OFF C",
                "EZH2i_ON A", "EZH2i_ON B", "EZH2i_ON C", 
                "EZH2i_OFF A", "EZH2i_OFF B", "EZH2i_OFF C", 
                "CHIR99021_ON A", "CHIR99021_ON B", "CHIR99021_ON C",
                "CHIR99021_OFF A", "CHIR99021_OFF B", "CHIR99021_OFF C",
                "Cisplatin(1)_ON A", "Cisplatin(1)_ON B", "Cisplatin(1)_ON C",
                "Cisplatin(2)_ON A", "Cisplatin(2)_ON B", "Cisplatin(2)_ON C",
                "Cisplatin_1weekOFF A", "Cisplatin_1weekOFF B", "Cisplatin_1weekOFF C",
                "Cisplatin_4weeksOFF A", "Cisplatin_4weeksOFF B", "Cisplatin_4weeksOFF C")


# Reorder the Category variable based on the custom order
plot.prop.df$description <- factor(plot.prop.df$description, levels = conditions)

plot.prop.df %>%
  ggplot( aes(x = description, y = proportion, fill = num_cell_states)) +
  geom_bar(stat = "summary", position = "stack", colour = "black", aes(x = description, y = proportion, fill = num_cell_states)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("plots/plasticty_lineages_summary_prop_new.pdf", width = 12, height = 5)

#Visualise all experimental arms
p1 <- plot.prop.df %>%
  filter(Condition == "POT" | Condition == "Untreated") %>%
  ggplot( aes(x = description, y = proportion, fill = num_cell_states)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", aes(x = description, y = proportion, fill = num_cell_states)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  ggtitle("Untreated Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p1, "plots/plasticty_lineages_summary_prop_UT_new.pdf", width = 7, height = 5)

p2 <- plot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "JQ1_ON" | Condition == "JQ1_OFF") %>%
  ggplot( aes(x = description, y = proportion, fill = num_cell_states)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", aes(x = description, y = proportion, fill = num_cell_states)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  ggtitle("BRD4i Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p2, "plots/plasticty_lineages_summary_prop_BRD4i_new.pdf", width = 7, height = 5)

p3 <- plot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "EZH2i_ON" | Condition == "EZH2i_OFF") %>%
  ggplot( aes(x = description, y = proportion, fill = num_cell_states)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", aes(x = description, y = proportion, fill = num_cell_states)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  ggtitle("EZH2i Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p3, "plots/plasticty_lineages_summary_prop_EZH2i_new.pdf", width = 7, height = 5)

p4 <- plot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "Cisplatin(1)_ON" | Condition == "Cisplatin_1weekOFF") %>%
  ggplot( aes(x = description, y = proportion, fill = num_cell_states)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", aes(x = description, y = proportion, fill = num_cell_states)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  ggtitle("Cisplatin Short-Recovery Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p4, "plots/plasticty_lineages_summary_prop_Cis1_new.pdf", width = 7, height = 5)

p5 <- plot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "Cisplatin(2)_ON" | Condition == "Cisplatin_4weeksOFF") %>%
  ggplot( aes(x = description, y = proportion, fill = num_cell_states)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", aes(x = description, y = proportion, fill = num_cell_states)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  ggtitle("Cisplatin Long-Recovery Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p5, "plots/plasticty_lineages_summary_prop_Cis2_new.pdf", width = 7, height = 5)

p6 <- plot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "CHIR99021_ON" | Condition == "CHIR99021_OFF") %>%
  ggplot( aes(x = description, y = proportion, fill = num_cell_states)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", aes(x = description, y = proportion, fill = num_cell_states)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  ggtitle("CHIR99021 Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p6, "plots/plasticty_lineages_summary_prop_CHIR99021_new.pdf", width = 7, height = 5)

p7 <- plot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "Cisplatin(2)_ON" | 
           Condition == "Cisplatin_1weekOFF" | Condition == "Cisplatin_4weeksOFF") %>%
  ggplot( aes(x = description, y = proportion, fill = num_cell_states)) +
  geom_bar(stat = "identity", position = "stack", colour = "black", aes(x = description, y = proportion, fill = num_cell_states)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  ggtitle("Cisplatin Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p7, "plots/plasticty_lineages_summary_prop_Cis_AL_newL.pdf", width = 7, height = 5)

p1 + p2 + p3 + p4 + p5 + p7
ggsave("plots/plasticty_lineages_summary_prop_all_new.pdf", width = 12, height = 7)

#Create boxplots to compare SD of each condition
#Calculate mean and SD for each condition
boxplot.prop.df <- plot.prop.df %>%
  dplyr::select(-description) %>%
  group_by(Condition, num_cell_states) %>%
  mutate(mean = mean(proportion),
         SD = sd(proportion)) %>%
  ungroup()

#Set plotting variables
samples <- c("POT", 
             "Untreated",
             "JQ1_ON",
             "JQ1_OFF",
             "EZH2i_ON", 
             "EZH2i_OFF",
             "CHIR99021_ON",
             "CHIR99021_OFF",
             "Cisplatin(1)_ON",
             "Cisplatin(2)_ON",
             "Cisplatin_1weekOFF",
             "Cisplatin_4weeksOFF")


# Reorder the Category variable based on the custom order
boxplot.prop.df$Condition <- factor(boxplot.prop.df$Condition, levels = samples)

p1 <- boxplot.prop.df %>%
  filter(Condition == "POT" | Condition == "Untreated") %>%
  ggplot(aes(x = Condition, y = proportion, fill = num_cell_states)) +
  geom_boxplot(color = "black", aes(ymin = mean - SD, ymax = mean + SD)) +
  geom_point(position = position_dodge(0.75)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  expand_limits(y=c(0,1)) +
  ggtitle("Untreated Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p1, "plots/plasticty_lineages_summary_prop_UT_boxplot_new.pdf", width = 7, height = 5)

p2 <- boxplot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "JQ1_ON" | Condition == "JQ1_OFF") %>%
  ggplot(aes(x = Condition, y = proportion, fill = num_cell_states)) +
  geom_boxplot(color = "black", aes(ymin = mean - SD, ymax = mean + SD)) +
  geom_point(position = position_dodge(0.75)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  expand_limits(y=c(0,1)) +
  ggtitle("BRD4i Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p2, "plots/plasticty_lineages_summary_prop_BRD4i_boxplot_new.pdf", width = 7, height = 5)

p3 <- boxplot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "EZH2i_ON" | Condition == "EZH2i_OFF") %>%
  ggplot(aes(x = Condition, y = proportion, fill = num_cell_states)) +
  geom_boxplot(color = "black", aes(ymin = mean - SD, ymax = mean + SD)) +
  geom_point(position = position_dodge(0.75)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  expand_limits(y=c(0,1)) +
  ggtitle("EZH2i Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p3, "plots/plasticty_lineages_summary_prop_EZH2i_boxplot_new.pdf", width = 7, height = 5)

p4 <- boxplot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "Cisplatin(1)_ON" | Condition == "Cisplatin_1weekOFF") %>%
  ggplot(aes(x = Condition, y = proportion, fill = num_cell_states)) +
  geom_boxplot(color = "black", aes(ymin = mean - SD, ymax = mean + SD)) +
  geom_point(position = position_dodge(0.75)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  expand_limits(y=c(0,1)) +
  ggtitle("Cisplatin Short-Recovery Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p4, "plots/plasticty_lineages_summary_prop_Cis1_boxplot_new.pdf", width = 7, height = 5)

p5 <- boxplot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "Cisplatin(2)_ON" | Condition == "Cisplatin_4weeksOFF") %>%
  ggplot(aes(x = Condition, y = proportion, fill = num_cell_states)) +
  geom_boxplot(color = "black", aes(ymin = mean - SD, ymax = mean + SD)) +
  geom_point(position = position_dodge(0.75)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  expand_limits(y=c(0,1)) +
  ggtitle("Cisplatin Long-Recovery Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p5, "plots/plasticty_lineages_summary_prop_Cis2_boxplot_new.pdf", width = 7, height = 5)

p6 <- boxplot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "CHIR99021_ON" | Condition == "CHIR99021_OFF") %>%
  ggplot(aes(x = Condition, y = proportion, fill = num_cell_states)) +
  geom_boxplot(color = "black", aes(ymin = mean - SD, ymax = mean + SD)) +
  geom_point(position = position_dodge(0.75)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  expand_limits(y=c(0,1)) +
  ggtitle("CHIR99021 Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p6, "plots/plasticty_lineages_summary_prop_CHIR99021_boxplot_new.pdf", width = 7, height = 5)

p7 <- boxplot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "Cisplatin_1weekOFF" | 
           Condition == "Cisplatin(2)_ON" | Condition == "Cisplatin_4weeksOFF") %>%
  ggplot(aes(x = Condition, y = proportion, fill = num_cell_states)) +
  geom_boxplot(color = "black", aes(ymin = mean - SD, ymax = mean + SD)) +
  geom_point(position = position_dodge(0.75)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  expand_limits(y=c(0,1)) +
  ggtitle("Cisplatin Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p7, "plots/plasticty_lineages_summary_prop_Cis_ALL_boxplot_new.pdf", width = 7, height = 5)

p1 + p2 + p3 + p4 + p5 + p7
ggsave("plots/plasticty_lineages_summary_prop_all_boxplot_new.pdf", width = 12, height = 7)

p1 + p2 + p3 + p4 + p5 + p6 + p7
ggsave("plots/plasticty_lineages_summary_prop_all+chir_boxplot_new.pdf", width = 12, height = 10)

#Conduct stats between conditions for visualisation
#Insert 0 values for conditions which do not appear
# Create all possible combinations of 'group' and 'condition'
complete_df <- expand.grid(description = unique(plot.prop.df$description),
                           num_cell_states = unique(plot.prop.df$num_cell_states))

# Merge with the original dataset
df_complete <- merge(complete_df, plot.prop.df, by = c("description", "num_cell_states"), all.x = TRUE)

# Replace NA values in 'value' column with 0
df_complete$proportion[is.na(df_complete$proportion)] <- 0

#Fill in lost values
df_complete$Condition[df_complete$description == "Cisplatin(1)_ON B"] <- "Cisplatin(1)_ON"
df_complete$Condition[df_complete$description == "Cisplatin(1)_ON C"] <- "Cisplatin(1)_ON"
df_complete$Condition[df_complete$description == "Cisplatin(2)_ON A"] <- "Cisplatin(2)_ON"
df_complete$Condition[df_complete$description == "Cisplatin(2)_ON C"] <- "Cisplatin(2)_ON"

#Conduct statistical test between Treatments at each timepoint
n_sample <- unique(df_complete$Condition)
n_group <- unique(df_complete$num_cell_states)

first = T
for (i in 1:length(n_sample)) {
  for(j in 1:length(n_sample)){
    for (k in 1:length(n_group)) {
      
      if(i > j){  
        
        tmp = t.test(df_complete$proportion[df_complete$Condition==n_sample[i] & df_complete$num_cell_states==n_group[k]],
                     df_complete$proportion[df_complete$Condition==n_sample[j] & df_complete$num_cell_states==n_group[k]])
        
        if(first){
          df_pval = data.frame(GroupA = n_sample[i],
                               GroupB = n_sample[j],
                               Cell_Group = n_group[k],
                               pvalue = tmp$p.value)
          first = F
        } else {
          df_pval = rbind.data.frame(df_pval,
                                     data.frame(GroupA = n_sample[i],
                                                GroupB = n_sample[j],
                                                Cell_Group = n_group[k],
                                                pvalue = tmp$p.value))
        }
      }
    }
  }
}

#Conduct stats per experimental arm
UT_df_pval <- df_pval %>%
  filter((GroupA == "POT" | GroupB == "POT") & (GroupA == "Untreated" | GroupB == "Untreated"))

results <- UT_df_pval %>%
  mutate(signif_value = ifelse(UT_df_pval$pvalue < 0.001, "***", 
                               ifelse(UT_df_pval$pvalue <= 0.01, "**",
                                      ifelse(UT_df_pval$pvalue <= 0.05, "*", "NA"))))

write.csv(results, "UT_cellular_cluster_stats_new.csv")

#Conduct stats per experimental arm
Cis_df_pval <- df_pval %>%
  filter(GroupB == "Untreated" & (GroupA == "Cisplatin(1)_ON" | GroupA == "Cisplatin_1weekOFF"))

results <- Cis_df_pval %>%
  mutate(signif_value = ifelse(Cis_df_pval$pvalue < 0.001, "***", 
                               ifelse(Cis_df_pval$pvalue <= 0.01, "**",
                                      ifelse(Cis_df_pval$pvalue <= 0.05, "*", "NA"))))

write.csv(results, "Cis_short_cellular_cluster_stats_new.csv")

#Conduct stats per experimental arm
Cis_df_pval <- df_pval %>%
  filter(GroupB == "Untreated" & (GroupA == "Cisplatin(2)_ON" | GroupA == "Cisplatin_4weeksOFF"))

results <- Cis_df_pval %>%
  mutate(signif_value = ifelse(Cis_df_pval$pvalue < 0.001, "***", 
                               ifelse(Cis_df_pval$pvalue <= 0.01, "**",
                                      ifelse(Cis_df_pval$pvalue <= 0.05, "*", "NA"))))

write.csv(results, "Cis_long_cellular_cluster_stats_new.csv")

#Conduct stats per experimental arm
Cis_df_pval <- df_pval %>%
  filter(GroupB == "Untreated" & (GroupA == "Cisplatin(2)_ON" | GroupA == "Cisplatin_4weeksOFF"
                                  | GroupA == "Cisplatin(1)_ON" | GroupA == "Cisplatin_1weekOFF"))

results <- Cis_df_pval %>%
  mutate(signif_value = ifelse(Cis_df_pval$pvalue < 0.001, "***", 
                               ifelse(Cis_df_pval$pvalue <= 0.01, "**",
                                      ifelse(Cis_df_pval$pvalue <= 0.05, "*", "NA"))))

write.csv(results, "Cis_all_cellular_cluster_stats_new.csv")

#Conduct stats per experimental arm
BRD4i_df_pval <- df_pval %>%
  filter(GroupB == "Untreated" & (GroupA == "JQ1_ON" | GroupA == "JQ1_OFF"))

results <- BRD4i_df_pval %>%
  mutate(signif_value = ifelse(BRD4i_df_pval$pvalue < 0.001, "***", 
                               ifelse(BRD4i_df_pval$pvalue <= 0.01, "**",
                                      ifelse(BRD4i_df_pval$pvalue <= 0.05, "*", "NA"))))

write.csv(results, "BRD4i_cellular_cluster_stats_new.csv")

#Conduct stats per experimental arm
EZH2i_df_pval <- df_pval %>%
  filter(GroupB == "Untreated" & (GroupA == "EZH2i_ON" | GroupA == "EZH2i_OFF"))

results <- EZH2i_df_pval %>%
  mutate(signif_value = ifelse(EZH2i_df_pval$pvalue < 0.001, "***", 
                               ifelse(EZH2i_df_pval$pvalue <= 0.01, "**",
                                      ifelse(EZH2i_df_pval$pvalue <= 0.05, "*", "NA"))))

write.csv(results, "EZH2i_cellular_cluster_stats_new.csv")

#CHIR99021
CHIR_df_pval <- df_pval %>%
  filter(GroupB == "Untreated" & (GroupA == "CHIR99021_ON" | GroupA == "CHIR99021_OFF"))

results <- CHIR_df_pval %>%
  mutate(signif_value = ifelse(CHIR_df_pval$pvalue < 0.001, "***", 
                               ifelse(CHIR_df_pval$pvalue <= 0.01, "**",
                                      ifelse(CHIR_df_pval$pvalue <= 0.05, "*", "NA"))))

write.csv(results, "CHIR_cellular_cluster_stats_new.csv")


#Repeat but with combined Cispaltin replicates
plot.prop.df$Condition <- as.character(plot.prop.df$Condition)
plot.prop.df$description <- as.character(plot.prop.df$description)
plot.prop.df$Condition[plot.prop.df$Condition == "Cisplatin(1)_ON"] <- "Cisplatin_ON"
plot.prop.df$Condition[plot.prop.df$Condition == "Cisplatin(2)_ON"] <- "Cisplatin_ON"

plot.prop.df$description[plot.prop.df$description == "Cisplatin(1)_ON A"] <- "Cisplatin_ON A"
plot.prop.df$description[plot.prop.df$description == "Cisplatin(1)_ON B"] <- "Cisplatin_ON B"
plot.prop.df$description[plot.prop.df$description == "Cisplatin(1)_ON C"] <- "Cisplatin_ON C"
plot.prop.df$description[plot.prop.df$description == "Cisplatin(2)_ON A"] <- "Cisplatin_ON A"
plot.prop.df$description[plot.prop.df$description == "Cisplatin(2)_ON B"] <- "Cisplatin_ON B"
plot.prop.df$description[plot.prop.df$description == "Cisplatin(2)_ON C"] <- "Cisplatin_ON C"

#Calculate mean and SD for each condition
boxplot.prop.df <- plot.prop.df %>%
  dplyr::select(-description) %>%
  group_by(Condition, num_cell_states) %>%
  mutate(mean = mean(proportion),
         SD = sd(proportion)) %>%
  ungroup()

#Set plotting variables
samples <- c("POT", 
             "Untreated",
             "JQ1_ON",
             "JQ1_OFF",
             "EZH2i_ON", 
             "EZH2i_OFF",
             "CHIR99021_ON",
             "CHIR99021_OFF",
             "Cisplatin_ON",
             "Cisplatin_1weekOFF",
             "Cisplatin_4weeksOFF")


# Reorder the Category variable based on the custom order
boxplot.prop.df$Condition <- factor(boxplot.prop.df$Condition, levels = samples)

p7 <- boxplot.prop.df %>%
  filter(Condition == "Untreated" | Condition == "Cisplatin_1weekOFF" | 
           Condition == "Cisplatin_ON" | Condition == "Cisplatin_4weeksOFF") %>%
  ggplot(aes(x = Condition, y = proportion, fill = num_cell_states)) +
  geom_boxplot(color = "black", aes(ymin = mean - SD, ymax = mean + SD)) +
  geom_point(position = position_dodge(0.75)) +
  scale_fill_viridis(discrete = TRUE, option = "D") +
  expand_limits(y=c(0,1)) +
  ggtitle("Cisplatin Arm") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(plot = p7, "plots/plasticty_lineages_summary_prop_v2_Cis_ALL_boxplot_new.pdf", width = 7, height = 5)


#Conduct stats between conditions for visualisation
#Insert 0 values for conditions which do not appear
# Create all possible combinations of 'group' and 'condition'
complete_df <- expand.grid(description = unique(plot.prop.df$description),
                           num_cell_states = unique(plot.prop.df$num_cell_states))

# Merge with the original dataset
df_complete <- merge(complete_df, plot.prop.df, by = c("description", "num_cell_states"), all.x = TRUE)

# Replace NA values in 'value' column with 0
df_complete$proportion[is.na(df_complete$proportion)] <- 0

#Fill in lost values
df_complete$Condition[df_complete$description == "Cisplatin_ON C"] <- "Cisplatin_ON"

#Conduct statistical test between Treatments at each timepoint
n_sample <- unique(df_complete$Condition)
n_group <- unique(df_complete$num_cell_states)

first = T
for (i in 1:length(n_sample)) {
  for(j in 1:length(n_sample)){
    for (k in 1:length(n_group)) {
      
      if(i > j){  
        
        tmp = t.test(df_complete$proportion[df_complete$Condition==n_sample[i] & df_complete$num_cell_states==n_group[k]],
                     df_complete$proportion[df_complete$Condition==n_sample[j] & df_complete$num_cell_states==n_group[k]])
        
        if(first){
          df_pval = data.frame(GroupA = n_sample[i],
                               GroupB = n_sample[j],
                               Cell_Group = n_group[k],
                               pvalue = tmp$p.value)
          first = F
        } else {
          df_pval = rbind.data.frame(df_pval,
                                     data.frame(GroupA = n_sample[i],
                                                GroupB = n_sample[j],
                                                Cell_Group = n_group[k],
                                                pvalue = tmp$p.value))
        }
      }
    }
  }
}

#Conduct stats per experimental arm
Cis_df_pval <- df_pval %>%
  filter((GroupB == "Untreated" | GroupB == "Cisplatin_ON" | GroupB == "Cisplatin_4weeksOFF" | GroupB == "Cisplatin_1weekOFF" ) & 
           (GroupA == "Cisplatin_ON" | GroupA == "Cisplatin_4weeksOFF" | GroupA == "Cisplatin_1weekOFF" | GroupA == "Untreated"))

results <- Cis_df_pval %>%
  mutate(signif_value = ifelse(Cis_df_pval$pvalue < 0.001, "***", 
                               ifelse(Cis_df_pval$pvalue <= 0.01, "**",
                                      ifelse(Cis_df_pval$pvalue <= 0.05, "*", "NA"))))

write.csv(results, "datafiles/SK-N-SH_only/clone_dynamics/Cis_all_v2_cellular_cluster_stats_new.csv")


## ggalluvial plots for each condition - POT Single data ----
#summary.df currently has the barcodes in each condition and how many states they appear in
#Assign replicate values correctly
analysis.df <- meta.df %>% filter(Condition == "POT" & !is.na(real_bc44))
length(unique(analysis.df$real_bc44[analysis.df$Condition == "POT" & (!(is.na(analysis.df$real_bc44)))])) #110 clones

#Summarise how many of these clones are in different cellular clusters
#Calculate the total clonal frequency within each replicate and condition, and therefore its proportion within each sample
test <- analysis.df %>%
  group_by(real_bc44) %>%
  mutate(clones_count = n(),
         total_clones = sum(clones_count),
         percentage = clones_count/total_clones * 100) %>% 
  ungroup() %>%
  dplyr::select(c(real_bc44, CellID, AMT.state, Condition))

#Calculate total frequency of each clone within POT and Untreated conditions
tmp <- test %>%
  group_by(real_bc44, Condition) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = Condition, values_from = AMT.state) %>%
  select(-row) %>%
  ungroup() %>%
  group_by(real_bc44, POT) %>%
  mutate(count = n()) %>%
  ungroup()

#Filter summary.df for condition of interest
barcodes <- summary.df %>% filter(Condition == "POT") %>% dplyr::select(-c(description, Condition))

#Check presence of all barcodes in each dataframe - should be same as clone number above; POT 110 clones
length(intersect(barcodes$real_bc44, tmp$real_bc44)) #81

#Merge two dataframes
intersect(colnames(tmp), colnames(barcodes))
state.df <- left_join(barcodes, tmp, by = c("real_bc44"), relationship = "many-to-many")

#Plot ggalluvial to visualise transitions
#Generate barcode palette
length(unique(state.df$real_bc44)) #81 clones

# Create a factor for real_bc44 based on num_cell_state
state.df <- state.df %>%
  arrange(num_cell_states) %>% # Order the rows by num_cell_states
  mutate(real_bc44 = factor(real_bc44, levels = unique(real_bc44))) 

#Create log value of counts
state.df$log_count <- log2(state.df$count)

# Add labels for num_cell_states to the left-hand side
state.df <- state.df %>%
  mutate(label = as.character(num_cell_states)) # Create labels for num_cell_states

#Remove all other columns
state.df <- state.df %>% dplyr::select(c(CellID, POT, log_count, count, real_bc44, label))

#Save state.df for downstream analysis so we do not have to regenerate it everytime
write.csv(state.df, "datafiles/POT_alluvial_data_new.csv")

state.df %>%
  ggplot( aes(axis1 = real_bc44, axis2 = POT, y = log_count)) +
  scale_x_discrete(limits = c("Cellecta Clone", "POT"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = POT)) +
  geom_stratum(width = 1/8, alpha=0.8, colour="black") +
  geom_text(aes(label = label), stat = "stratum", nudge_x = -0.125) + # Add labels to the left
  scale_fill_manual(values = amt.cols) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))
ggsave("plots/pot_alluvial_1_new.pdf", dpi=700, width=5, height=5)

state.df %>%
  ggplot( aes(axis1=real_bc44, axis2 = POT, y = log_count)) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"POT"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill = POT), width = 1/8, alpha=0.8, colour="black") +
  geom_text(aes(label = label), stat = "stratum", nudge_x = -0.125) + # Add labels to the left
  scale_fill_manual(values = amt.cols) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "none")
ggsave("plots/pot_alluvial_2_new.pdf", dpi=700, width=5, height=5)

state.df %>%
  ggplot( aes(axis1 = real_bc44, axis2 = POT, y = log_count)) +
  scale_x_discrete(limits = c("Cellecta Clone" ,"POT"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=real_bc44), width = 1/8, alpha=0.8, colour="black") +
  geom_text(aes(label = label), stat = "stratum", nudge_x = -0.125) + # Add labels to the left
  scale_fill_manual(values=ut_palette) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "none")
ggsave("plots/pot_alluvial_3_new.pdf", dpi=700, width=5, height=5)

state.df %>%
  ggplot( aes(axis1 = real_bc44, axis2 = label, y = log_count)) +
  scale_x_discrete(limits = c("Cellecta Clone", "Label"), expand = c(.2, .05)) +
  geom_stratum(aes(fill = label), width = 1/8, alpha=0.8, colour="black") +
  geom_text(aes(label = label), stat = "stratum", nudge_x = -0.125) + # Add labels to the left
  scale_fill_viridis(discrete = TRUE, option = "D") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "none")
ggsave("plots/pot_alluvial_label_bar_new.pdf", dpi=700, width=5, height=5)


## ggalluvial plots for each condition - Condition Full data ----

#### N.B. alluvial plots were generated for each individual condition and replicate as follows
# Only the example of Cisplatin(1)_ON A is demonstrated below

#summary.df currently has the barcodes in each condition and how many states they appear in
#Assign replicate values correctly
analysis.df <- meta.df %>% filter(Condition == "Cisplatin(1)_ON" & !is.na(real_bc44))
length(unique(analysis.df$real_bc44[analysis.df$Condition == "Cisplatin(1)_ON" & (!(is.na(analysis.df$real_bc44)))])) #68 clones

#Summarise how many of these clones are in different cellular clusters
#Calculate the total clonal frequency within each replicate and condition, and therefore its proportion within each sample
test <- analysis.df %>%
  group_by(real_bc44) %>%
  mutate(clones_count = n(),
         total_clones = sum(clones_count),
         percentage = clones_count/total_clones * 100) %>% 
  ungroup() %>%
  dplyr::select(c(real_bc44, CellID, AMT.state, Condition, replicate))

#Isolate only those clones which are in a single seurat cluster in POT
test2A <- test %>% filter(replicate == "A") 
test2B <- test %>% filter(replicate == "B") 
test2C <- test %>% filter(replicate == "C") 

tmp_A <- test2A %>%
  group_by(real_bc44, replicate) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = replicate, values_from = AMT.state) %>%
  select(-row) %>%
  group_by(real_bc44, A) %>%
  mutate(count = n())

tmp_B <- test2B %>%
  group_by(real_bc44, replicate) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = replicate, values_from = AMT.state) %>%
  select(-row) %>%
  group_by(real_bc44, B) %>%
  mutate(count = n())

tmp_C <- test2C %>%
  group_by(real_bc44, replicate) %>%
  mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = replicate, values_from = AMT.state) %>%
  select(-row) %>%
  group_by(real_bc44, C) %>%
  mutate(count = n())

#Visualise replicate A
#Filter summary.df for condition of interest
barcodes <- summary.df %>% filter(description == "Cisplatin(1)_ON A") %>% dplyr::select(-c(description, Condition))

#Check presence of all barcodes in each dataframe - should be same as clone number above; 110 clones
length(intersect(barcodes$real_bc44, tmp_A$real_bc44)) #43

#Merge two dataframes
intersect(colnames(tmp_A), colnames(barcodes))
state.df <- left_join(barcodes, tmp_A, by = c("real_bc44"), relationship = "many-to-many")

#Plot ggalluvial to visualise transitions
#Generate barcode palette
length(unique(state.df$real_bc44)) 

# Create a factor for real_bc44 based on num_cell_state
state.df <- state.df %>%
  arrange(num_cell_states) %>% # Order the rows by num_cell_states
  mutate(real_bc44 = factor(real_bc44, levels = unique(real_bc44))) 

#Create log value of counts
state.df$log_count <- log2(state.df$count)

# Add labels for num_cell_states to the left-hand side
state.df <- state.df %>%
  mutate(label = as.character(num_cell_states)) # Create labels for num_cell_states

#Remove all other columns
state.df <- state.df %>% dplyr::select(c(CellID, A, log_count, count, real_bc44, label))

#Save state.df for downstream analysis so we do not have to regenerate it everytime
write.csv(state.df, "datafiles/Cisplatin(1)_ON_A_alluvial_data_new.csv")

state.df %>%
  ggplot( aes(axis1 = real_bc44, axis2 = A, y = log_count)) +
  scale_x_discrete(limits = c("Cellecta Clone", "Cisplatin(1)_ON A"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = A)) +
  geom_stratum(width = 1/8, alpha=0.8, colour="black") +
  geom_text(aes(label = label), stat = "stratum", nudge_x = -0.125) + # Add labels to the left
  scale_fill_manual(values = amt.cols) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))
ggsave("plots/Cisplatin(1)_ON_A_alluvial_1_new.pdf", dpi=700, width=5, height=5)

state.df %>%
  ggplot( aes(axis1=real_bc44, axis2 = A, y = log_count)) +
  scale_x_discrete(limits = c("Cellecta Clone", "Cisplatin(1)_ON A"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill = A), width = 1/8, alpha=0.8, colour="black") +
  geom_text(aes(label = label), stat = "stratum", nudge_x = -0.125) + # Add labels to the left
  scale_fill_manual(values = amt.cols) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "none")
ggsave("plots/Cisplatin(1)_ON_A_alluvial_2_new.pdf", dpi=700, width=5, height=5)

state.df %>%
  ggplot( aes(axis1 = real_bc44, axis2 = A, y = log_count)) +
  scale_x_discrete(limits = c("Cellecta Clone", "Cisplatin(1)_ON A"), expand = c(.2, .05)) +
  geom_alluvium() +
  geom_stratum(aes(fill=real_bc44), width = 1/8, alpha=0.8, colour="black") +
  geom_text(aes(label = label), stat = "stratum", nudge_x = -0.125) + # Add labels to the left
  scale_fill_manual(values=ut_palette) +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "none")
ggsave("plots/Cisplatin(1)_ON_A_alluvial_3_new.pdf", dpi=700, width=5, height=5)

state.df %>%
  ggplot( aes(axis1 = real_bc44, axis2 = label, y = log_count)) +
  scale_x_discrete(limits = c("Cellecta Clone", "Label"), expand = c(.2, .05)) +
  geom_stratum(aes(fill = label), width = 1/8, alpha=0.8, colour="black") +
  geom_text(aes(label = label), stat = "stratum", nudge_x = -0.125) + # Add labels to the left
  scale_fill_viridis(discrete = TRUE, option = "D") +
  theme_cowplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size=8),
        legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "none")
ggsave("Cisplatin(1)_ON_A_alluvial_label_bar_new.pdf", dpi=700, width=5, height=5)




# Summary and stats of alluvial ----
#Load in all dataframes of information plotted on ggalluvial plots
alluvial.list <- list.files("datafiles", pattern="_alluvial_data_new.csv", full.names = TRUE)

# Initialize a list to store dataframes
df_list <- list()

# Loop through the files and process each
for (x in seq_along(alluvial.list)) {
  # Extract filename without path and remove ".csv"
  file_name <- gsub("_alluvial_data_new.csv$", "", basename(alluvial.list[x])) 
  
  # Ensure the name is a valid R variable name
  file_name <- make.names(file_name)
  
  # Read the CSV file
  df <- read.csv(alluvial.list[x], sep = ",", header = TRUE, stringsAsFactors = FALSE)
  
  # Add a new column with the dataframe name
  df$name <- file_name
  
  # Store in list
  df_list[[file_name]] <- df
}

# Merge all dataframes into one
merged_df <- bind_rows(df_list, .id = "source")
merged_df$source <- NULL
merged_df$X <- NULL


#Extract info per replicate or POT so we can remove excess information
merged_A.df <- merged_df %>%
  dplyr::select(c(A, count, real_bc44, label, name)) %>%
  filter(!is.na(A)) %>%
  distinct()

merged_B.df <- merged_df %>%
  dplyr::select(c(B, count, real_bc44, label, name)) %>%
  filter(!is.na(B)) %>%
  distinct()

merged_C.df <- merged_df %>%
  dplyr::select(c(C, count, real_bc44, label, name)) %>%
  filter(!is.na(C)) %>%
  distinct()

#Redefine label for 1 and more than 1
merged_A.df$label[merged_A.df$label == 1] <- "=1"
merged_A.df$label[merged_A.df$label == 2] <- ">1"
merged_A.df$label[merged_A.df$label == 3] <- ">1"

merged_B.df$label[merged_B.df$label == 1] <- "=1"
merged_B.df$label[merged_B.df$label == 2] <- ">1"
merged_B.df$label[merged_B.df$label == 3] <- ">1"

merged_C.df$label[merged_C.df$label == 1] <- "=1"
merged_C.df$label[merged_C.df$label == 2] <- ">1"
merged_C.df$label[merged_C.df$label == 3] <- ">1"

#Extract only condition, reomve replicates
merged_A.df$Condition <- sub("^(.*)_(.*?)$", "\\1", merged_A.df$name)
merged_B.df$Condition <- sub("^(.*)_(.*?)$", "\\1", merged_B.df$name)
merged_C.df$Condition <- sub("^(.*)_(.*?)$", "\\1", merged_C.df$name)

#For each replicate/condition summarise the amount of cellID which are in each num_cell_states/label condition
#We do not need the barcode information now as we just need to summarise the number of cells which are deemed "plastic" at that moment
sum_A.df <- merged_A.df %>%
  group_by(Condition, label) %>%
  summarise(freq = sum(count)) %>%
  ungroup()

sum_B.df <- merged_B.df %>%
  group_by(Condition, label) %>%
  summarise(freq = sum(count)) %>%
  ungroup()

sum_C.df <- merged_C.df %>%
  group_by(Condition, label) %>%
  summarise(freq = sum(count)) %>%
  ungroup()

#Merge all
df_list <- list(sum_A.df, sum_B.df, sum_C.df)
sum.df <- do.call("rbind", df_list)
rm(sum_A.df, sum_B.df, sum_C.df, merged_A.df, merged_B.df, merged_C.df)

#Aggregate across replicates
sum.df <- sum.df %>%
  group_by(Condition, label) %>%
  summarise(freq_all = sum(freq)) %>% ungroup()

#Conduct statistical test between conditions - chi squared test
#Create contingency table with rows containing treatments and columns containing groups
library(tidyr)
library(tibble)
library(dplyr)

# Convert tibble to matrix
sum_matrix <- sum.df %>%
  pivot_wider(names_from = label, values_from = freq_all) %>%
  column_to_rownames(var = "Condition") %>%
  as.matrix()

# Conduct the chi-squared test - UT vs each other condition
# Null hypothesis that plasticity groups and treatment are independent
# p < 0.05 there for condition and plasticity groups are dependent
# Define "Untreated" as BRD4i_OFF
n_samples <- unique(rownames(sum_matrix))

first = T
for (i in 1:length(n_samples)) {
  for (j in 1:length(n_samples)) {
    
    if (i > j) {
      
      #Isolate only conditions being used
      subset_data <- sum_matrix[c(i, j), ] 
      
      # Perform Chi-Square test
      chi_test <- chisq.test(subset_data)
      
      if(first){
        df_pval = data.frame(GroupB = n_samples[i],
                             GroupA = n_samples[j],
                             Chi_Square_Stat = chi_test$statistic,
                             pvalue = chi_test$p.value)
        first = F
      } else {
        df_pval = rbind.data.frame(df_pval,
                                   data.frame(GroupB = n_samples[i],
                                              GroupA = n_samples[j],
                                              Chi_Square_Stat = chi_test$statistic,
                                              pvalue = chi_test$p.value))
      }
    }
  }
}

# Apply Bonferroni correction for multiple comparisons
df_pval$p_adjusted <- p.adjust(df_pval$pvalue, method = "bonferroni")

#Visualise pvalues of each comparison
#Extract pvalues from each comparison for plotting
corr.df <- df_pval %>%
  select(c(GroupA, GroupB, p_adjusted))

#Remove rownames
rownames(corr.df) <- NULL

#Calculate ratio between group =1 and group >1
sum_matrix.df <- as.data.frame(sum_matrix)
sum_matrix.df$Condition <- rownames(sum_matrix.df)
rownames(sum_matrix.df) <- NULL
sum_matrix.df$total <- sum_matrix.df$`=1` + sum_matrix.df$`>1`
sum_matrix.df$ratio <- sum_matrix.df$`>1` / sum_matrix.df$total

#Calculate change up or down from UT
sum_matrix.df <- sum_matrix.df %>%
  mutate(change = ratio - ratio[Condition == "Untreated"])

#Remove unneeded columns and data
tmp <- sum_matrix.df %>%
  select(c(Condition, change))

untreated.df <- corr.df %>%
  filter(GroupB == "Untreated")

#Merge tmp and untreated.df based on Condition (tmp) and GroupA (untreated.df)
merged_df <- left_join(tmp, untreated.df, by = c("Condition" = "GroupA"))

#Remove untreated NA column and filter only significant p-values (p_adjusted < 0.05)
filtered_df <- merged_df %>%
  filter(!is.na(GroupB))

#Define size based on the sign of change
filtered_df <- filtered_df %>%
  mutate(size_var = ifelse(change >= 0, "positive", "negative")) 

filtered_df <- filtered_df %>%
  mutate(size_var = ifelse(p_adjusted >= 0.05, "not-significant", size_var))

#Visualise data 
filtered_df %>%
  mutate(Condition = fct_relevel(Condition, "POT", "Cisplatin1ON", "Cisplatin2ON",
                                 "Cisplatin1weeksOFF", "Cisplatin4weeksOFF",
                                 "CHIR99021_ON", "CHIR99021_OFF",
                                 "BRD4i_ON", "BRD4i_OFF",
                                 "EZH2i_ON", "EZH2i_OFF")) %>%
  ggplot( aes(x = Condition, y = GroupB, fill = size_var, size = -log10(p_adjusted))) +
  geom_point(shape = 21, color = "black") +   
  scale_fill_manual(values = c("red", "grey", "blue")) +
  scale_size(range = c(3, 12)) +
  coord_fixed() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  labs(title = "Significant Condition Comparisons (p < 0.05)",
       x = "Condition", y = "",
       size = "- log 2 Adjusted P-value")
ggsave("chi-squared_results_plasticity_all.pdf", width = 12, height = 5)


#### N.B. plots were generated per condition
## Only Cisplatin 1 and associated samples are described here

filtered_df %>%
  filter(Condition == "POT" | Condition == "Cisplatin1ON" | Condition == "Cisplatin1weeksOFF") %>%
  mutate(Condition = fct_relevel(Condition, "POT", "Cisplatin1ON", "Cisplatin1weeksOFF")) %>%
  ggplot( aes(x = Condition, y = GroupB, fill = size_var, size = -log10(p_adjusted))) +
  geom_point(shape = 21, color = "black") +   
  scale_fill_manual(values = c("grey", "blue")) +
  scale_size(range = c(3, 12)) +
  coord_fixed() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  labs(title = "Significant Condition Comparisons (p < 0.05)",
       x = "Condition", y = "",
       size = "- log 2 Adjusted P-value")
ggsave("chi-squared_results_plasticity_Cis_short.pdf", width = 7, height = 5)


#Save corr.df and filtered.df for later use
write.csv(corr.df, "correlation_chi-squared_data_plasticity.csv", row.names = FALSE)
write.csv(filtered_df, "UTvsall_chi-squared_data_plasticity.csv", row.names = FALSE)
write.csv(df_pval, "chi-squared_raw_plasticity.csv", row.names = FALSE)





## Correlation between clone frequency and phenotypic plasticity ----
## Load packages and dependencies
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

#Define colour palettes
amt.cols <- c("ADRN" = "#990099",
              "intermediate" = "lightgrey",
              "MES" = "#F37735")


#Load in RDS object with barcodes for processing
nb.seurat <- readRDS("nb_seurat_Wbarcodes.RDS")

meta.df <- as.data.frame(nb.seurat@meta.data)

#Edit errors in naming
meta.df$Condition <- as.character(meta.df$Condition)
meta.df$description <- as.character(meta.df$description)
meta.df$Condition[meta.df$Condition == "Cisplatin_1weeksOFF"] <- "Cisplatin_1weekOFF"

meta.df$description[meta.df$description == "Cisplatin_1weeksOFF C"] <- "Cisplatin_1weekOFF C"
meta.df$description[meta.df$description == "Cisplatin_1weeksOFF B"] <- "Cisplatin_1weekOFF B"
meta.df$description[meta.df$description == "Cisplatin_4weekOFF B"] <- "Cisplatin_4weeksOFF B"
meta.df$description[meta.df$description == "Cisplatin_4weekOFF C"] <- "Cisplatin_4weeksOFF C"

meta.df$replicate[meta.df$replicate == " A"] <- "A"
meta.df$replicate[meta.df$replicate == " B"] <- "B"
meta.df$replicate[meta.df$replicate == " C"] <- "C"

#Add in POT replicates
meta.df$description[meta.df$batch == "multiplex5" & meta.df$Condition == "POT"] <- "POT A"
meta.df$description[meta.df$batch == "multiplex11" & meta.df$Condition == "POT"] <- "POT B"
meta.df$description[meta.df$batch == "multiplex12" & meta.df$Condition == "POT"] <- "POT C"

#Isolate CellID, real_bc44 and num_cell_states
#Remove unnecessary variables
subset.meta.df <- meta.df %>% dplyr::select(c(real_bc44, CellID, AMT.state, Condition, description, replicate))

#Remove information without barcodes
subset.meta.df <- subset.meta.df %>% filter(!is.na(real_bc44))

#Subset for barcodes 
subset.meta.df <- subset.meta.df %>%
  group_by(real_bc44, description) %>%
  mutate(count = n()) %>%
  filter(count >= 3) %>%
  select(-c(count))

# Count unique real_bc44 per description
result <- subset.meta.df %>%
  group_by(description) %>%
  summarise(unique_real_bc44_count = n_distinct(real_bc44))

#Compute how many cell states are present for each clone in each condition
summary.df <- subset.meta.df %>%
  group_by(real_bc44, description, Condition) %>%
  summarize(num_cell_states = n_distinct(AMT.state)) %>%
  ungroup() 

#Combine num_cell_state information with remainder of meta data
sub.meta.df <- meta.df %>% dplyr::select(c(CellID, real_bc44, Condition, description, replicate))
sub.meta.df <- sub.meta.df %>% filter(!is.na(real_bc44)) #29285 cells
rownames(sub.meta.df) <- NULL

#Summarise the count of each clone in each description
tmp <- sub.meta.df %>%
  group_by(real_bc44, description, Condition, replicate) %>%
  summarise(total_clone_count = n())

intersect(colnames(tmp), colnames(summary.df))
state.df <- merge(summary.df, tmp, by = c("real_bc44", "Condition", "description"), relationship = "many-to-many")

#Conduct Spearmans correlation
#Isolate dataframe for each condition
# Function to compute correlation for each group
# Define function with error handling
compute_correlation <- function(subset_data) {
  list(Spearman = cor(subset_data$total_clone_count, subset_data$num_cell_states, method = "spearman"))
}

cor_results <- lapply(split(state.df, state.df$description), compute_correlation)
print(cor_results)

#Log2 clone frequency for easier visuals
state.df <- state.df %>%
  mutate(log_count = log(total_clone_count))

# Correlation plot of raw data, faceted by condition
library("ggpubr")
ggscatter(state.df, x = "num_cell_states", y = "log_count", 
          color = "replicate",
          add = "reg.line", conf.int = TRUE, conf.int.level = 0.95,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Cell State Occupancy", ylab = "Clone Frequency") +
  facet_wrap(~Condition)
ggsave("plots/scClonal_dynamics/freq_occupancy_correlation_new.pdf", width = 12, height = 12)

## Summary for thesis plots ----
#Remove cells with no cellecta barcode
barcode.df <- subset.meta.df %>%
  filter(!is.na(real_bc44)) %>%
  select(c(real_bc44, Condition, replicate, description))

rownames(barcode.df) <- NULL

#Compute barcode information
sum.barcode.df <- barcode.df %>%
  group_by(description, Condition, real_bc44) %>%
  mutate(N = n()) %>% ungroup() %>%
  group_by(description, Condition) %>%
  mutate(N_tot = n(),
         Proportion = N/N_tot,
         Percentage = Proportion * 100) %>% ungroup() %>%
  distinct()

tmp <- sum.barcode.df %>%
  group_by(description, Condition) %>%
  summarise(count = n_distinct(real_bc44)) %>% ungroup()

#Visualise barcode information
tmp %>%
  mutate(Condition = fct_relevel(Condition,
                                 "POT", "Untreated",
                                 "Cisplatin(1)_ON", "Cisplatin(2)_ON", "Cisplatin_1weekOFF", "Cisplatin_4weeksOFF",
                                 "CHIR99021_ON", "CHIR99021_OFF",
                                 "JQ1_ON", "JQ1_OFF",
                                 "EZH2i_ON", "EZH2i_OFF")) %>%
  ggplot( aes(x = Condition, y = count)) +
  geom_boxplot(colour = "black") +
  geom_point()+
  expand_limits(y = c(0, 80)) +
  ggtitle("Clone Frequency Across Experimental Samples") +
  xlab("Experimental Sample") +
  ylab("Number of Unique Cellecta Clones")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 1, hjust = 1),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))
ggsave("clones_across_scRNA_boxplot_new.pdf", width = 7, height = 5)  

tmp %>%
  filter(Condition == "POT" | Condition == "Untreated") %>%
  mutate(Condition = fct_relevel(Condition,
                                 "POT", "Untreated")) %>%
  ggplot( aes(x = Condition, y = count)) +
  geom_boxplot(colour = "black") +
  geom_point()+
  expand_limits(y = c(0, 80)) +
  ggtitle("Clone Frequency Across Experimental Samples") +
  xlab("Experimental Sample") +
  ylab("Number of Unique Cellecta Clones")+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 1, hjust = 1),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))
ggsave("POT_clones_across_scRNA_boxplot_new.pdf", width = 3, height = 5)  

#Conduct stats between populations
stats.df <- tmp
n_sample <- unique(stats.df$Condition)

first = T
for (i in 1:length(n_sample)) {
  for (j in 1:length(n_sample)) {
    if(i > j) {
      
      tmp = t.test(stats.df$count[stats.df$Condition == n_sample[i]],
                   stats.df$count[stats.df$Condition == n_sample[j]])
      
      if(first) {
        df_pval = data.frame(GroupA = n_sample[i],
                             GroupB = n_sample [j],
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

df_pval <- df_pval %>% filter((GroupA == "POT" | GroupB == "POT"))

results <- df_pval %>%
  mutate(signif_value = ifelse(df_pval$pvalue < 0.001, "***",
                               ifelse(df_pval$pvalue <= 0.01, "**",
                                      ifelse(df_pval$pvalue <= 0.05, "*", "NA"))))
write.csv(results, "stats_barcodes_between_samples_new.csv")


#Visualise clone barcodes as bar chart for each sample
set.seed(12)
length(unique(sum.barcode.df$real_bc44)) #382
P382 = createPalette(382, c("#ff0000", "#00ff00", "#0000ff"), M = 382)
P382 <- as.vector(t(matrix(P382)))
names(P382) <- unique(as.character(sum.barcode.df$real_bc44), row.names = FALSE)

#Plot each bar set individually
sum.barcode.df %>%
  ggplot( aes(x = replicate, y = Proportion, fill = reorder(real_bc44, -Proportion))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Condition, ncol = 6) +
  scale_fill_manual(values = P382) +
  ylab("Cellecta Barcode (%)") +
  xlab(" ") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        strip.background = element_rect(color = "white", fill = "white", linewidth = 0.25,
                                        linetype = "solid"))
ggsave("all_bar_new.pdf", width = 12, height = 5)

### N.B. plots were generated for each condition individually
# Only POT and untreated samples are described below 

sum.barcode.df %>%
  filter(Condition == "POT") %>%
  ggplot( aes(x = replicate, y = Proportion, fill = reorder(real_bc44, -Proportion))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Condition, scales = "free", ncol = 6) +
  expand_limits(y = c(0,1)) +
  scale_fill_manual(values = P382) +
  ylab("Cellecta Barcode (%)") +
  xlab(" ") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        strip.background = element_rect(color = "white", fill = "white", linewidth = 0.25,
                                        linetype = "solid"))
ggsave("POT_bar.pdf", width = 3, height = 5)

sum.barcode.df %>%
  filter(Condition == "Untreated") %>%
  ggplot( aes(x = replicate, y = Proportion, fill = reorder(real_bc44, -Proportion))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Condition, scales = "free", ncol = 6) +
  expand_limits(y = c(0,1)) +
  scale_fill_manual(values = P382) +
  ylab("Cellecta Barcode (%)") +
  xlab(" ") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.position = "none",
        strip.background = element_rect(color = "white", fill = "white", linewidth = 0.25,
                                        linetype = "solid"))
ggsave("UT_bar.pdf", width = 3, height = 5)

