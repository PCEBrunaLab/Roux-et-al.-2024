#ClonalDynamics.R
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

#Load in barcode summary file
all_barcodes <- read.csv("datafiles/all_barcodes.tsv")
barcodes = unique(all_barcodes$real_bc44) #8515 unique barcodes

## Redefine description variables 
all_barcodes$recovery <- ifelse(grepl("_R", all_barcodes$sample_id), TRUE, FALSE)
all_barcodes$sample_id <- str_replace(all_barcodes$sample_id,"_R","")

all_barcodes[c('pool', 'sample_number', 'condition', 'replicate', 'triplicate', 'sample_number_2')] <- str_split_fixed(all_barcodes$sample_id, '_', 6)
all_barcodes <- all_barcodes %>% select(-c(sample_number_2, pool))

all_barcodes <- all_barcodes %>% mutate(triplicate = ifelse(condition == "POT", replicate, triplicate))
all_barcodes <- all_barcodes %>% mutate(replicate = ifelse(condition == "POT", NA, replicate))

#Some files were incorrectly labelled before sending to TPU, edit names of these files
all_barcodes$condition[all_barcodes$sample_number == "0049" | all_barcodes$sample_number == "0050" | all_barcodes$sample_number == "0051"] <- "JQ1"
all_barcodes$condition[all_barcodes$sample_number == "0052" | all_barcodes$sample_number == "0053" | all_barcodes$sample_number == "0054"] <- "Cis"

#Save datafile
write.csv(all_barcodes, "datafiles/full_barcodes.csv", row.names = FALSE)

#Calculate percentage representation of barcodes within samples
all_barcodes <- all_barcodes %>%
  group_by(condition, recovery, real_bc44, replicate) %>%
  mutate(count = sum(N)) %>% ungroup() %>%
  group_by(condition, replicate, recovery) %>%
  mutate(total_count = sum(count)) %>% ungroup() %>%
  mutate(barcode_percentage = (count/total_count)*100)

#Alter sample descriptions for plotting ease
all_barcodes <- all_barcodes %>% mutate(full_sample = ifelse(recovery == TRUE, paste0(condition, " ", "R", " ", replicate), 
                              paste0(condition, " ", "", replicate)))
all_barcodes <- all_barcodes %>% mutate(sample = ifelse(recovery == TRUE, paste0(condition, " ", "R"), paste0(condition)))

#Summarise the number of unique barcodes within each sample
barcode_summary <- all_barcodes %>% select(c(real_bc44, full_sample, sample)) %>%
  group_by(full_sample, sample) %>% summarise(count = n_distinct(real_bc44))
barcode_summary$sample <- factor(barcode_summary$sample, 
                                 levels=c("POT", "UT", "UT R", "Cis", "Cis R", "JQ1", "JQ1 R"))

#Plot Extended Figure 3b
options(scipen = 999)
ggplot(data=barcode_summary, aes(x=sample, y=count)) +
  geom_boxplot() +
  ggtitle("Clone Frequency Across Experimental Samples") + 
  xlab("Experimental Sample")+
  ylab("Number of Unique Cellecta Clones")+
  expand_limits(y=c(0,3000)) +
  scale_y_continuous(breaks = seq(0, 3000, 500)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))

#Remove technical triplicate column and remove duplicated rows
all_barcodes <- all_barcodes %>% mutate(replicate = ifelse(condition == "POT", paste0(triplicate), replicate))
test <- all_barcodes %>% select(-c(triplicate, N, sample_id, sample_number))
barcodes.df <- test %>% dplyr::distinct()

#Define clonal reductions within samples compared to POT
#Define proportion instead of percentage as per proper shannon calculations
all_barcodes <- all_barcodes %>% mutate(proportion = barcode_percentage/100)

shannon_summary <- all_barcodes %>% group_by(full_sample) %>%
  mutate(shannon = (-1)* (sum(proportion * log(proportion))),
         richness = n_distinct(real_bc44)) %>% ungroup() %>%
  select(c(shannon, sample, richness)) %>%
  distinct()

#Visualise shannon diversity
shannon_summary$sample <- factor(shannon_summary$sample, 
                                 levels=c("POT", "UT", "UT R", "Cis", "Cis R", "JQ1", "JQ1 R"))

#Plot Figure 2b
ggplot(shannon_summary, aes(x=sample, y=shannon)) +
  geom_bar(stat="summary", fill = "white", colour = "black") +
  geom_point()+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))

#Conduct pair-wise t-test between sample and UT to identify significant reductions in population
## Remove POT as this does not have biological replicates
shannon_summary <- shannon_summary %>% filter(sample != "POT")
n_sample <- unique(shannon_summary$sample)
first = T
for (i in 1:length(n_sample)) {
  for(j in 1:length(n_sample)){
    if(i > j){  
      
      tmp = t.test(shannon_summary$shannon[shannon_summary$sample==n_sample[i]],
                   shannon_summary$shannon[shannon_summary$sample==n_sample[j]])
      
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
df_pval$fdpvaluer<= 0.05
df_pval <- df_pval %>%
  mutate(signif = ifelse(df_pval$pvalue <= 0.05, "TRUE", "FALSE"),
         signif_value = ifelse(df_pval$pvalue <= 0.001, "***", 
                               ifelse(df_pval$pvalue <= 0.01, "**",
                                      ifelse(df_pval$pvalue <= 0.05, "*", "NA"))))

#Plot Figure 2c
all_barcodes %>%
  group_by(full_sample,triplicate,real_bc44, sample)%>%
  summarise(abundance = sum(N)) %>%
  ungroup() %>% group_by(full_sample, sample) %>%
  arrange(desc(abundance)) %>%
  mutate(cumulative_prop = cumsum(abundance)/sum(abundance),
         id = row_number()) %>%
  ungroup() %>% group_by(sample, id) %>%
  summarise(min = min(cumulative_prop),
            max = max(cumulative_prop),
            avg = mean(cumulative_prop)) %>%
  ungroup() %>%
  ggplot(aes(x=id,col=sample))+
  geom_ribbon(aes(ymin=min, ymax=max, fill=sample),alpha=0.5,lty=0) +
  geom_line(aes(y=avg))+
  scale_colour_manual(values = c("#4D8AC6", "#9B62A7", "#FFCCCC", "#DF4828", "#DDAA3C", "#A6BE54", "#60AB9E"),
                      breaks = c("POT", "UT", "UT R", "Cis", "Cis R", "JQ1", "JQ1 R")) +
  scale_fill_manual(values = c("#4D8AC6", "#9B62A7", "#FFCCCC", "#DF4828", "#DDAA3C", "#A6BE54", "#60AB9E"),
                    breaks = c("POT", "UT", "UT R", "Cis", "Cis R", "JQ1", "JQ1 R")) +
  scale_x_log10()+
  theme_bw() + 
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))

#Plot Venn diagrams to look at overlap of cellecta barcodes between biological replciates
## Technical replicates in case of POT samples

#Define the barcodes present in each sample condition
POT.1.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$replicate == "1" & barcodes.df$sample == "POT"])
POT.2.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$replicate == "2" & barcodes.df$sample == "POT"])
POT.3.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$replicate == "3" & barcodes.df$sample == "POT"])
UT.A.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "UT A"])
UT.B.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "UT B"])
UT.C.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "UT C"])
UT.R.A.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "UT R A"])
UT.R.B.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "UT R B"])
UT.R.C.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "UT R C"])
Cis.A.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "Cis A"])
Cis.B.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "Cis B"])
Cis.C.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "Cis C"])
Cis.R.A.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "Cis R A"])
Cis.R.B.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "Cis R B"])
Cis.R.C.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "Cis R C"])
JQ1.A.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "JQ1 A"])
JQ1.B.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "JQ1 B"])
JQ1.C.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "JQ1 C"])
JQ1.R.A.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "JQ1 R A"])
JQ1.R.B.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "JQ1 R B"])
JQ1.R.C.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "JQ1 R C"])

#Plot Extended Figure 3c
pot.samples <- list(`1`=POT.1.barcodes ,`2`=POT.2.barcodes, `3`=POT.3.barcodes)
ggVennDiagram(pot.samples, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white", high = "blue") +
  scale_color_manual(values = c("black", "black", "black")) +
  ggtitle("Common Barcodes Between POT Triplicates")

ut.samples <- list(A=UT.A.barcodes ,B=UT.B.barcodes, C=UT.C.barcodes)
ggVennDiagram(ut.samples, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white", high = "blue") +
  scale_color_manual(values = c("black", "black", "black")) +
  ggtitle("Common Barcodes Between UT Replciates")

ut.r.samples <- list(A=UT.R.A.barcodes ,B=UT.R.B.barcodes, C=UT.R.C.barcodes)
ggVennDiagram(ut.r.samples, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white", high = "blue") +
  scale_color_manual(values = c("black", "black", "black")) +
  ggtitle("Common Barcodes Between UT OFF Replciates")

cis.samples <- list(A=Cis.A.barcodes ,B=Cis.B.barcodes, C=Cis.C.barcodes)
ggVennDiagram(cis.samples, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white", high = "blue") +
  scale_color_manual(values = c("black", "black", "black")) +
  ggtitle("Common Barcodes Between Cisplatin Replciates")

cis.r.samples <- list(A=Cis.R.A.barcodes ,B=Cis.R.B.barcodes, C=Cis.R.C.barcodes)
ggVennDiagram(cis.r.samples, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white", high = "blue") +
  scale_color_manual(values = c("black", "black", "black")) +
  ggtitle("Common Barcodes Between Cisplatin OFF Replciates")

jq1.samples <- list(A=JQ1.A.barcodes ,B=JQ1.B.barcodes, C=JQ1.C.barcodes)
ggVennDiagram(jq1.samples, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white", high = "blue") +
  scale_color_manual(values = c("black", "black", "black")) +
  ggtitle("Common Barcodes Between JQ1 Replciates")

jq1.r.samples <- list(A=JQ1.R.A.barcodes ,B=JQ1.R.B.barcodes, C=JQ1.R.C.barcodes)
ggVennDiagram(jq1.r.samples, label_alpha = 0) +
  ggplot2::scale_fill_gradient(low="white", high = "blue") +
  scale_color_manual(values = c("black", "black", "black")) +
  ggtitle("Common Barcodes Between JQ1 OFF Replciates")

#Isolate barcodes which are present in POT and process only these going forward
POT.barcodes <- unique(barcodes.df$real_bc44[barcodes.df$full_sample == "POT NA"]) #1437 clones
filtered.barcode.df <- barcodes.df %>%filter(real_bc44 %in% POT.barcodes)

#Define colour palette for plotting 1437 unique cellecta clones
set.seed(12)
P1437 = createPalette(1437, c("#ff0000", "#00ff00", "#0000ff"), M=1437)
P1437 <- as.vector(t(matrix(P1437)))
names(P1437) = unique(filtered.barcode.df$real_bc44)

filtered.barcode.df <- filtered.barcode.df %>% mutate(sample = ifelse(recovery == TRUE, paste0(condition, " ", "R"), paste0(condition))) 
filtered.barcode.df <- filtered.barcode.df %>% mutate(replicate = ifelse(condition == "POT", NA, replicate))

#Isolate samples which are not the POT
no.pot.filtered.df <- filtered.barcode.df %>% filter(condition != "POT")

#Calculate percentage representation of barcodes within POT only
pot.filtered.barcode.df <- filtered.barcode.df %>%
  filter(condition == "POT") %>% group_by(real_bc44) %>%
  mutate(sum_perc = sum(barcode_percentage)) %>% ungroup() %>%
  group_by(real_bc44) %>%
  mutate(barcode_percentage = sum_perc) %>% select(-c(sum_perc))

#Merge this informaiton back with original samples
test <- rbind(no.pot.filtered.df, pot.filtered.barcode.df)

#Plot Extended Figure 3a
test %>%
  mutate(sample = factor(sample, levels=c("POT", "UT", "UT R", "Cis", "Cis R", "JQ1", "JQ1 R"))) %>%
  ggplot( aes(x=replicate, y=barcode_percentage, fill=reorder(real_bc44, barcode_percentage))) +
  geom_bar(stat="summary", position="stack") +
  facet_wrap(~sample, scales = "free", ncol = 7)+
  scale_fill_manual(values = P1437) +
  ylab("Cellecta Barcode (%)")+
  xlab(" ")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=8),
        legend.position = "none",
        aspect.ratio = 1,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Calculate percentage reduction in each clonal population for each sample compared to POT
reduction.df <- filtered.barcode.df %>%
  group_by(full_sample, sample) %>%
  summarise(barcodes = n_distinct(real_bc44)) %>%
  ungroup() %>%
  mutate(representation = barcodes/barcodes[full_sample == "POT NA"] * 100,
         reduction = 100 - representation)

#Plot Extended Figure 3d
reduction.df %>%
  mutate(sample = factor(sample, levels=c("POT", "UT", "UT R", "Cis", "Cis R", "JQ1", "JQ1 R"))) %>%
  ggplot( aes(x=sample, y=reduction)) +
  geom_boxplot(colour="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Experimental Sample")+
  ylab("Cellecta Clonal Reduction (%)")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10))

#Introduce estimated clone frequency
cell_numbers.df <- read.csv("datafiles/cell_numbers.csv", header = TRUE)
cell_numbers.df$full_sample <- str_replace(cell_numbers.df$full_sample,"_"," ")
cell_numbers.df$full_sample <- str_replace(cell_numbers.df$full_sample,"_"," ")
cell_numbers.df$full_sample <- str_replace(cell_numbers.df$full_sample,"POT","POT NA")

#Calculate estimated clone frequency (taking into account number of cells inputted for library prep)
test <- left_join(filtered.barcode.df, cell_numbers.df)
filtered.barcode.df <- test %>% mutate(estimated_clone_frequency = (barcode_percentage * number_cells_per_sample)/100) 

#Remove unnessecary columns
filtered.barcode.df <- filtered.barcode.df %>%
  select(-c(cell_number, conc, dna_ng, volume, total_dna, perc_dna_per_sample, number_cells_per_sample))

#Find top 75 clones in each sample
top.75 <- filtered.barcode.df %>% filter(sample == "POT") %>%
  arrange(desc(estimated_clone_frequency)) %>% distinct(real_bc44) %>%
  head(75)

#Identify cloens from POT sample and order these clones based on representation in POT sample
POT.perc.barcodes.df <- filtered.barcode.df %>% filter(sample == "POT") %>%
  arrange(desc(estimated_clone_frequency)) %>% distinct(real_bc44)
POT.perc.barcodes <- unique(POT.perc.barcodes.df$real_bc44)

#Plot Figure 2d Upper panle
library(ggforce)
filtered.barcode.df %>%
  filter(real_bc44 %in% top.75$real_bc44) %>%
  mutate(sample = factor(sample, levels=c("POT", "UT", "UT R", "Cis", "Cis R", "JQ1", "JQ1 R"))) %>%
  mutate(real_bc44 = factor(real_bc44, levels=POT.perc.barcodes)) %>%
  ggplot( aes(x=replicate, y=desc(real_bc44), size = estimated_clone_frequency, colour=real_bc44)) +
  geom_point(alpha=0.7) +
  facet_grid(~sample, scales = "free") +
  scale_colour_manual(values = P1437)+
  labs(y = "Cellecta Clone") +
  scale_size_continuous(range = c(0, 12), breaks = seq(100,2100,250)) +
  guides(color = FALSE, size=guide_legend(ncol=8, title="Estimated \nClone Frequency")) +
  xlab("") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 8, vjust = -1.2),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        aspect.ratio = 4,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Select subset of clones to look at more closely between 15 and 25 top represented
top.15 <- filtered.barcode.df %>%
  filter(sample == "POT") %>% arrange(desc(barcode_percentage)) %>%
  distinct(real_bc44) %>% head(15)

top.25 <- filtered.barcode.df %>%
  filter(sample == "POT") %>% arrange(desc(barcode_percentage)) %>%
  distinct(real_bc44) %>% head(25)

#Plot Figure 2d Lower panel
filtered.barcode.df %>%
  filter(real_bc44 %in% top.25$real_bc44) %>%
  filter(!(real_bc44 %in% top.15$real_bc44)) %>%
  mutate(sample = factor(sample, levels=c("POT", "UT", "UT R", "Cis", "Cis R", "JQ1", "JQ1 R"))) %>%
  mutate(real_bc44 = factor(real_bc44, levels=POT.perc.barcodes)) %>%
  ggplot( aes(x=replicate, y=desc(real_bc44), size = estimated_clone_frequency, colour=real_bc44)) +
  geom_point(alpha=0.7) +
  facet_grid(~sample, scales = "free") +
  scale_colour_manual(values = P1437)+
  labs(y = "Cellecta Clone") +
  scale_size_continuous(range = c(0, 5), breaks = seq(100,2100,250)) +
  guides(color = FALSE, size=guide_legend(ncol=8, title="Estimated \nClone Frequency")) +
  xlab("") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 8, vjust = -1.2),
        strip.text.x = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        aspect.ratio = 1.5,
        strip.background = element_rect(color="white", fill="white", linewidth=0.25, linetype="solid"))

#Isolate and plot clones of interest to visualise change of population overtime
bubble.filtered.df <- filtered.barcode.df %>% filter(real_bc44 %in% top.25$real_bc44) %>%
  filter(!(real_bc44 %in% top.15$real_bc44))

bubble.filtered.df$timepoint[bubble.filtered.df$sample == "POT"] <- "1"
bubble.filtered.df$timepoint[bubble.filtered.df$sample == "UT" | bubble.filtered.df$sample == "Cis" | bubble.filtered.df$sample == "JQ1"] <- "2"
bubble.filtered.df$timepoint[bubble.filtered.df$sample == "UT R" | bubble.filtered.df$sample == "Cis R" | bubble.filtered.df$sample == "JQ1 R"] <- "3"

bubble.filtered.df$timepoint <- as.numeric(bubble.filtered.df$timepoint)

#Isolate Untreated data
bubble.ut.df <- bubble.filtered.df %>% filter(sample == "POT" | sample == "UT" | sample == "UT R")

bubble.ut.df$replicate[bubble.ut.df$sample == "POT"] <- "A"
bubble.ut.A.df <- bubble.ut.df %>% filter(replicate == "A")
bubble.ut.df$replicate[bubble.ut.df$sample == "POT"] <- "B"
bubble.ut.B.df <- bubble.ut.df %>% filter(replicate == "B")
bubble.ut.df$replicate[bubble.ut.df$sample == "POT"] <- "C"
bubble.ut.C.df <- bubble.ut.df %>% filter(replicate == "C")
ut.df <- do.call("rbind", list(bubble.ut.A.df, bubble.ut.B.df, bubble.ut.C.df))

#Plot Figure 2e top panel
p1 <- ggplot(ut.df %>%
         mutate(sample = factor(sample, levels=c("POT", "UT", "UT R"))) %>%
         mutate(real_bc44 = factor(real_bc44, levels=POT.perc.barcodes)),
       aes(x=sample, y=estimated_clone_frequency, color=real_bc44, group=replicate)) +
  geom_line(aes(x=sample, y=estimated_clone_frequency), size=1.25) +
  geom_point(size = 2)+
  scale_y_continuous(limits=c(0,100)) +
  scale_colour_manual(values = P1437)+
  facet_wrap(~real_bc44, nrow = 1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10),
        strip.text.x = element_blank(),
        legend.position = "none")

#Isolate Cisplatin data
bubble.cis.df <- bubble.filtered.df %>% filter(sample == "POT" | sample == "Cis" | sample == "Cis R") 

bubble.cis.df$replicate[bubble.cis.df$sample == "POT"] <- "A"
bubble.cis.A.df <- bubble.cis.df %>% filter(replicate == "A")
bubble.cis.df$replicate[bubble.cis.df$sample == "POT"] <- "B"
bubble.cis.B.df <- bubble.cis.df %>% filter(replicate == "B")
bubble.cis.df$replicate[bubble.cis.df$sample == "POT"] <- "C"
bubble.cis.C.df <- bubble.cis.df %>% filter(replicate == "C")
cis.df <- do.call("rbind", list(bubble.cis.A.df, bubble.cis.B.df, bubble.cis.C.df))

#Plot Figure 2e middle panel
p2 <- ggplot(cis.df %>%
               mutate(sample = factor(sample, levels=c("POT", "Cis", "Cis R"))) %>%
               mutate(real_bc44 = factor(real_bc44, levels=POT.perc.barcodes)),
             aes(x=sample, y=estimated_clone_frequency, color=real_bc44, group=replicate)) +
  geom_line(aes(x=sample, y=estimated_clone_frequency), size=1.25) +
  geom_point(size = 2)+
  scale_y_continuous(limits=c(0,100)) +
  scale_colour_manual(values = P1437)+
  facet_wrap(~real_bc44, nrow = 1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10),
        strip.text.x = element_blank(),
        legend.position = "none")


#Isolate JQ1 data
bubble.jq1.df <- bubble.filtered.df %>%filter(sample == "POT" | sample == "JQ1" | sample == "JQ1 R") 

bubble.jq1.df$replicate[bubble.jq1.df$sample == "POT"] <- "A"
bubble.jq1.A.df <- bubble.jq1.df %>% filter(replicate == "A")
bubble.jq1.df$replicate[bubble.jq1.df$sample == "POT"] <- "B"
bubble.jq1.B.df <- bubble.jq1.df %>% filter(replicate == "B")
bubble.jq1.df$replicate[bubble.jq1.df$sample == "POT"] <- "C"
bubble.jq1.C.df <- bubble.jq1.df %>% filter(replicate == "C")
jq1.df <- do.call("rbind", list(bubble.jq1.A.df, bubble.jq1.B.df, bubble.jq1.C.df))

#Plor Figure 2e bottom panel
p3 <- ggplot(jq1.df %>%
               mutate(sample = factor(sample, levels=c("POT", "JQ1", "JQ1 R"))) %>%
               mutate(real_bc44 = factor(real_bc44, levels=POT.perc.barcodes)),
             aes(x=sample, y=estimated_clone_frequency, color=real_bc44, group=replicate)) +
  geom_line(aes(x=sample, y=estimated_clone_frequency), size=1.25) +
  geom_point(size = 2)+
  scale_y_continuous(limits=c(0,100)) +
  scale_colour_manual(values = P1437)+
  facet_wrap(~real_bc44, nrow = 1) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        title = element_text(size = 10),
        strip.text.x = element_blank(),
        legend.position = "none")

#Plot entire Figure 2e
p1 / p2 / P3
