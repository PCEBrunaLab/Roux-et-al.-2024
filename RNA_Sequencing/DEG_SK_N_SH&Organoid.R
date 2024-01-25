################################################################################
# Differential Expression Analysis Pipeline for Multiple Projects
# This script performs differential expression analysis using the DESeq2 package
# in R. It has been adapted for use in two separate projects: one involving 
# organoid data and another involving cell line data. The analysis is designed 
# to identify genes that are differentially expressed under various experimental 
# conditions.
#
# Sample Metadata:
# - Metadata for each project is loaded from respective CSV files, which include 
#   treatment labels and other relevant covariates for the samples.
#
# Input Files:
# - Organoid project: 'organoids_samples.csv' and 'organoids_merged_counts.csv'
# - Cell line project: 'cellline_samples.csv' and 'cellline_merged_counts.csv'
#
# Experimental Design:
# - The design formula within the DESeqDataSetFromMatrix function should be adjusted 
#   according to the specific structure of each experiment. It is essential to ensure 
#   that the factor levels  and contrasts in the 'results()' 
#   function calls match the desired comparisons for each project.
#
# Output:
# - Diagnostic plots: Quality control plots like sample clustering, PCA plots,
#   heatmaps etc. These allow to visualize relationships between samples
#   and confirm the data meets quality thresholds
#
# - Differential Expression: Lists of genes with significant differential expression 
#   between the conditions of interest. Visualizations like volcano plots or heatmaps
#   to highlight the differences.
#
# - Functional Analysis: DOT plots visualizing pathways/ontology enrichment results 
#   from tools like GO or KEGG.
#
# Usage Notes:
# - To switch between projects, ensure that the correct metadata and count data files 
#   are loaded, and adjust the contrasts and annotations accordingly.
# - The script sections are marked with comments indicating where project-specific 
#   adjustments can be made.
#

################################################################################

# set to work directory
setwd("/Users/asadr/Documents/BulkRNAseq/Results/DEG/Ceci/NewAnalysis")

# Install packages
#required.packages <- c("DESeq2", "vsn", "DEGreport")
#if (!all(required.packages %in% installed.packages()[,"Package"])){
#  BiocManager::install(required.packages[!required.packages %in% installed.packages()[,"Package"]])
#}

# load required libraries
library(DESeq2)
library(vsn)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(DEGreport)
library(pheatmap)
library(tidyverse)
library(RColorBrewer)

# load data.
meta <- read.csv("Samples.csv", header=TRUE)
count <- read.csv("merged_counts.csv",
                            header=TRUE, row.names = 1)

# check to make sure that all rows labels in meta are columns in data
all(colnames(count) == rownames(meta))
all(colnames(count) %in% meta$Sample.Name)

# Run DESeqDataSetFromMatrix function
design <- ~ State
dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = design)
dds <- DESeq(dds)

# Pre-filtering the dataset
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

# at least 3 samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 3

# The variance stabilizing transformation and the rlog
lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
meanSdPlot(cts, ranks = FALSE)
# Set the dimensions and resolution of your image
png("meanSdPlot.png", width = 800, height = 800, res = 150)
# Close the device to save the plot to the file
dev.off()

# And for logarithm-transformed counts
log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)
# Set the dimensions and resolution of your image
png("meanSdPlot_logarithm-transformed.png", width = 800, height = 800, res = 150)
# Close the device to save the plot to the file
dev.off()

# Both vst and rlog return a DESeqTransform object which is based on the SummarizedExperiment class. 
#The transformed values are no longer counts, and are stored in the assay slot.
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

# Sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$State, vsd$Treatment, sep = " - " )
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# Set the font family for the entire R session
par(family = "serif")

# Generate your plot
png("heatmap_variance_stabilizing_transformed_3.png", width = 1000, height = 1000)
pheatmap(sampleDistMatrix, 
         col = colors, 
         fontsize = 12,  
         fontsize_row = 12,  
         fontsize_col = 12,  
         show_colnames = TRUE,  
         show_rownames = TRUE)  
dev.off()

# Reset the font family to default
par(family = "sans")

library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
png("heatmap_Poisson_Distance.png")
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$State, dds$Treatment, sep=" - " )
colnames(samplePoisDistMatrix) <- rownames(sampleDistMatrix)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)
dev.off()


# Perform PCA on Data using VST data
pca_data <- plotPCA(vsd, intgroup = c("State", "Treatment"), returnData = TRUE)

# Plot with ggplot
library(ggplot2)
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment,   shape = State)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  labs(x = "PC1", y = "PC2")

# Call ggsave on a separate line
ggsave("PCA.png", p, dpi = 300)


p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment, shape = State)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(pca_data)), hjust = 1.5, vjust = 1.5)  +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  labs(x = "PC1", y = "PC2")

ggsave("PCA_Name.png", p, dpi = 300, width = 12, height = 8, units = "in")


# Perform UMAP
vst_data <- assay(vst(dds))
umap_results <- umap::umap(t(vst_data))

# Convert UMAP results to a data frame
umap_data <- as.data.frame(umap_results$layout)

# Add in metadata and sample names
umap_data <- cbind(umap_data, colData(dds))

# Create UMAP plot with sample labels
p1 <- ggplot(umap_data, aes(x = V1, y = V2, color = Treatment, shape = State)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample.Name), size = 3, vjust = 1, hjust = 0.5) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, family = "Times"),
    axis.title = element_text(size = 14, face = "bold", family = "Times"),
    legend.text = element_text(size = 10, family = "Times"),
    legend.title = element_text(size = 12, family = "Times"),
    plot.title = element_text(family = "Times")
  ) +
  labs(x = "UMAP1", y = "UMAP2", title = "UMAP plot with sample labels")

# Save the first plot
ggsave("UMAP_with_labels.png", plot = p1, dpi = 600, width = 8, height = 8, units = "in")

# Create another UMAP plot without sample labels
p2 <- ggplot(umap_data, aes(x = V1, y = V2, color = Treatment, shape = State)) +
  geom_point(size = 3) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, family = "Times"),
    axis.title = element_text(size = 14, face = "bold", family = "Times"),
    legend.text = element_text(size = 10, family = "Times"),
    legend.title = element_text(size = 12, family = "Times"),
    plot.title = element_text(family = "Times")
  ) +
  labs(x = "UMAP1", y = "UMAP2", title = "UMAP plot")

# Save the second plot
ggsave("UMAP_without_labels.png", plot = p2, dpi = 600, width = 8, height = 8, units = "in")

###################################################################
## Differential expression analysis between groups

# Combine 'Treatment' and 'State' into a new column 'Treatment_State'
meta$Treatment_State <- factor(paste0(meta$Treatment, "_", meta$State))


# Create DESeqDataSet with the updated design
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ Treatment_State)

# Run DESeq
dds <- DESeq(dds)

rld <- rlog(dds)  # 'dds' is your DESeqDataSet
rld_counts <- assay(rld)


# Inspect the new levels of 'Treatment_State'
levels(dds$Treatment_State)

# Now you can get the results for 'JQ1_MES' vs 'UT7d_MES'
res_JQ1_MESvsUT7d_MES <- results(dds, contrast=c("Treatment_State","JQ1_MES","UT7d_MES"))

# Order by adjusted p-value
res_JQ1_MESvsUT7d_MES <- res_JQ1_MESvsUT7d_MES[order(res_JQ1_MESvsUT7d_MES$padj),]

# Print the top significant genes
head(res_JQ1_MESvsUT7d_MES)

#	the row names are ensembl gene identifiers
res_JQ1_MESvsUT7d_MES$ENSG <- row.names(res_JQ1_MESvsUT7d_MES)

#	 then use biomaRt to look up these ENSG identifiers in ensembl
library(biomaRt)
mart <- useMart(host="https://feb2014.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

#	biomaRt uses a neat function called getBM
annotation_JQ1_MESvsUT7d_MES <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters="ensembl_gene_id", values = row.names(res_JQ1_MESvsUT7d_MES), mart=mart)
#	add some simple names
names(annotation_JQ1_MESvsUT7d_MES) <- c("ENSG", "symbol", "chrom", "start", "end")

# Convert DESeqResults to data frame
res_JQ1_MESvsUT7d_MES_df <- as.data.frame(res_JQ1_MESvsUT7d_MES)

# Merge with annotation
res_JQ1_MESvsUT7d_MES_df <- merge(annotation_JQ1_MESvsUT7d_MES, res_JQ1_MESvsUT7d_MES_df, by="ENSG", all.x=F, all.y=T)

# and tidy it up a bit
res_JQ1_MESvsUT7d_MES_df <- res_JQ1_MESvsUT7d_MES_df[!is.na(res_JQ1_MESvsUT7d_MES_df$symbol),]
res_JQ1_MESvsUT7d_MES_df <- res_JQ1_MESvsUT7d_MES_df[order(res_JQ1_MESvsUT7d_MES_df$padj),]
res_JQ1_MESvsUT7d_MES_df<- res_JQ1_MESvsUT7d_MES_df[!duplicated(res_JQ1_MESvsUT7d_MES_df$symbol),]
head(res_JQ1_MESvsUT7d_MES_df)

sig_genes_res_JQ1_MESvsUT7d_MES <- subset(res_JQ1_MESvsUT7d_MES_df, padj < 0.05 & abs(log2FoldChange) > 1)
head(sig_genes_res_JQ1_MESvsUT7d_MES)

write.csv(sig_genes_res_JQ1_MESvsUT7d_MES, file = "sig_genes_res_JQ1_MESvsUT7d_MES.csv", row.names = TRUE)

# Volcano Plot
# Create a dataframe from res_ADRNvsMES and add gene names as a new column

g <- ggplot(res_JQ1_MESvsUT7d_MES_df, aes(x = log2FoldChange, y = -log10(padj), color = factor(ifelse(padj < 0.05, sign(log2FoldChange), "Nonsignificant")))) +
  geom_point(alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("darkblue", "darkred", "grey"), 
                     labels = c("Down", "Up", "Non Sig"),
                     breaks = c(-1, 1, "Nonsignificant")) +
  labs(x = "log2 fold change", y = "-log10 p-value", title = "JQ1MESvsUT7dMES") +
  theme_bw(base_size = 12) + 
  theme(legend.position = "right",
        legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

ggsave("volcano_JQ1MESvsUT7dME.png", g, dpi = 300, width = 6, height = 6)

### Correct the sample name code
# Sort p-values and select the 30th smallest value
sorted_padj <- sort(res_JQ1_MESvsUT7d_MES_df$padj, na.last = TRUE)
padj_threshold <- sorted_padj[20]

# Define threshold for high expression
log2FC_threshold <- 2  # Change this to your desired threshold for high expression
# Select genes with padj less than the threshold
sig_genes <- res_JQ1_MESvsUT7d_MES_df[res_JQ1_MESvsUT7d_MES_df$padj < padj_threshold, ]

# Create a subset for very differentially expressed genes
very_diff_genes <- res_JQ1_MESvsUT7d_MES_df[abs(res_JQ1_MESvsUT7d_MES_df$log2FoldChange) > log2FC_threshold & res_JQ1_MESvsUT7d_MES_df$padj < padj_threshold, ]

# Plot
g <- ggplot(res_JQ1_MESvsUT7d_MES_df, aes(x = log2FoldChange, y = -log10(padj), color = factor(ifelse(padj < 0.05, sign(log2FoldChange), "Nonsignificant")))) +
  geom_point(alpha = 0.6, size = 0.5) +
  geom_text(data = very_diff_genes, aes(label = symbol), vjust = -0.5, hjust = 0.5, size = 3, check_overlap = TRUE) +
  scale_color_manual(values = c("darkblue", "darkred", "grey"), 
                     labels = c("Down", "Up", "Non Sig"),
                     breaks = c(-1, 1, "Nonsignificant")) +
  labs(x = "log2 fold change", y = "-log10 p-value", title = "JQ1MESvsUT7dMES") +
  theme_bw(base_size = 12) + 
  theme(legend.position = "right",
        legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

ggsave("volcano_JQ1MESvsUT7dMES_SamoleName.png", g, dpi = 300, width = 6, height = 6)

# HeatMap Plot

# Get the Ensembl IDs for the top 50 significant genes based on padj
significant_genes_ensembl <- res_JQ1_MESvsUT7d_MES_df$ENSG[order(res_JQ1_MESvsUT7d_MES_df$padj)][1:50]

# Get the symbols for these genes
significant_genes_symbols <- res_JQ1_MESvsUT7d_MES_df$symbol[order(res_JQ1_MESvsUT7d_MES_df$padj)][1:50]

# Get the counts for the significant genes
rld_counts_sig <- rld_counts[significant_genes_ensembl, ]

all(significant_genes_ensembl %in% rownames(rld_counts))

# Replace the rownames with the gene symbols
rownames(rld_counts_sig) <- significant_genes_symbols

# Select the samples from the ADRN and MES groups
selected_columns <- which(meta$Treatment_State %in% c("JQ1_MES", "UT7d_MES"))

# Subset the count data for these columns
rld_counts_sig_selected <- rld_counts_sig[, selected_columns]

# Subset the metadata for these columns
meta_selected <- meta[selected_columns, ]

sorted_indices <- order(meta_selected$Treatment_State)
is.numeric(sorted_indices)
any(is.na(sorted_indices))
max(sorted_indices) <= ncol(rld_counts_sig_selected)

rownames(meta_selected) <- meta_selected$Sample.Name
meta_selected <- meta_selected[ colnames(rld_counts_sig_selected[, sorted_indices]), ]

# Drop unused levels
meta_selected$Treatment_State <- droplevels(meta_selected$Treatment_State)

# Create annotation dataframe
annotation <- data.frame(Treatment_State = meta_selected$Treatment_State)
rownames(annotation) <- rownames(meta_selected)

# Make sure the annotation dataframe and heatmap data have the same order
annotation <- annotation[colnames(rld_counts_sig_selected[, sorted_indices]), , drop = FALSE]

# Generate heatmap
heatmap <- pheatmap(rld_counts_sig_selected[, sorted_indices],
                    cluster_rows = T,
                    show_rownames = T,
                    annotation = annotation,
                    border_color = NA,
                    fontsize = 10,
                    scale = "row",
                    fontsize_row = 8,
                    height = 20)


png("heatmap_JQ1MESvsUT7dMES.png", width = 7, height = 7, units = "in", res = 300)
heatmap
dev.off()

# Load the necessary libraries
library(clusterProfiler)

# Create a vector of your significant genes' symbols
geneList <- sig_genes_res_JQ1_MESvsUT7d_MES$symbol

# Run the GO analysis
go <- enrichGO(gene     = geneList, 
               OrgDb    = "org.Hs.eg.db", 
               keyType  = "SYMBOL",
               ont      = "BP", # for Biological Process. Change to "CC" for Cellular Component or "MF" for Molecular Function as needed
               pAdjustMethod = "BH", # adjust method for p-value
               qvalueCutoff  = 0.05, # cutoff for adjusted p-value
               readable      = TRUE) # convert entrezID to gene names


## Output results from GO analysis to a table
cluster_summary <- data.frame(go)
# Save the GO results to a CSV file
write.csv(cluster_summary, file = "GO_results_MESvsUT7d_MES.csv")

# Open a PDF device
pdf("GO_dotplot.pdf", width = 10, height = 10)

# Create the dotplot
dotplot(go, showCategory = 20)

# Close the device
dev.off()

# Create a dotplot of the KEGG analysis results
GO_plot <- dotplot(go, showCategory = 20)

# Save the plot to a PNG file
ggplot2::ggsave(filename = "GO_dotplot.png", plot = GO_plot, width = 10, height = 10, dpi = 300)

# Load the required package
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
geneList_entrez <- mapIds(org.Hs.eg.db, keys = geneList, column = "ENTREZID", keytype = "SYMBOL")

# Remove any NAs that resulted from genes that couldn't be converted
geneList_entrez <- geneList_entrez[!is.na(geneList_entrez)]

# Run the KEGG analysis again using the Entrez IDs
kegg <- enrichKEGG(gene         = geneList_entrez, 
                   organism     = 'hsa', 
                   keyType      = 'kegg', 
                   pAdjustMethod = 'BH', 
                   qvalueCutoff = 0.05)



# Save the KEGG results to a CSV file
write.csv(as.data.frame(kegg), file = "KEGG_results_MESvsUT7d_MES.csv")

# Open a PDF device
pdf("KEGG_dotplot.pdf", width = 10, height = 10)
# Dotplot for KEGG analysis
dotplot(kegg, showCategory=20)
# Close the device
dev.off()

# Create a dotplot of the KEGG analysis results
kegg_plot <- dotplot(kegg, showCategory = 20)

# Save the plot to a PNG file
ggplot2::ggsave(filename = "kegg_dotplot.png", plot = kegg_plot, width = 10, height = 10, dpi = 300)

######################################################################################

# Now can get the results for 'JQ1_Int' vs 'UT7d_Int'
res_JQ1_IntvsUT7d_Int <- results(dds, contrast=c("Treatment_State","JQ1_Int","UT7d_Int"))

# Order by adjusted p-value
res_JQ1_IntvsUT7d_Int <- res_JQ1_IntvsUT7d_Int[order(res_JQ1_IntvsUT7d_Int$padj),]

# Print the top significant genes
head(res_JQ1_IntvsUT7d_Int)

#	the row names are ensembl gene identifiers
res_JQ1_IntvsUT7d_Int$ENSG <- row.names(res_JQ1_IntvsUT7d_Int)

#	 then use biomaRt to look up these ENSG identifiers in ensembl
library(biomaRt)
mart <- useMart(host="https://feb2014.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

#	biomaRt uses a neat function called getBM
annotation_res_JQ1_IntvsUT7d_Int <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters="ensembl_gene_id", values = row.names(res_JQ1_IntvsUT7d_Int), mart=mart)
#	add some simple names
names(annotation_res_JQ1_IntvsUT7d_Int) <- c("ENSG", "symbol", "chrom", "start", "end")

# Convert DESeqResults to data frame
res_JQ1_IntvsUT7d_Int_df <- as.data.frame(res_JQ1_IntvsUT7d_Int)

# Merge with annotation
res_JQ1_IntvsUT7d_Int_df <- merge(annotation_res_JQ1_IntvsUT7d_Int, res_JQ1_IntvsUT7d_Int_df, by="ENSG", all.x=F, all.y=T)

# and tidy it up a bit
res_JQ1_IntvsUT7d_Int_df <- res_JQ1_IntvsUT7d_Int_df[!is.na(res_JQ1_IntvsUT7d_Int_df$symbol),]
res_JQ1_IntvsUT7d_Int_df <- res_JQ1_IntvsUT7d_Int_df[order(res_JQ1_IntvsUT7d_Int_df$padj),]
res_JQ1_IntvsUT7d_Int_df<- res_JQ1_IntvsUT7d_Int_df[!duplicated(res_JQ1_IntvsUT7d_Int_df$symbol),]
head(res_JQ1_IntvsUT7d_Int_df)

sig_genes_res_JQ1_IntvsUT7d_Int <- subset(res_JQ1_IntvsUT7d_Int_df, padj < 0.05 & abs(log2FoldChange) > 1)
head(sig_genes_res_JQ1_IntvsUT7d_Int)

write.csv(sig_genes_res_JQ1_IntvsUT7d_Int, file = "sig_genes_res_JQ1_IntvsUT7d_Int.csv", row.names = TRUE)

# Volcano Plot
# Create a dataframe from res_ADRNvsMES and add gene names as a new column

g <- ggplot(res_JQ1_IntvsUT7d_Int_df, aes(x = log2FoldChange, y = -log10(padj), color = factor(ifelse(padj < 0.05, sign(log2FoldChange), "Nonsignificant")))) +
  geom_point(alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("darkblue", "darkred", "grey"), 
                     labels = c("Down", "Up", "Non Sig"),
                     breaks = c(-1, 1, "Nonsignificant")) +
  labs(x = "log2 fold change", y = "-log10 p-value", title = "JQ1IntvsUT7dInt") +
  theme_bw(base_size = 12) + 
  theme(legend.position = "right",
        legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

ggsave("volcano_JQ1IntvsUT7dInt.png", g, dpi = 300, width = 6, height = 6)

### Correct the sample name code
# Sort p-values and select the 30th smallest value
sorted_padj <- sort(res_JQ1_IntvsUT7d_Int_df$padj, na.last = TRUE)
padj_threshold <- sorted_padj[20]

# Define threshold for high expression
log2FC_threshold <- 2  # Change this to your desired threshold for high expression
# Select genes with padj less than the threshold
sig_genes <- res_JQ1_IntvsUT7d_Int_df[res_JQ1_IntvsUT7d_Int_df$padj < padj_threshold, ]

# Create a subset for very differentially expressed genes
very_diff_genes <- res_JQ1_IntvsUT7d_Int_df[abs(res_JQ1_IntvsUT7d_Int_df$log2FoldChange) > log2FC_threshold & res_JQ1_IntvsUT7d_Int_df$padj < padj_threshold, ]

# Plot
g <- ggplot(res_JQ1_IntvsUT7d_Int_df, aes(x = log2FoldChange, y = -log10(padj), color = factor(ifelse(padj < 0.05, sign(log2FoldChange), "Nonsignificant")))) +
  geom_point(alpha = 0.6, size = 0.5) +
  geom_text(data = very_diff_genes, aes(label = symbol), vjust = -0.5, hjust = 0.5, size = 3, check_overlap = TRUE) +
  scale_color_manual(values = c("darkblue", "darkred", "grey"), 
                     labels = c("Down", "Up", "Non Sig"),
                     breaks = c(-1, 1, "Nonsignificant")) +
  labs(x = "log2 fold change", y = "-log10 p-value", title = "JQ1IntvsUT7dInt") +
  theme_bw(base_size = 12) + 
  theme(legend.position = "right",
        legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

ggsave("volcano_JQ1IntvsUT7dInt_SamoleName.png", g, dpi = 300, width = 6, height = 6)

# HeatMap Plot

# Get the Ensembl IDs for the top 50 significant genes based on padj
significant_genes_ensembl <- res_JQ1_IntvsUT7d_Int_df$ENSG[order(res_JQ1_IntvsUT7d_Int_df$padj)][1:50]

# Get the symbols for these genes
significant_genes_symbols <- res_JQ1_IntvsUT7d_Int_df$symbol[order(res_JQ1_IntvsUT7d_Int_df$padj)][1:50]

# Get the counts for the significant genes
rld_counts_sig <- rld_counts[significant_genes_ensembl, ]

all(significant_genes_ensembl %in% rownames(rld_counts))

# Replace the rownames with the gene symbols
rownames(rld_counts_sig) <- significant_genes_symbols

# Select the samples from the ADRN and MES groups
selected_columns <- which(meta$Treatment_State %in% c("JQ1_Int", "UT7d_Int"))

# Subset the count data for these columns
rld_counts_sig_selected <- rld_counts_sig[, selected_columns]

# Subset the metadata for these columns
meta_selected <- meta[selected_columns, ]

sorted_indices <- order(meta_selected$Treatment_State)
is.numeric(sorted_indices)
any(is.na(sorted_indices))
max(sorted_indices) <= ncol(rld_counts_sig_selected)

rownames(meta_selected) <- meta_selected$Sample.Name
meta_selected <- meta_selected[ colnames(rld_counts_sig_selected[, sorted_indices]), ]

# Drop unused levels
meta_selected$Treatment_State <- droplevels(meta_selected$Treatment_State)

# Create annotation dataframe
annotation <- data.frame(Treatment_State = meta_selected$Treatment_State)
rownames(annotation) <- rownames(meta_selected)

# Make sure the annotation dataframe and heatmap data have the same order
annotation <- annotation[colnames(rld_counts_sig_selected[, sorted_indices]), , drop = FALSE]

# Generate heatmap
heatmap <- pheatmap(rld_counts_sig_selected[, sorted_indices],
                    cluster_rows = T,
                    show_rownames = T,
                    annotation = annotation,
                    border_color = NA,
                    fontsize = 10,
                    scale = "row",
                    fontsize_row = 8,
                    height = 20)

png("heatmap_JQ1IntvsUT7dInt.png", width = 7, height = 7, units = "in", res = 300)
heatmap
dev.off()

# Load the necessary libraries
library(clusterProfiler)

# Create a vector of your significant genes' symbols
geneList <- sig_genes_res_JQ1_IntvsUT7d_Int$symbol

# Run the GO analysis
go <- enrichGO(gene     = geneList, 
               OrgDb    = "org.Hs.eg.db", 
               keyType  = "SYMBOL",
               ont      = "BP", # for Biological Process. Change to "CC" for Cellular Component or "MF" for Molecular Function as needed
               pAdjustMethod = "BH", # adjust method for p-value
               qvalueCutoff  = 0.05, # cutoff for adjusted p-value
               readable      = TRUE) # convert entrezID to gene names


## Output results from GO analysis to a table
cluster_summary <- data.frame(go)
# Save the GO results to a CSV file
write.csv(cluster_summary, file = "GO_results_JQ1_IntvsUT7d_Int.csv")

# Open a PDF device
pdf("GO_dotplot.pdf", width = 10, height = 10)

# Create the dotplot
dotplot(go, showCategory = 20)

# Close the device
dev.off()

# Create a dotplot of the KEGG analysis results
GO_plot <- dotplot(go, showCategory = 20)

# Save the plot to a PNG file
ggplot2::ggsave(filename = "GO_dotplot.png", plot = GO_plot, width = 10, height = 10, dpi = 300)

# Load the required package
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
geneList_entrez <- mapIds(org.Hs.eg.db, keys = geneList, column = "ENTREZID", keytype = "SYMBOL")

# Remove any NAs that resulted from genes that couldn't be converted
geneList_entrez <- geneList_entrez[!is.na(geneList_entrez)]

# Run the KEGG analysis again using the Entrez IDs
kegg <- enrichKEGG(gene         = geneList_entrez, 
                   organism     = 'hsa', 
                   keyType      = 'kegg', 
                   pAdjustMethod = 'BH', 
                   qvalueCutoff = 0.05)



# Save the KEGG results to a CSV file
write.csv(as.data.frame(kegg), file = "KEGG_results_JQ1_IntvsUT7d_Int.csv")

# Open a PDF device
pdf("KEGG_dotplot.pdf", width = 10, height = 10)
# Dotplot for KEGG analysis
dotplot(kegg, showCategory=20)
# Close the device
dev.off()

# Create a dotplot of the KEGG analysis results
kegg_plot <- dotplot(kegg, showCategory = 20)

# Save the plot to a PNG file
ggplot2::ggsave(filename = "kegg_dotplot.png", plot = kegg_plot, width = 10, height = 10, dpi = 300)

#########################################################################
# Now can get the results for 'JQ1_ADRN' vs 'UT7d_ADRN'
res_JQ1_ADRNvsUT7d_ADRN <- results(dds, contrast = c("Treatment_State", "JQ1_ADRN", "UT7d_ADRN"))

# Order by adjusted p-value
res_JQ1_ADRNvsUT7d_ADRN <- res_JQ1_ADRNvsUT7d_ADRN[order(res_JQ1_ADRNvsUT7d_ADRN$padj),]

# Print the top significant genes
head(res_JQ1_ADRNvsUT7d_ADRN)

#	the row names are ensembl gene identifiers
res_JQ1_ADRNvsUT7d_ADRN$ENSG <- row.names(res_JQ1_ADRNvsUT7d_ADRN)

#	 then use biomaRt to look up these ENSG identifiers in ensembl
library(biomaRt)
mart <- useMart(host="https://feb2014.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

#	biomaRt uses a neat function called getBM
annotation_res_JQ1_ADRNvsUT7d_ADRN <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters="ensembl_gene_id", values = row.names(res_JQ1_ADRNvsUT7d_ADRN), mart=mart)
#	add some simple names
names(annotation_res_JQ1_ADRNvsUT7d_ADRN) <- c("ENSG", "symbol", "chrom", "start", "end")

# Convert DESeqResults to data frame
res_JQ1_ADRNvsUT7d_ADRN_df <- as.data.frame(res_JQ1_ADRNvsUT7d_ADRN)

# Merge with annotation
res_JQ1_ADRNvsUT7d_ADRN_df <- merge(annotation_res_JQ1_ADRNvsUT7d_ADRN, res_JQ1_ADRNvsUT7d_ADRN_df, by="ENSG", all.x=F, all.y=T)

# and tidy it up a bit
res_JQ1_ADRNvsUT7d_ADRN_df <- res_JQ1_ADRNvsUT7d_ADRN_df[!is.na(res_JQ1_ADRNvsUT7d_ADRN_df$symbol),]
res_JQ1_ADRNvsUT7d_ADRN_df <- res_JQ1_ADRNvsUT7d_ADRN_df[order(res_JQ1_ADRNvsUT7d_ADRN_df$padj),]
res_JQ1_ADRNvsUT7d_ADRN_df<- res_JQ1_ADRNvsUT7d_ADRN_df[!duplicated(res_JQ1_ADRNvsUT7d_ADRN_df$symbol),]
head(res_JQ1_ADRNvsUT7d_ADRN_df)

sig_genes_res_JQ1_ADRNvsUT7d_ADRN <- subset(res_JQ1_ADRNvsUT7d_ADRN_df, padj < 0.05 & abs(log2FoldChange) > 1)
head(sig_genes_res_JQ1_ADRNvsUT7d_ADRN)

write.csv(sig_genes_res_JQ1_ADRNvsUT7d_ADRN, file = "sig_genes_res_JQ1_ADRNvsUT7d_ADRN.csv", row.names = TRUE)

# Volcano Plot
# Create a dataframe from res_ADRNvsMES and add gene names as a new column

g <- ggplot(res_JQ1_ADRNvsUT7d_ADRN_df, aes(x = log2FoldChange, y = -log10(padj), color = factor(ifelse(padj < 0.05, sign(log2FoldChange), "Nonsignificant")))) +
  geom_point(alpha = 0.6, size = 0.5) +
  scale_color_manual(values = c("darkblue", "darkred", "grey"), 
                     labels = c("Down", "Up", "Non Sig"),
                     breaks = c(-1, 1, "Nonsignificant")) +
  labs(x = "log2 fold change", y = "-log10 p-value", title = "JQ1ADRNvsUT2ADRN") +
  theme_bw(base_size = 12) + 
  theme(legend.position = "right",
        legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

ggsave("volcano_JQ1ADRNvsUT2ADRN.png", g, dpi = 300, width = 6, height = 6)

### Correct the sample name code
# Sort p-values and select the 30th smallest value
sorted_padj <- sort(res_JQ1_ADRNvsUT7d_ADRN_df$padj, na.last = TRUE)
padj_threshold <- sorted_padj[20]

# Define threshold for high expression
log2FC_threshold <- 2  # Change this to your desired threshold for high expression
# Select genes with padj less than the threshold
sig_genes <- res_JQ1_ADRNvsUT7d_ADRN_df[res_JQ1_ADRNvsUT7d_ADRN_df$padj < padj_threshold, ]

# Create a subset for very differentially expressed genes
very_diff_genes <- res_JQ1_ADRNvsUT7d_ADRN_df[abs(res_JQ1_ADRNvsUT7d_ADRN_df$log2FoldChange) > log2FC_threshold & res_JQ1_ADRNvsUT7d_ADRN_df$padj < padj_threshold, ]

# Plot
g <- ggplot(res_JQ1_ADRNvsUT7d_ADRN_df, aes(x = log2FoldChange, y = -log10(padj), color = factor(ifelse(padj < 0.05, sign(log2FoldChange), "Nonsignificant")))) +
  geom_point(alpha = 0.6, size = 0.5) +
  geom_text(data = very_diff_genes, aes(label = symbol), vjust = -0.5, hjust = 0.5, size = 3, check_overlap = TRUE) +
  scale_color_manual(values = c("darkblue", "darkred", "grey"), 
                     labels = c("Down", "Up", "Non Sig"),
                     breaks = c(-1, 1, "Nonsignificant")) +
  labs(x = "log2 fold change", y = "-log10 p-value", title = "JQ1ADRNvsUT2ADRN") +
  theme_bw(base_size = 12) + 
  theme(legend.position = "right",
        legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

ggsave("volcano_JQ1ADRNvsUT7dADRN_SamoleName.png", g, dpi = 300, width = 6, height = 6)

# HeatMap Plot

# Get the Ensembl IDs for the top 50 significant genes based on padj
significant_genes_ensembl <- res_JQ1_ADRNvsUT7d_ADRN_df$ENSG[order(res_JQ1_ADRNvsUT7d_ADRN_df$padj)][1:50]

# Get the symbols for these genes
significant_genes_symbols <- res_JQ1_ADRNvsUT7d_ADRN_df$symbol[order(res_JQ1_ADRNvsUT7d_ADRN_df$padj)][1:50]

# Get the counts for the significant genes
rld_counts_sig <- rld_counts[significant_genes_ensembl, ]

all(significant_genes_ensembl %in% rownames(rld_counts))

# Replace the rownames with the gene symbols
rownames(rld_counts_sig) <- significant_genes_symbols

# Select the samples from the ADRN and MES groups
selected_columns <- which(meta$Treatment_State %in% c("JQ1_ADRN", "UT7d_ADRN"))

# Subset the count data for these columns
rld_counts_sig_selected <- rld_counts_sig[, selected_columns]

# Subset the metadata for these columns
meta_selected <- meta[selected_columns, ]

sorted_indices <- order(meta_selected$Treatment_State)
is.numeric(sorted_indices)
any(is.na(sorted_indices))
max(sorted_indices) <= ncol(rld_counts_sig_selected)


rownames(meta_selected) <- meta_selected$Sample.Name
meta_selected <- meta_selected[ colnames(rld_counts_sig_selected[, sorted_indices]), ]

# Drop unused levels
meta_selected$Treatment_State <- droplevels(meta_selected$Treatment_State)

# Create annotation dataframe
annotation <- data.frame(Treatment_State = meta_selected$Treatment_State)
rownames(annotation) <- rownames(meta_selected)

# Make sure the annotation dataframe and heatmap data have the same order
annotation <- annotation[colnames(rld_counts_sig_selected[, sorted_indices]), , drop = FALSE]

# Generate heatmap
heatmap <- pheatmap(rld_counts_sig_selected[, sorted_indices],
                    cluster_rows = T,
                    show_rownames = T,
                    annotation = annotation,
                    border_color = NA,
                    fontsize = 10,
                    scale = "row",
                    fontsize_row = 8,
                    height = 20)

png("heatmap_JQ1ADRNvsUT7dADRN.png", width = 7, height = 7, units = "in", res = 300)
heatmap
dev.off()

# Load the necessary libraries
library(clusterProfiler)

# Create a vector of your significant genes' symbols
geneList <- sig_genes_res_JQ1_ADRNvsUT7d_ADRN$symbol

# Run the GO analysis
go <- enrichGO(gene     = geneList, 
               OrgDb    = "org.Hs.eg.db", 
               keyType  = "SYMBOL",
               ont      = "BP", # for Biological Process. Change to "CC" for Cellular Component or "MF" for Molecular Function as needed
               pAdjustMethod = "BH", # adjust method for p-value
               qvalueCutoff  = 0.05, # cutoff for adjusted p-value
               readable      = TRUE) # convert entrezID to gene names


## Output results from GO analysis to a table
cluster_summary <- data.frame(go)
# Save the GO results to a CSV file
write.csv(cluster_summary, file = "GO_results_res_JQ1_ADRNvsUT7d_ADRN.csv")

# Open a PDF device
pdf("GO_dotplot.pdf", width = 10, height = 10)

# Create the dotplot
dotplot(go, showCategory = 20)

# Close the device
dev.off()

# Create a dotplot of the KEGG analysis results
GO_plot <- dotplot(go, showCategory = 20)

# Save the plot to a PNG file
ggplot2::ggsave(filename = "GO_dotplot.png", plot = GO_plot, width = 10, height = 10, dpi = 300)

# Load the required package
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
geneList_entrez <- mapIds(org.Hs.eg.db, keys = geneList, column = "ENTREZID", keytype = "SYMBOL")

# Remove any NAs that resulted from genes that couldn't be converted
geneList_entrez <- geneList_entrez[!is.na(geneList_entrez)]

# Run the KEGG analysis again using the Entrez IDs
kegg <- enrichKEGG(gene         = geneList_entrez, 
                   organism     = 'hsa', 
                   keyType      = 'kegg', 
                   pAdjustMethod = 'BH', 
                   qvalueCutoff = 0.05)



# Save the KEGG results to a CSV file
write.csv(as.data.frame(kegg), file = "KEGG_results_res_JQ1_ADRNvsUT7d_ADRN.csv")

# Open a PDF device
pdf("KEGG_dotplot.pdf", width = 10, height = 10)
# Dotplot for KEGG analysis
dotplot(kegg, showCategory=20)
# Close the device
dev.off()

# Create a dotplot of the KEGG analysis results
kegg_plot <- dotplot(kegg, showCategory = 20)

# Save the plot to a PNG file
ggplot2::ggsave(filename = "kegg_dotplot.png", plot = kegg_plot, width = 10, height = 10, dpi = 300)


# Save the entire workspace
save.image(file = "analysis_workspace.RData")

