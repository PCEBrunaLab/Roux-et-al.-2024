###NeuroPlaS###

##Differential Expression##

# Load necessary libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(survival)


# Load the Seurat object
nb.harmony.seurat<- readRDS("nb_seurat_Wbarcodes_manuscript.RDS")

# Create a new column in the metadata to define intermediate and non-intermediate cells
nb.harmony.seurat$cell_group <- ifelse(nb.harmony.seurat$AMT.state == "intermediate", "Intermediate", "Non-Intermediate")

# Set the identity of each cell to the new cell_group column
Idents(nb.harmony.seurat) <- nb.harmony.seurat$cell_group

table(nb.harmony.seurat$cell_group)

# Perform differential expression analysis between intermediate and non-intermediate cells
intermediate_vs_non_intermediate <- FindMarkers(nb.harmony.seurat, ident.1 = "Intermediate", ident.2 = "Non-Intermediate", 
                                                min.pct = 0.05, logfc.threshold = 0.2, test.use = "wilcox", only.pos = FALSE)

# Filter genes
sig_intermediate_vs_non_intermediate <- intermediate_vs_non_intermediate %>% 
  filter(p_val_adj < 0.05) %>% 
  arrange(desc(avg_log2FC))

# Save the significant intermediate markers to a CSV file
write.csv(sig_intermediate_vs_non_intermediate, file = "significant_intermediate_vs_non_intermediate_markers_MAST.csv", row.names = TRUE)


# Print some diagnostics
print(paste("Total significant genes:", nrow(sig_intermediate_vs_non_intermediate)))
print(head(sig_intermediate_vs_non_intermediate))

# Create a column to categorize genes as up-regulated, down-regulated, or non-differential
sig_intermediate_vs_non_intermediate <- sig_intermediate_vs_non_intermediate %>%
  mutate(sig = case_when(
    avg_log2FC > 1 & p_val_adj < 0.05 ~ "Up-Regulated",
    avg_log2FC < -1 & p_val_adj < 0.05 ~ "Down-Regulated",
    TRUE ~ "Non-Differential"
  ))

# Adjust p-values to avoid log(0) issues
sig_intermediate_vs_non_intermediate <- sig_intermediate_vs_non_intermediate %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, .Machine$double.eps, p_val_adj))

# Create a new color palette with the desired colors
my_colors <- c("#0000FF", "#FF0000", "#9E9E9E") # Blue, Red, Gray
names(my_colors) <- c("Down-Regulated", "Up-Regulated", "Non-Differential")

# Calculate the maximum y-value for -log10(p_val_adj)
max_y_value <- max(-log10(sig_intermediate_vs_non_intermediate$p_val_adj), na.rm = TRUE)
print(paste("Max y-value (before buffer):", max_y_value))

# Add a buffer to the maximum y-value
max_y_value <- max_y_value + 1
print(paste("Max y-value (after buffer):", max_y_value))

volcano_plot <- ggplot(sig_intermediate_vs_non_intermediate, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(size = 1.5) +
  xlab("Log2 Fold Change") +
  ylab("-log10 p-value") +
  ggtitle("Volcano Plot: Intermediate vs Non-Intermediate") +
  theme_bw() +
  scale_color_manual(values = my_colors, name = "Legend") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(1e-6), linetype = "dashed", color = "gray") +
  coord_cartesian(ylim = c(0, max_y_value)) +
  theme(legend.position = "right")

# Save the plot as PNG and PDF
ggsave("volcano_plot.png", volcano_plot, width = 8, height = 6, dpi = 300)
ggsave("volcano_plot.pdf", volcano_plot, width = 8, height = 6)


#########################################################################################################################################
##NeuroPlaS in the PDXs, PDOs and Cell line##

# M.Dyer Cell State
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)

# Function to load gene list and calculate expression scores
prepare_seurat_object <- function(seurat_obj, gene_list_path) {
  # Load the gene list from CSV file
  gene_list <- read_csv(gene_list_path)
  
  # Convert gene list to a vector
  genes <- gene_list$gene
  
  # Calculate module scores for the gene list
  seurat_obj <- AddModuleScore(seurat_obj, features = list(genes), name = "NeuroPlaS.Sig", assay = "SCT")
  
  return(seurat_obj)
}

# Function to calculate stats and generate statistical report
generate_stats_report <- function(plot_data) {
  # Calculate means, medians, and other stats for each group
  stats <- plot_data %>%
    group_by(Dyer_state) %>%
    summarise(
      mean_score = mean(NeuroPlaS_score, na.rm = TRUE),
      median_score = median(NeuroPlaS_score, na.rm = TRUE),
      sd_score = sd(NeuroPlaS_score, na.rm = TRUE),
      count = n()
    )
  
  # Calculate p-values using Wilcoxon test
  states <- unique(plot_data$Dyer_state)
  p_values <- data.frame(
    comparison = character(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(length(states) - 1)) {
    for (j in (i + 1):length(states)) {
      group1 <- plot_data$NeuroPlaS_score[plot_data$Dyer_state == states[i]]
      group2 <- plot_data$NeuroPlaS_score[plot_data$Dyer_state == states[j]]
      
      wilcox_test <- wilcox.test(group1, group2)
      p_val <- wilcox_test$p.value
      
      p_values <- rbind(p_values, data.frame(
        comparison = paste(states[i], "vs", states[j]),
        p_value = p_val
      ))
    }
  }
  
  # Adjust p-values for multiple testing using Benjamini-Hochberg method
  p_values$adjusted_p_value <- p.adjust(p_values$p_value, method = "BH")
  
  return(list(stats = stats, p_values = p_values))
}

# Function to create violin plot and calculate p-values
create_violin_plot <- function(seurat_obj, output_prefix) {
  # Extract necessary data
  plot_data <- data.frame(
    NeuroPlaS_score = seurat_obj$NeuroPlaS.Sig1,
    Dyer_state = factor(seurat_obj$Dyer.state, levels = c("ADRN", "MES", "SYMP"))
  )
  
  # Generate statistical report
  stats_report <- generate_stats_report(plot_data)
  
  # Create the violin plot
  p <- ggplot(plot_data, aes(x = Dyer_state, y = NeuroPlaS_score, fill = Dyer_state)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c("#990099", "#F37735", "#009999")) +
    scale_colour_manual(values = c("#990099", "#F37735", "#009999"), labels = c("ADRN", "MES", "SYMP")) +
    scale_x_discrete(labels = c("ADRN", "MES", "SYMP")) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5)
    ) +
    labs(title = paste("NeuroPlaS Score by Dyer State -", output_prefix), 
         x = "Dyer State", y = "NeuroPlaS Score")
  
  # Function to get significance symbols
  get_significance_symbol <- function(p_value) {
    if (p_value < 0.001) return("***")
    else if (p_value < 0.01) return("**")
    else if (p_value < 0.05) return("*")
    else return("ns")
  }
  
  # Add adjusted p-value annotations
  y_max <- max(plot_data$NeuroPlaS_score, na.rm = TRUE)
  y_range <- diff(range(plot_data$NeuroPlaS_score, na.rm = TRUE))
  
  for (i in 1:nrow(stats_report$p_values)) {
    comparison <- strsplit(stats_report$p_values$comparison[i], " vs ")[[1]]
    x1 <- which(levels(plot_data$Dyer_state) == comparison[1])
    x2 <- which(levels(plot_data$Dyer_state) == comparison[2])
    
    p <- p +
      annotate("text", x = mean(c(x1, x2)), y = y_max + (0.05 * i * y_range), 
               label = get_significance_symbol(stats_report$p_values$adjusted_p_value[i]),
               size = 5) +
      annotate("segment", x = x1, xend = x2, 
               y = y_max + (0.03 * i * y_range), yend = y_max + (0.03 * i * y_range))
  }
  
  # Save the plot
  ggsave(paste0(output_prefix, "_NeuroPlaS_Violin_Plot_Dyer_State.png"), plot = p, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_prefix, "_NeuroPlaS_Violin_Plot_Dyer_State.pdf"), plot = p, width = 8, height = 6)
  
  # Save statistical report as TXT file
  report_file <- paste0(output_prefix, "_NeuroPlaS_Statistical_Report.txt")
  sink(report_file)
  
  cat("NeuroPlaS Statistical Report\n")
  cat("============================\n\n")
  
  cat("NeuroPlaS Score Statistics:\n")
  print(stats_report$stats)
  
  cat("\nWilcoxon Test P-values (unadjusted and adjusted):\n")
  print(stats_report$p_values)
  
  sink()
  
  # Print statistical report to console
  cat("\nNeuroPlaS Score Statistics:\n")
  print(stats_report$stats)
  cat("\nWilcoxon Test P-values (unadjusted and adjusted):\n")
  print(stats_report$p_values)
  
  return(p)
}

# Main execution
# Load processed Seurat objects
nb_pdx  <- readRDS("/Users/asadr/Library/CloudStorage/OneDrive-SharedLibraries-TheInstituteofCancerResearch/The PCE Team - PCE Lab/Projects/scRNA-seq_new_manuscript_models_2024/datafiles/pdx_seurat_dyer.rds")
nb_organoids <- readRDS("/Users/asadr/Library/CloudStorage/OneDrive-SharedLibraries-TheInstituteofCancerResearch/The PCE Team - PCE Lab/Projects/scRNA-seq_new_manuscript_models_2024/datafiles/organoid_seurat_dyer.rds")
nb_cell_lines <- readRDS("/Users/asadr/Library/CloudStorage/OneDrive-SharedLibraries-TheInstituteofCancerResearch/The PCE Team - PCE Lab/Projects/scRNA-seq_new_manuscript_models_2024/datafiles/cell_lines_seurat_dyer.rds")

# Path to the gene list CSV file
gene_list_path <- "significant_intermediate_vs_non_intermediate_markers_FC1.csv"

# Prepare Seurat objects with updated NeuroPlaS scores
nb_organoids <- prepare_seurat_object(nb_organoids, gene_list_path)
nb_pdx <- prepare_seurat_object(nb_pdx, gene_list_path)
nb_cell_lines  <- prepare_seurat_object(nb_cell_lines , gene_list_path)

# Create violin plots for each model separately
plot_organoids <- create_violin_plot(nb_organoids, "Organoids")
plot_pdx <- create_violin_plot(nb_pdx, "PDX")
plot_cell_lines <- create_violin_plot(nb_cell_lines , "cell_lines")

# Save individual plots if needed
ggsave("Organoids_NeuroPlaS_Violin_Plot_Dyer_State.png", plot = plot_organoids, width = 8, height = 6, dpi = 300)
ggsave("Organoids_NeuroPlaS_Violin_Plot_Dyer_State.pdf", plot = plot_organoids, width = 8, height = 6)

ggsave("PDX_NeuroPlaS_Violin_Plot_Dyer_State.png", plot = plot_pdx, width = 8, height = 6, dpi = 300)
ggsave("PDX_NeuroPlaS_Violin_Plot_Dyer_State.pdf", plot = plot_pdx, width = 8, height = 6)

ggsave("cell_lines_NeuroPlaS_Violin_Plot_Dyer_State.png", plot = plot_cell_lines , width = 8, height = 6, dpi = 300)
ggsave("cell_lines_NeuroPlaS_Violin_Plot_Dyer_State.pdf", plot = plot_cell_lines , width = 8, height = 6)

# NeuroPlaS Analysis for AMT State

# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)

# Function to load gene list and calculate expression scores
prepare_seurat_object <- function(seurat_obj, gene_list_path) {
  gene_list <- read_csv(gene_list_path)
  genes <- gene_list$gene
  seurat_obj <- AddModuleScore(seurat_obj, features = list(genes), name = "NeuroPlaS.Sig", assay = "SCT")
  return(seurat_obj)
}

# Function to calculate stats and generate statistical report
generate_stats_report <- function(plot_data) {
  stats <- plot_data %>%
    group_by(AMT_state) %>%
    summarise(
      mean_score = mean(NeuroPlaS_score, na.rm = TRUE),
      median_score = median(NeuroPlaS_score, na.rm = TRUE),
      sd_score = sd(NeuroPlaS_score, na.rm = TRUE),
      count = n()
    )
  
  states <- unique(plot_data$AMT_state)
  p_values <- data.frame(
    comparison = character(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(length(states) - 1)) {
    for (j in (i + 1):length(states)) {
      group1 <- plot_data$NeuroPlaS_score[plot_data$AMT_state == states[i]]
      group2 <- plot_data$NeuroPlaS_score[plot_data$AMT_state == states[j]]
      
      wilcox_test <- wilcox.test(group1, group2)
      p_val <- wilcox_test$p.value
      
      p_values <- rbind(p_values, data.frame(
        comparison = paste(states[i], "vs", states[j]),
        p_value = p_val
      ))
    }
  }
  
  p_values$adjusted_p_value <- p.adjust(p_values$p_value, method = "BH")
  
  return(list(stats = stats, p_values = p_values))
}

# Function to create violin plot and calculate p-values
create_violin_plot <- function(seurat_obj, output_prefix) {
  plot_data <- data.frame(
    NeuroPlaS_score = seurat_obj$NeuroPlaS.Sig1,
    AMT_state = factor(seurat_obj$AMT.state, levels = c("ADRN", "MES", "intermediate"))
  )
  
  stats_report <- generate_stats_report(plot_data)
  
  p <- ggplot(plot_data, aes(x = AMT_state, y = NeuroPlaS_score, fill = AMT_state)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c("#990099", "#F37735", "lightgrey")) +
    scale_x_discrete(labels = c("ADRN", "MES", "intermediate")) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5)
    ) +
    labs(title = paste("NeuroPlaS Score by AMT State -", output_prefix), 
         x = "AMT State", y = "NeuroPlaS Score")
  
  get_significance_symbol <- function(p_value) {
    if (p_value < 0.001) return("***")
    else if (p_value < 0.01) return("**")
    else if (p_value < 0.05) return("*")
    else return("ns")
  }
  
  y_max <- max(plot_data$NeuroPlaS_score, na.rm = TRUE)
  y_range <- diff(range(plot_data$NeuroPlaS_score, na.rm = TRUE))
  
  for (i in 1:nrow(stats_report$p_values)) {
    comparison <- strsplit(stats_report$p_values$comparison[i], " vs ")[[1]]
    x1 <- which(levels(plot_data$AMT_state) == comparison[1])
    x2 <- which(levels(plot_data$AMT_state) == comparison[2])
    
    p <- p +
      annotate("text", x = mean(c(x1, x2)), y = y_max + (0.05 * i * y_range), 
               label = get_significance_symbol(stats_report$p_values$adjusted_p_value[i]),
               size = 5) +
      annotate("segment", x = x1, xend = x2, 
               y = y_max + (0.03 * i * y_range), yend = y_max + (0.03 * i * y_range))
  }
  
  ggsave(paste0(output_prefix, "_NeuroPlaS_Violin_Plot_AMT_State.png"), plot = p, width = 8, height = 6, dpi = 300)
  ggsave(paste0(output_prefix, "_NeuroPlaS_Violin_Plot_AMT_State.pdf"), plot = p, width = 8, height = 6)
  
  report_file <- paste0(output_prefix, "_NeuroPlaS_Statistical_Report_AMT_State.txt")
  sink(report_file)
  
  cat("NeuroPlaS Statistical Report\n")
  cat("============================\n\n")
  
  cat("NeuroPlaS Score Statistics:\n")
  print(stats_report$stats)
  
  cat("\nWilcoxon Test P-values (unadjusted and adjusted):\n")
  print(stats_report$p_values)
  
  sink()
  
  cat("\nNeuroPlaS Score Statistics:\n")
  print(stats_report$stats)
  cat("\nWilcoxon Test P-values (unadjusted and adjusted):\n")
  print(stats_report$p_values)
  
  return(p)
}

# Main execution
# Load processed Seurat objects
nb_organoids <- readRDS("/Users/asadr/Library/CloudStorage/OneDrive-SharedLibraries-TheInstituteofCancerResearch/The PCE Team - PCE Lab/Projects/scRNA-seq_new_manuscript_models_2024/datafiles/Organoids ONLY/nb_seurat_AMT_non-batch_organoids.rds")
nb_pdx <- readRDS("/Users/asadr/Library/CloudStorage/OneDrive-SharedLibraries-TheInstituteofCancerResearch/The PCE Team - PCE Lab/Projects/scRNA-seq_new_manuscript_models_2024/datafiles/PDX ONLY/nb_seurat_AMT_non-batch_pdx.rds")
nb_cell_lines <- readRDS("/Users/asadr/Library/CloudStorage/OneDrive-SharedLibraries-TheInstituteofCancerResearch/The PCE Team - PCE Lab/Projects/scRNA-seq_new_manuscript_models_2024/datafiles/cell_lines_seurat_dyer.rds")


# Path to the gene list CSV file
gene_list_path <- "significant_intermediate_vs_non_intermediate_markers_FC1.csv"

# Prepare Seurat objects with updated NeuroPlaS scores
nb_organoids <- prepare_seurat_object(nb_organoids, gene_list_path)
nb_pdx <- prepare_seurat_object(nb_pdx, gene_list_path)
nb_cell_lines <- prepare_seurat_object(nb_cell_lines, gene_list_path)

# Create violin plots for AMT state
plot_organoids <- create_violin_plot(nb_organoids, "Organoids")
plot_pdx <- create_violin_plot(nb_pdx, "PDX")
plot_cell_lines <- create_violin_plot(nb_cell_lines, "cell_lines")

# Save individual plots
ggsave("Organoids_NeuroPlaS_Violin_Plot_AMT_State.png", plot = plot_organoids, width = 8, height = 6, dpi = 300)
ggsave("Organoids_NeuroPlaS_Violin_Plot_AMT_State.pdf", plot = plot_organoids, width = 8, height = 6)

ggsave("PDX_NeuroPlaS_Violin_Plot_AMT_State.png", plot = plot_pdx, width = 8, height = 6, dpi = 300)
ggsave("PDX_NeuroPlaS_Violin_Plot_AMT_State.pdf", plot = plot_pdx, width = 8, height = 6)

ggsave("cell_lines_NeuroPlaS_Violin_Plot_AMT_State.png", plot = plot_cell_lines, width = 8, height = 6, dpi = 300)
ggsave("cell_lines_NeuroPlaS_Violin_Plot_AMT_State.pdf", plot = plot_cell_lines, width = 8, height = 6)
##########################################################################################################################
##NeuroPlaS in public dataset (GOSH and PMC)##

library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(pals)

process_seurat_object <- function(seurat_obj, object_name) {
  # Update the Seurat object
  seurat_obj <- UpdateSeuratObject(seurat_obj)
  
  # If the error persists, manually add the 'images' slot
  if (!"images" %in% slotNames(seurat_obj)) {
    seurat_obj@images <- list()
  }
  
  # Set the default assay to RNA
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Define plot themes and colors
  umap.theme <- function() {
    theme_bw() +
      theme(aspect.ratio = 1,
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(colour = "#16161D", linewidth = 0.8),
            axis.ticks = element_line(colour = "#16161D", linewidth = 0.8),
            legend.title = element_blank(),
            legend.text = element_text(size = 14),
            legend.key.size = unit(1.2, "cm"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14))
  }
  
  # Load gene lists and compute module scores
  hybrid_df <- read.csv("significant_intermediate_vs_non_intermediate_markers_FC1.csv")
  hybrid_sig <- list(intersect(hybrid_df[["gene"]], rownames(seurat_obj)))
  seurat_obj <- AddModuleScore(seurat_obj, features = hybrid_sig, assay = "RNA", seed = 12345, name = "NeuroPlaS.Sig")
  
  # Load van Groningen signatures
  vanGroningen.df <- as.data.frame(read_xlsx(path = "vanGroningen_2017.xlsx", col_names = FALSE))
  colnames(vanGroningen.df) <- c("Gene", "Signature")
  vanGroningen.sig <- lapply(unique(vanGroningen.df$Signature), function(x) {
    intersect(vanGroningen.df[vanGroningen.df$Signature == x, "Gene"], rownames(seurat_obj))
  })
  names(vanGroningen.sig) <- unique(vanGroningen.df$Signature)
  
  # Add van Groningen signatures
  seurat_obj <- AddModuleScore(seurat_obj, features = vanGroningen.sig[c("MES", "ADRN")], 
                               assay = "RNA", seed = 12345, 
                               name = c("MES.Sig", "ADRN.Sig"))
  
  # Load Dyer signatures from the Excel file
  dyer.df <- as.data.frame(read_xlsx(path = "transcription_factors.xlsx", col_names = FALSE))
  colnames(dyer.df) <- c("Gene", "Signature")
  dyer.sig <- lapply(unique(dyer.df$Signature), function(x){
    intersect(dyer.df[dyer.df$Signature == x, "Gene"], rownames(seurat_obj))
  })
  names(dyer.sig) <- unique(dyer.df$Signature)
  
  # Add Dyer signatures
  seurat_obj <- AddModuleScore(seurat_obj, features = dyer.sig[c("mesenchymal", "sympathoblast", "adrenergic")], 
                               assay = "RNA", seed = 12345, 
                               name = c("MES.Sig.Dyer", "SYM.Sig.Dyer", "ADRN.Sig.Dyer"))
  
  # Print added module score columns
  print("Added module score columns:")
  added_columns <- grep("Sig", colnames(seurat_obj@meta.data), value = TRUE)
  print(added_columns)
  
  # Function to remove duplicates and numeric suffixes
  clean_feature_names <- function(features) {
    # Remove numbers at the end of feature names
    cleaned_features <- gsub("\\d+$", "", features)
    # Remove duplicates
    unique(cleaned_features)
  }
  
  # Clean feature names
  unique_features <- clean_feature_names(added_columns)
  print("Unique features to be plotted:")
  print(unique_features)
  
  # Generate and save UMAP feature plots
  plot_feature <- function(feature, color_palette) {
    matching_columns <- grep(paste0("^", feature, "\\d*$"), colnames(seurat_obj@meta.data), value = TRUE)
    
    if (length(matching_columns) == 0) {
      warning(paste("Feature not found:", feature))
      return(NULL)
    }
    
    feature_to_plot <- matching_columns[1]
    
    FeaturePlot(seurat_obj, features = feature_to_plot, order = TRUE, min.cutoff = 0, max.cutoff = 0.5,
                pt.size = 0.5) +
      scale_colour_gradientn(colours = color_palette) +
      umap.theme() +
      ggtitle(feature) +
      theme(plot.title = element_text(size = 24, face = "bold"),
            legend.text = element_text(size = 14),
            legend.key.size = unit(1.2, "cm"))
  }
  
  # Define color palettes
  color_palettes <- list(
    MES = colorRampPalette(c("white", "#F37735"))(100),  # Orange for MES
    ADRN = colorRampPalette(c("white", "#990099"))(100),  # Purple for ADRN
    NeuroPlaS = colorRampPalette(c("white", "blue"))(100),  # Blue for NeuroPlaS
    SYM = pals::brewer.greens(100)  # Green for SYM
  )
  
  # Create plots
  plot_list <- list()
  for (feature in unique_features) {
    if (grepl("^MES", feature)) {
      color_palette <- color_palettes[["MES"]]
    } else if (grepl("^ADRN", feature)) {
      color_palette <- color_palettes[["ADRN"]]
    } else if (grepl("^NeuroPlaS", feature)) {
      color_palette <- color_palettes[["NeuroPlaS"]]
    } else if (grepl("^SYM", feature)) {
      color_palette <- color_palettes[["SYM"]]
    } else {
      next  # Skip other features
    }
    
    plot <- plot_feature(feature, color_palette)
    if (!is.null(plot)) {
      plot_list[[feature]] <- plot
    }
  }
  
  # Arrange plots
  combined_umap <- wrap_plots(plot_list, ncol = 3)
  
  # Create directories if they don't exist
  dir.create(file.path("plots", "sig", object_name), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("plots", "sig", object_name, "individual"), recursive = TRUE, showWarnings = FALSE)
  
  # Save combined plot
  ggsave(paste0("plots/sig/", object_name, "/umaps_signatures_combined.pdf"), plot = combined_umap, width = 30, height = 22)
  ggsave(paste0("plots/sig/", object_name, "/umaps_signatures_combined.png"), plot = combined_umap, width = 30, height = 22, dpi = 300)
  
  # Save individual plots
  for (feature in names(plot_list)) {
    ggsave(
      filename = paste0("plots/sig/", object_name, "/individual/umap_", feature, ".pdf"),
      plot = plot_list[[feature]],
      width = 10,
      height = 8
    )
  }
  
  # Generate and save UMAP with cell state as legend
  umap_with_cell_state <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_state", label = TRUE, pt.size = 0.5) +
    umap.theme()
  
  # Save UMAP plot
  ggsave(paste0("plots/sig/", object_name, "/umap_with_cell_state.pdf"), plot = umap_with_cell_state, width = 10, height = 8)
  ggsave(paste0("plots/sig/", object_name, "/umap_with_cell_state.png"), plot = umap_with_cell_state, width = 10, height = 8)
  
  return(seurat_obj)
}

# Process nb_GOSH
nb_GOSH <- readRDS("nb_GOSH.rds")
nb_GOSH <- process_seurat_object(nb_GOSH, "nb_GOSH")

# Process nb_PMC
nb_PMC <- readRDS("nb_PMC.rds")
nb_PMC <- process_seurat_object(nb_PMC, "nb_PMC")
####################################################################################################################################
##NeuroPlaS Analysis on DTPs##

# Updated NeuroPlaS Analysis 
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(readr)

# Function to load gene list and calculate expression scores
prepare_seurat_object <- function(seurat_obj, gene_list_path) {
  # Load the gene list from CSV file
  gene_list <- read_csv(gene_list_path)
  
  # Convert gene list to a vector
  genes <- gene_list$gene
  
  # Calculate module scores for the gene list
  seurat_obj <- AddModuleScore(seurat_obj, features = list(genes), name = "NeuroPlaS", assay = "SCT")
  
  return(seurat_obj)
}

# Function to calculate stats and generate statistical report
generate_stats_report <- function(plot_data) {
  # Calculate means, medians, and other stats for each group
  stats <- plot_data %>%
    group_by(Status) %>%
    summarise(
      mean_score = mean(NeuroPlaS_score, na.rm = TRUE),
      median_score = median(NeuroPlaS_score, na.rm = TRUE),
      sd_score = sd(NeuroPlaS_score, na.rm = TRUE),
      count = n()
    )
  
  # Calculate p-values using Wilcoxon test
  states <- unique(plot_data$Status)
  p_values <- data.frame(
    comparison = character(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:(length(states) - 1)) {
    for (j in (i + 1):length(states)) {
      group1 <- plot_data$NeuroPlaS_score[plot_data$Status == states[i]]
      group2 <- plot_data$NeuroPlaS_score[plot_data$Status == states[j]]
      
      wilcox_test <- wilcox.test(group1, group2)
      p_val <- wilcox_test$p.value
      
      p_values <- rbind(p_values, data.frame(
        comparison = paste(states[i], "vs", states[j]),
        p_value = p_val
      ))
    }
  }
  
  # Adjust p-values for multiple testing using Benjamini-Hochberg method
  p_values$adjusted_p_value <- p.adjust(p_values$p_value, method = "BH")
  
  return(list(stats = stats, p_values = p_values))
}

# Function to create violin plot and calculate p-values
create_violin_plot <- function(plot_data, output_prefix) {
  # Generate statistical report
  stats_report <- generate_stats_report(plot_data)
  
  # Set the desired order of groups
  plot_data$Status <- factor(plot_data$Status, levels = c("POT Non-Survivor", "POT Survivor", "Untreated Non-Survivor", "Untreated Survivor", "Cisplatin Survivor"))
  
  # Create the violin plot
  p <- ggplot(plot_data, aes(x = Status, y = NeuroPlaS_score, fill = Status)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = c("POT Non-Survivor" = "#Fe8690", "POT Survivor" = "#Fe0d22", 
                                 "Untreated Non-Survivor" = "#43acff", "Untreated Survivor" = "#2e78b2", 
                                 "Cisplatin Survivor" = "#9A6A20")) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(title = paste("NeuroPlaS Score by Status -", output_prefix), 
         x = "Status", y = "NeuroPlaS Score")
  
  # Save the plot
  ggsave(paste0(output_prefix, "_NeuroPlaS_Violin_Plot.png"), plot = p, width = 12, height = 8, dpi = 300)
  ggsave(paste0(output_prefix, "_NeuroPlaS_Violin_Plot.pdf"), plot = p, width = 12, height = 8)
  
  # Save statistical report as TXT file
  report_file <- paste0(output_prefix, "_NeuroPlaS_Statistical_Report.txt")
  sink(report_file)
  
  cat("NeuroPlaS Statistical Report\n")
  cat("============================\n\n")
  
  cat("NeuroPlaS Score Statistics:\n")
  print(stats_report$stats)
  
  cat("\nWilcoxon Test P-values (unadjusted and adjusted):\n")
  print(stats_report$p_values)
  
  sink()
  
  # Print statistical report to console
  cat("\nNeuroPlaS Score Statistics:\n")
  print(stats_report$stats)
  cat("\nWilcoxon Test P-values (unadjusted and adjusted):\n")
  print(stats_report$p_values)
  
  return(p)
}

# Main analysis
# Load the processed Seurat object
nb.seurat <- readRDS("nb_seurat_Wbarcodes_manuscript.RDS")

# Prepare Seurat object
gene_list_path <- "significant_intermediate_vs_non_intermediate_markers_FC1.csv"
nb.seurat <- prepare_seurat_object(nb.seurat, gene_list_path)

# Extract metadata and create plot data
all_cells.df <- nb.seurat@meta.data %>% 
  filter(Condition %in% c("Cisplatin_4weeksOFF", "POT", "Untreated")) %>% 
  filter(!is.na(real_bc44))

# Assign groups based on condition and survival status
cis_barcodes <- all_cells.df %>% 
  filter(Condition == "Cisplatin_4weeksOFF")
cis.barcodes <- unique(cis_barcodes$real_bc44)

all_cells.df$group <- case_when(
  all_cells.df$Condition == "Cisplatin_4weeksOFF" ~ "Cisplatin Survivor",
  all_cells.df$Condition == "POT" & all_cells.df$real_bc44 %in% cis.barcodes ~ "POT Survivor",
  all_cells.df$Condition == "POT" & !(all_cells.df$real_bc44 %in% cis.barcodes) ~ "POT Non-Survivor",
  all_cells.df$Condition == "Untreated" & all_cells.df$real_bc44 %in% cis.barcodes ~ "Untreated Survivor",
  all_cells.df$Condition == "Untreated" & !(all_cells.df$real_bc44 %in% cis.barcodes) ~ "Untreated Non-Survivor"
)

# Create a subset of the Seurat object for all cells
all_cells_subset <- subset(nb.seurat, cells = all_cells.df$CellID)

# Create plot data
plot_data <- data.frame(
  NeuroPlaS_score = all_cells_subset$NeuroPlaS1,
  Status = all_cells.df$group
)

# Create violin plot and generate statistical report
violin_plot <- create_violin_plot(plot_data, "All_Conditions")

# Display the plot
print(violin_plot)
####################################################################################################################################

##Survival Analysis##
library(readxl)
library(biomaRt)
library(GSVA)
library(GSEABase)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(cowplot)

# Read the middle gene list from the CSV file
middle_genes_new <- read.csv("significant_intermediate_vs_non_intermediate_markers_FC1.csv")

# Load the expression data
expr_data <- read.table("GSE62564_SEQC_NB_RNA-Seq_log2RPM.txt", header = TRUE, row.names = 1)

# Map gene symbols
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
refseq_ids <- rownames(expr_data)
genes_info <- getBM(attributes = c('refseq_mrna', 'hgnc_symbol'), 
                    filters = 'refseq_mrna', 
                    values = refseq_ids, 
                    mart = mart)

# Create a named vector for mapping, handling duplicates
refseq_to_symbol <- setNames(genes_info$hgnc_symbol, genes_info$refseq_mrna)
duplicates <- duplicated(refseq_to_symbol) | duplicated(refseq_to_symbol, fromLast = TRUE)
refseq_to_symbol[duplicates] <- make.unique(refseq_to_symbol[duplicates])

# Map RefSeq IDs to gene symbols
new_row_names <- sapply(rownames(expr_data), function(id) {
  if(id %in% names(refseq_to_symbol)) refseq_to_symbol[id] else id
})

# Assign new row names to expr_data
rownames(expr_data) <- new_row_names

# Prepare gene sets for ssGSEA
gene_set <- as.list(middle_genes_new$gene)
names(gene_set) <- middle_genes_new$gene
gene_sets <- lapply(seq_along(gene_set), function(i) {
  GeneSet(geneIds = gene_set[[i]], setName = names(gene_set)[i])
})
gmt_gene_set <- GeneSetCollection(gene_sets, setName = "HYB_Genes_Signature")

# Prepare expression matrix (using full dataset)
expr_matrix <- as.matrix(expr_data)

# Create a parameter object for ssGSEA
ssgsea_params <- ssgseaParam(exprData = expr_matrix,
                             geneSets = gmt_gene_set,
                             minSize = 1,
                             maxSize = Inf,
                             alpha = 0.25,
                             normalize = TRUE)

# Run GSVA with the parameter object
ssgsea_scores <- gsva(ssgsea_params)

# Standardize scores
set.seed(123)
standardized_scores <- scale(t(ssgsea_scores))

# Perform K-means clustering
kmeans_results <- kmeans(standardized_scores, centers = 2, nstart = 25)
cluster_assignments <- kmeans_results$cluster

# Viewing the samples in each cluster
clusters <- split(names(cluster_assignments), cluster_assignments)

# Calculate the median of the ssGSEA scores for each cluster to identify high and low groups
cluster_medians <- sapply(split(ssgsea_scores, cluster_assignments), function(cluster) median(unlist(cluster)))
print(cluster_medians)
cat("\nMedian ssGSEA Score for Cluster 1:", cluster_medians[1], "\n")
cat("Median ssGSEA Score for Cluster 2:", cluster_medians[2], "\n")

# Determine which cluster is high and which is low based on median ssGSEA score
if (cluster_medians[1] > cluster_medians[2]) {
  cat("Cluster 1 is the High group and Cluster 2 is the Low group\n")
  high_group <- clusters[[1]]
  low_group <- clusters[[2]]
  cluster_assignments <- ifelse(cluster_assignments == 1, "HYB_high", "HYB_low")
} else {
  cat("Cluster 2 is the High group and Cluster 1 is the Low group\n")
  high_group <- clusters[[2]]
  low_group <- clusters[[1]]
  cluster_assignments <- ifelse(cluster_assignments == 2, "HYB_high", "HYB_low")
}

# Printing the high and low group samples
cat("High Group Samples:\n")
print(high_group)
cat("\nLow Group Samples:\n")
print(low_group)



# Preparing for Survival Analysis
# Load necessary packages
library(tidyverse)
library(survival)
library(survminer)

# Prepare the clinical data
#clinical_data <- read.table("clinical_data.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Check structure
#str(clinical_data)


# Transpose the data so that columns become rows
#clinical_data_transposed <- t(clinical_data[-1, ])
#clinical_data_transposed <- as.data.frame(clinical_data_transposed)

# Extract sample names which is the first row of clinical_data
#sample_names <- clinical_data[1, -1]  # -1 to exclude the first column which is 'sample_name'

# Truncate sample_names to match the number of columns
#sample_names <- sample_names[1:ncol(clinical_data_transposed)]

# Assign the sample names as column names to the transposed data
#names(clinical_data_transposed) <- sample_names

# Remove the first row from the data frame
#clinical_data_transposed <- clinical_data_transposed[-1, ]

#write.table(clinical_data_transposed, file = "clinical_data_transposed.txt", sep = "\t", quote = FALSE)
clinical_data_transposed <- read.csv("clinical_data_transposed.csv")


clinical_data_clean <- clinical_data_transposed %>%
  mutate(
    age_at_diagnosis = as.numeric(sub("age at diagnosis: ", "", age_at_diagnosis)),
    MYCN_status = as.numeric(ifelse(sub("mycn status: ", "", MYCN_status) == "N/A", NA, sub("mycn status: ", "", MYCN_status))),
    risk_group = as.numeric(sub("high risk: ", "", risk_group)),
    inss_stage = as.factor(sub("inss stage: ", "", inss_stage)),
    progression = as.numeric(sub("progression: ", "", progression)),
    death_from_disease = as.numeric(sub("death from disease: ", "", death_from_disease)), 
    OS = as.numeric(OS),
    EFS = as.numeric(EFS),
    EFS_bin = as.numeric(EFS_bin)
  )

# Check the first few rows of the cleaned data
head(clinical_data_clean)

# Print the structure of the cleaned data frame
str(clinical_data_clean)

# Update 'clusters_df' with the correct cluster assignments and labels
clusters_df <- data.frame(sample_name = colnames(ssgsea_scores), cluster_group = cluster_assignments)
clusters_df$cluster_group <- factor(clusters_df$cluster_group, levels = c("HYB_low", "HYB_high"), labels = c("HYB Low", "HYB High"))

# Merge this with clinical data
clinical_data_merged <- merge(clinical_data_clean, clusters_df, by = "sample_name", all.x = TRUE)

# Check for any NA entries that may cause issues in survival analysis
print(sum(is.na(clinical_data_merged$cluster_group)))

# Create the survival object using 'OS' as the time variable
surv_object <- with(clinical_data_merged, Surv(time = OS, event = death_from_disease))

# Fit Kaplan-Meier model
km_fit <- survfit(surv_object ~ cluster_group, data = clinical_data_merged)

# Create the Kaplan-Meier plot with ggsurvplot including the risk table
km_plot <- ggsurvplot(
  km_fit,
  data = clinical_data_merged,
  risk.table = TRUE, # Include the risk table
  pval = TRUE, # Include the p-value of the log-rank test
  conf.int = FALSE, # Do not include the confidence interval
  xlab = "Time, day",
  ylab = "Survival Probability",
  ggtheme = theme_minimal() + theme(text = element_text(size = 16),  # Adjust overall text size
                                    axis.title = element_text(size = 18),  # Adjust axis title size
                                    plot.title = element_text(size = 20),  # Adjust plot title size
                                    legend.title = element_text(size = 16),  # Adjust legend title size
                                    legend.text = element_text(size = 14)),
  risk.table.height = 0.2, # Height of the risk table
  legend.title = "Cluster Group",
  palette = c("darkgreen", "gray"),
  legend.labs = c("HYB Low", "HYB High") # Make sure legends are correctly labeled
)

# Save the complete plot including the risk table - PNG Format
png(file = "kaplan_meier_plot_NeuroPlaS.png", width = 10, height = 8, units = "in", res = 300)
print(km_plot) # Print the complete ggsurvplot object
dev.off() # Close the PNG device

pdf("kaplan_meier_plot_NeuroPlaS.pdf", width = 10, height = 8)
combined_plot <- plot_grid(km_plot$plot, km_plot$table, ncol = 1, rel_heights = c(0.7, 0.3))
print(combined_plot)
dev.off()


# MYCN Status analysis with HYB clusters
clinical_data_merged_combined <- clinical_data_merged[!is.na(clinical_data_merged$MYCN_status), ]
clinical_data_merged_combined$MYCN_status_binary <- ifelse(clinical_data_merged_combined$MYCN_status == 0, "Non-amplified", "Amplified")
clinical_data_merged_combined$group_mycn <- paste(clinical_data_merged_combined$cluster_group, clinical_data_merged_combined$MYCN_status_binary, sep = "_")

surv_object_combined <- with(clinical_data_merged_combined, Surv(time = OS, event = death_from_disease))
km_fit_combined <- survfit(surv_object_combined ~ group_mycn, data = clinical_data_merged_combined)
legend_labs <- levels(factor(clinical_data_merged_combined$group_mycn))

combined_km_plot <- ggsurvplot(
  km_fit_combined,
  data = clinical_data_merged_combined,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time, day",
  ylab = "Survival Probability",
  ggtheme = theme_minimal() + theme(text = element_text(size = 14),  # Adjust overall text size
                                    axis.title = element_text(size = 16),  # Adjust axis title size
                                    plot.title = element_text(size = 16),  # Adjust plot title size
                                    legend.title = element_text(size = 14),  # Adjust legend title size
                                    legend.text = element_text(size = 13)),
  risk.table.height = 0.2,
  legend.title = "Group",
  palette = rainbow(length(legend_labs)),
  legend.labs = legend_labs
)

# Save combined Kaplan-Meier plot - PNG Format
ggsave("combined_kaplan_meier_plot_with_mycn_NeuroPlaS.png", plot = combined_km_plot$plot, width = 10, height = 8, units = "in", dpi = 300)
# Save combined Kaplan-Meier plot - PDF Format
ggsave("combined_kaplan_meier_plot_with_mycn_NeuroPlaS.pdf", plot = combined_km_plot$plot, width = 10, height = 8)


# Perform pairwise log-rank tests
pairwise_comparisons <- pairwise_survdiff(Surv(OS, death_from_disease) ~ group_mycn, data = clinical_data_merged_combined)
print(pairwise_comparisons)

# Extract p-values for all pairwise comparisons
pairwise_p_values <- pairwise_comparisons$p.value
print(pairwise_p_values)

# Multivariate Cox regression analysis Low as reference group
clinical_data_merged_combined$cluster_group <- relevel(factor(clinical_data_merged_combined$cluster_group), ref = "HYB Low")
clinical_data_merged_combined$MYCN_status <- factor(clinical_data_merged_combined$MYCN_status, levels = c(0, 1), labels = c("Non-amplified", "Amplified"))
cox_formula <- Surv(OS, death_from_disease) ~ cluster_group + MYCN_status + inss_stage
cox_model <- coxph(cox_formula, data = clinical_data_merged_combined)

# Print summary of the Cox regression model
summary(cox_model)

# Perform necessary calculations (assuming this is already done)

# Redirect output to a text file
sink("statistical_results.txt")

# Print pairwise comparisons
cat("Pairwise Comparisons using Log-Rank test:\n\n")
print(pairwise_comparisons)

# Extract and print p-values for all pairwise comparisons
cat("\nExtracted Pairwise P-Values:\n\n")
pairwise_p_values <- pairwise_comparisons$p.value
print(pairwise_p_values)

# Print summary of the Cox regression model
cat("\nCox Regression Model Summary:\n\n")
cox_formula <- Surv(OS, death_from_disease) ~ cluster_group + MYCN_status + inss_stage
cox_model <- coxph(cox_formula, data = clinical_data_merged_combined)
print(summary(cox_model))

# Stop redirecting output
sink()

# Customize and plot the forest plot for Cox regression
forest_plot <- ggforest(
  cox_model,
  data = clinical_data_merged_combined,
  main = "Multivariate Cox Regression",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 0.7,
  refLabel = "Reference",
  noDigits = 2
)

# Add p-values to the forest plot
p_values <- summary(cox_model)$coefficients[, "Pr(>|z|)"]
forest_plot <- forest_plot +
  geom_text(
    aes(x = 0.78, y = seq_len(nrow(summary(cox_model)$coefficients)), label = sprintf("%.3f", p_values)),
    hjust = 0,
    vjust = 0.5,
    size = 3,
    color = "black"
  ) +
  labs(x = "Hazard Ratio (95% CI)")

# Save the forest plot - PNG Format
ggsave("cox_regression_forest_plot_NeuroPlaS.png", plot = forest_plot, width = 8, height = 6, dpi = 300)
ggsave("cox_regression_forest_plot_NeuroPlaS.pdf", plot = forest_plot, width = 8, height = 6, dpi = 300)
# Check proportional hazards assumption
cox.zph(cox_model)

# Subgroup based on HYB cluster and MYCN
clinical_data_merged_combined$combined_group <- paste(clinical_data_merged_combined$cluster_group, clinical_data_merged_combined$MYCN_status_binary, sep = "_")

# Relevel the combined_group factor to set "HYB Low_Non-amplified" as the reference level
clinical_data_merged_combined$combined_group <- relevel(factor(clinical_data_merged_combined$combined_group), ref = "HYB Low_Non-amplified") # 

# Fit the Cox regression model 
cox_formula <- Surv(OS, death_from_disease) ~ combined_group
cox_model_combined <- coxph(cox_formula, data = clinical_data_merged_combined)

# create a forest plot 
forest_plot_combined <- ggforest(cox_model_combined, data = clinical_data_merged_combined, main = "Multivariate Cox Regression", cpositions = c(0.02, 0.22, 0.4), fontsize = 0.7,  noDigits = 2, refLabel = "") 
ggsave("cox_regression_forest_plot_NeuroPlaS.png", plot = forest_plot_combined, width = 8, height = 6, dpi = 300)
ggsave("cox_regression_forest_plot_NeuroPlaS.pdf", plot = forest_plot_combined, width = 8, height = 6, dpi = 300)

table(clinical_data_merged_combined$combined_group, clinical_data_merged_combined$death_from_disease)
