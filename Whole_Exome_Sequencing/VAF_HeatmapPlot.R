# Load required libraries
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# Set working directory
setwd("/Users/asadr/Documents/WES/maftools/new/Primary_Results/Input_files/Paper")

# List your CSV files
csv_files <- list.files(path = "/Users/asadr/Documents/WES/maftools/new/Primary_Results/Input_files/Paper", 
                        full.names = TRUE, pattern = "*.csv")

# Assuming the file names contain the cell state information
labels <- rep(c("ADRN", "MES", "Mixed"), times = c(3, 4, 5))

# Function to read and add a cell state column
read_and_label <- function(file_path, cell_state) {
  read.csv(file_path) %>%
    mutate(Cell_State = cell_state)
}

df_list <- mapply(read_and_label, csv_files, labels, SIMPLIFY = FALSE)

# Combine all data frames into one
combined_df <- bind_rows(df_list)
write.csv(combined_df, file = "combined_data.csv", row.names = FALSE)

# Create a unique identifier for each mutation
combined_df <- combined_df %>%
  unite("Mutation_ID", Hugo_Symbol, Start_Position, End_Position, sep = "_")

# Pivot the data to have one row per Mutation_ID and one column per sample
heatmap_data <- combined_df %>%
  group_by(Mutation_ID, Tumor_Sample_Barcode) %>%
  summarise(VAF = mean(VAF, na.rm = TRUE)) %>%
  pivot_wider(names_from = Tumor_Sample_Barcode, values_from = VAF)

# Convert to matrix and set row names
heatmap_matrix <- as.matrix(heatmap_data[,-1]) 
rownames(heatmap_matrix) <- heatmap_data$Mutation_ID

# Replacing NAs with 0
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Heatmap color settings
breaks <- c(0, 0.1, 1)
colors <- c("white", "blue", "green")
color_fun <- colorRamp2(breaks, colors)

# Assuming you have a vector of column names as they appear in the heatmap:
sample_names <- colnames(heatmap_matrix)

# Assign groups based on sample naming convention:
sample_groups <- ifelse(grepl("ADRN|Mixed", sample_names), "Plastic", "fixed")

# Create a data frame or factor for the annotation:
annotation_df <- data.frame(Group = factor(sample_groups, levels = c("Fixed", "Plastic")))
annotation_colors <- list(Group = c("Plastic" = "red", "Fixed" = "white"))
ha <- HeatmapAnnotation(df = annotation_df, col = annotation_colors)

# Draw the heatmap
png("my_heatmap_filter_paper.png", width = 8, height = 10, units = 'in', res = 300)
vaf <- Heatmap(heatmap_matrix,
               top_annotation = ha,
               name = "VAF",
               col = color_fun,  
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               show_column_dend = TRUE,
               show_row_names = FALSE,
               show_row_dend = TRUE)
dev.off()

# Use pdf() instead of png() to create a PDF file
pdf("my_heatmap_filter.pdf", width = 8, height = 12)
plot(vaf)
dev.off()
