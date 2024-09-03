## NB manuscript violin plots

## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

## [ Load dependencies ] ----

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

## [ Plot functions ] ----

## Function to calculate stats and generate statistical report
generate.stats.report <- function(df, sample, score) {
  stats <- df %>%
    group_by(!!sym(sample)) %>%
    summarise(
      mean_score = mean(!!sym(score), na.rm = TRUE),
      median_score = median(!!sym(score), na.rm = TRUE),
      sd_score = sd(!!sym(score), na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    )
  ## Calculate p-values using Wilcoxon test
  states <- unique(df[[sample]])
  p_values <- data.frame(
    comparison = character(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  for (i in 1:(length(states) - 1)) {
    for (j in (i + 1):length(states)) {
      group1 <- df[[score]][df[[sample]] == states[i]]
      group2 <- df[[score]][df[[sample]] == states[j]]
      wilcox_test <- wilcox.test(group1, group2)
      p_val <- wilcox_test$p.value
      p_values <- rbind(p_values, data.frame(
        comparison = paste(states[i], "vs", states[j]),
        p_value = p_val
      ))
    }
  }
  ## Adjust p-values for multiple testing using Benjamini-Hochberg method
  p_values$adjusted_p_value <- p.adjust(p_values$p_value, method = "BH")
  return(list(stats = stats, p_values = p_values))
}

## Function to create violin plot and calculate p-values
create.violin.plot <- function(seurat.obj, output.prefix, name, sample, score, cols) {
  ## Extract data from Seurat object
  plot.df <- data.frame(Sample = seurat.obj@meta.data[[sample]],
                        Score = seurat.obj@meta.data[[score]])
  
  stats.report <- generate.stats.report(plot.df, "Sample", "Score")
  
  p <- ggplot(plot.df, aes(x = Sample, y = Score, fill = Sample)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = cols) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    labs(title = paste("Dyer score by sample -", name), x = "", y = name)
  # ## Function to get significance symbols
  # get.significance.symbol <- function(p_value) {
  #   if (p_value < 0.001) return("***")
  #   else if (p_value < 0.01) return("**")
  #   else if (p_value < 0.05) return("*")
  #   else return("ns")
  # }
  # ## Add adjusted p-value annotations
  # y.max <- max(plot.df$Score, na.rm = TRUE)
  # y.range <- diff(range(plot.df$Score, na.rm = TRUE))
  # for (i in 1:nrow(stats.report$p_values)) {
  #   comparison <- strsplit(stats.report$p_values$comparison[i], " vs ")[[1]]
  #   x1 <- which(levels(plot.df$Sample) == comparison[1])
  #   x2 <- which(levels(plot.df$Sample) == comparison[2])
  #   p <- p +
  #     annotate("text", x = mean(c(x1, x2)), y = y.max + (0.5 * i * y.range),
  #              label = get.significance.symbol(stats.report$p_values$adjusted_p_value[i]),
  #              size = 5) +
  #     annotate("segment", x = x1, xend = x2,
  #              y = y.max + (0.45 * i * y.range), yend = y.max + (0.45 * i * y.range))
  # }
  ## Save the plot
  ggsave(paste0("plots/dyer/", output.prefix, "_dyer_sample_violin.png"), plot = p, width = 8, height = 6, dpi = 300)
  ggsave(paste0("plots/dyer/", output.prefix, "_dyer_sample_violin.pdf"), plot = p, width = 8, height = 6)
  ## Save statistical report as TXT file
  report_file <- paste0("stats/", output.prefix, "_dyer_stats_report.txt")
  sink(report_file)
  cat("Dyer Score Statistical Report\n")
  cat("=============================\n\n")
  cat("Dyer Score Statistics:\n")
  print(stats.report$stats)
  cat("\nWilcoxon Test P-values (unadjusted and adjusted):\n")
  print(stats.report$p_values)
  sink()
  # Print statistical report to console
  cat("\nDyer Score Statistics:\n")
  print(stats.report$stats)
  cat("\nWilcoxon Test P-values (unadjusted and adjusted):\n")
  print(stats.report$p_values)
  return(p)
}

## [ NB039 violins ] ----

organoid.seurat <- readRDS("data/organoid_seurat_dyer.rds")

nb039.seurat <- subset(organoid.seurat, subset = Model == "NB039")

nb039.seurat$Sample_Type <- factor(nb039.seurat$Sample_Type,
                                   levels = c("NB039_untreated", "NB039_cisplatin", "NB039_cisplatin_recovery"))

nb039.cols <- c("NB039_untreated" = "#43ACFF",
                "NB039_cisplatin" = "#DDAA33",
                "NB039_cisplatin_recovery" = "#9A6A20")

## Create violin plots for each state by condition
mes.plot <- create.violin.plot(nb039.seurat, output.prefix = "nb039_mes", name = "NB039 MES", sample = "Sample_Type", score = "MES.Dyer", cols = nb039.cols)
adrn.plot <- create.violin.plot(nb039.seurat, output.prefix = "nb039_adrn", name = "NB039 ADRN", sample = "Sample_Type", score = "ADRN.Dyer", cols = nb039.cols)
symp.plot <- create.violin.plot(nb039.seurat, output.prefix = "nb039_symp", name = "NB039 SYMP", sample = "Sample_Type", score = "SYMP.Dyer", cols = nb039.cols)

## [ NB067 violins ] ----

organoid.seurat <- readRDS("data/organoid_seurat_dyer.rds")

nb067.seurat <- subset(organoid.seurat, subset = Model == "NB067")

nb067.seurat$Sample_Type <- factor(nb067.seurat$Sample_Type,
                                   levels = c("NB067_untreated", "NB067_cisplatin", "NB067_cisplatin_recovery"))

nb067.cols <- c("NB067_untreated" = "#43ACFF",
                "NB067_cisplatin" = "#DDAA33",
                "NB067_cisplatin_recovery" = "#9A6A20")

## Create violin plots for each state by condition
mes.plot <- create.violin.plot(nb067.seurat, output.prefix = "nb067_mes", name = "NB067 MES", sample = "Sample_Type", score = "MES.Dyer", cols = nb067.cols)
adrn.plot <- create.violin.plot(nb067.seurat, output.prefix = "nb067_adrn", name = "NB067 ADRN", sample = "Sample_Type", score = "ADRN.Dyer", cols = nb067.cols)
symp.plot <- create.violin.plot(nb067.seurat, output.prefix = "nb067_symp", name = "NB067 SYMP", sample = "Sample_Type", score = "SYMP.Dyer", cols = nb067.cols)

## [ PDX violins ] ----

pdx.seurat <- readRDS("data/pdx_seurat_dyer.rds")

pdx.seurat$Sample_Type <- factor(pdx.seurat$Sample_Type,
                                 levels = c("GRNB5_untreated", "GRNB5_cisplatin", "GRNB5_cisplatin_recovery"))

pdx.cols <- c("GRNB5_untreated" = "#43ACFF",
              "GRNB5_cisplatin" = "#DDAA33",
              "GRNB5_cisplatin_recovery" = "#9A6A20")

## Create violin plots for each state by condition
mes.plot <- create.violin.plot(pdx.seurat, output.prefix = "pdx_mes", name = "PDX MES", sample = "Sample_Type", score = "MES.Dyer", cols = pdx.cols)
adrn.plot <- create.violin.plot(pdx.seurat, output.prefix = "pdx_adrn", name = "PDX ADRN", sample = "Sample_Type", score = "ADRN.Dyer", cols = pdx.cols)
symp.plot <- create.violin.plot(pdx.seurat, output.prefix = "pdx_symp", name = "PDX SYMP", sample = "Sample_Type", score = "SYMP.Dyer", cols = pdx.cols)

amt.cols <- c("ADRN" = "#990099",
              "intermediate" = "lightgrey",
              "MES" = "#F37735")

vg.plot <- create.violin.plot(pdx.seurat, output.prefix = "pdx_vg_symp", name = "PDX van Groningen", sample = "AMT.state", score = "SYMP.Dyer", cols = amt.cols)

## [ Cell lines violins ] ----

lines.seurat <- readRDS("data/cell_lines_seurat_dyer.rds")

lines.seurat$Sample_Type <- factor(lines.seurat$Sample_Type,
                                   levels = c("SK-N-SH_untreated", "SK-N-SH_cisplatin", "SK-N-SH_cisplatin_recovery",
                                              "SH-EP_untreated", "SH-EP_cisplatin", "SH-EP_cisplatin_recovery",
                                              "SH-SY5Y_untreated", "SH-EP_cisplatin", "SH-EP_cisplatin_recovery"))

lines.cols <- c("SK-N-SH_untreated" = "#43ACFF",
                "SK-N-SH_cisplatin" = "#DDAA33",
                "SK-N-SH_cisplatin_recovery" = "#9A6A20",
                "SH-EP_untreated" = "#43ACFF",
                "SH-EP_cisplatin" = "#DDAA33",
                "SH-EP_cisplatin_recovery" = "#9A6A20",
                "SH-SY5Y_untreated" = "#43ACFF",
                "SH-SY5Y_cisplatin" = "#DDAA33",
                "SH-SY5Y_cisplatin_recovery" = "#9A6A20")

amt.cols <- c("ADRN" = "#990099",
              "intermediate" = "lightgrey",
              "MES" = "#F37735")

sknsh.seurat <- subset(lines.seurat, subset = Model == "SK-N-SH")
shep.seurat <- subset(lines.seurat, subset = Model == "SH-EP")
shsy5y.seurat <- subset(lines.seurat, subset = Model == "SH-SY5Y")

## Create violin plots for each state by condition
## SK-N-SH 
sknsh.mes.plot <- create.violin.plot(sknsh.seurat, output.prefix = "sknsh_mes", name = "SK-N-SH MES", sample = "Sample_Type", score = "MES.Dyer", cols = lines.cols)
sknsh.adrn.plot <- create.violin.plot(sknsh.seurat, output.prefix = "sknsh_adrn", name = "SK-N-SH ADRN", sample = "Sample_Type", score = "ADRN.Dyer", cols = lines.cols)
sknsh.symp.plot <- create.violin.plot(sknsh.seurat, output.prefix = "sknsh_symp", name = "SK-N-SH SYMP", sample = "Sample_Type", score = "SYMP.Dyer", cols = lines.cols)

sknsh.vg.plot <- create.violin.plot(sknsh.seurat, output.prefix = "sknsh_vg_symp", name = "SK-N-SH van Groningen", sample = "AMT.state", score = "SYMP.Dyer", cols = amt.cols)

## SH-EP
shep.mes.plot <- create.violin.plot(shep.seurat, output.prefix = "shep_mes", name = "SH-EP MES", sample = "Sample_Type", score = "MES.Dyer", cols = lines.cols)
shep.adrn.plot <- create.violin.plot(shep.seurat, output.prefix = "shep_adrn", name = "SH-EP ADRN", sample = "Sample_Type", score = "ADRN.Dyer", cols = lines.cols)
shep.symp.plot <- create.violin.plot(shep.seurat, output.prefix = "shep_symp", name = "SH-EP SYMP", sample = "Sample_Type", score = "SYMP.Dyer", cols = lines.cols)

shep.vg.plot <- create.violin.plot(shep.seurat, output.prefix = "shep_vg_symp", name = "SH-EP van Groningen", sample = "AMT.state", score = "SYMP.Dyer", cols = amt.cols)

## SH-SY5Y
shsy5y.mes.plot <- create.violin.plot(shsy5y.seurat, output.prefix = "shsy5y_mes", name = "SH-SY5Y MES", sample = "Sample_Type", score = "MES.Dyer", cols = lines.cols)
shsy5y.adrn.plot <- create.violin.plot(shsy5y.seurat, output.prefix = "shsy5y_adrn", name = "SH-SY5Y ADRN", sample = "Sample_Type", score = "ADRN.Dyer", cols = lines.cols)
shsy5y.symp.plot <- create.violin.plot(shsy5y.seurat, output.prefix = "shsy5y_symp", name = "SH-SY5Y SYMP", sample = "Sample_Type", score = "SYMP.Dyer", cols = lines.cols)

shsy5y.vg.plot <- create.violin.plot(shsy5y.seurat, output.prefix = "shsy5y_vg_symp", name = "SH-SY5Y van Groningen", sample = "AMT.state", score = "SYMP.Dyer", cols = amt.cols)

## [ Barcoded cells violins ] ----

nb.seurat <- readRDS("../Updated NB analysis 2/data/nb_seurat_dyer.rds")

nb.seurat$Sample_Type <- factor(nb.seurat$Condition,
                                levels = c("POT", "Untreated", "Cisplatin(1)_ON", "Cisplatin(2)_ON",
                                           "Cisplatin_1weeksOFF", "Cisplatin_4weeksOFF",
                                           "EZH2i_ON", "EZH2i_OFF", "JQ1_ON", "JQ1_OFF"))

nb.cols <- c("POT" = "#004488",
             "Untreated" = "#43ACFF",
             "Cisplatin(1)_ON" = "#DDAA33",
             "Cisplatin(2)_ON" = "#DDAA33",
             "Cisplatin_1weeksOFF" = "#C37F26",
             "Cisplatin_4weeksOFF" = "#9A6A20",
             "EZH2i_ON" = "#CC6677",
             "EZH2i_OFF" = "#80404A",
             "JQ1_ON" = "#44AA99",
             "JQ1_OFF" = "#2B6A60")

## Create violin plots for each state by condition
mes.plot <- create.violin.plot(nb.seurat, output.prefix = "nb_mes", name = "NB barcoded MES", sample = "Condition", score = "MES.Dyer", cols = nb.cols)
adrn.plot <- create.violin.plot(nb.seurat, output.prefix = "nb_adrn", name = "NB barcoded ADRN", sample = "Condition", score = "ADRN.Dyer", cols = nb.cols)
symp.plot <- create.violin.plot(nb.seurat, output.prefix = "nb_symp", name = "NB barcoded SYMP", sample = "Condition", score = "SYMP.Dyer", cols = nb.cols)

amt.cols <- c("ADRN" = "#990099",
              "intermediate" = "lightgrey",
              "MES" = "#F37735")

vg.plot <- create.violin.plot(nb.seurat, output.prefix = "nb_vg_symp", name = "NB barcoded van Groningen", sample = "AMT.state", score = "SYMP.Dyer", cols = amt.cols)
