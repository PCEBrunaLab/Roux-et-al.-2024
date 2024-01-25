# QDNAseq Analysis for Chromosomal Aberrations
# Source: https://rpubs.com/mike88bell/854300

# Set the working directory to where your BAM files are located
setwd("/Users/asadr/Documents/WGS/Bam_Files/Picard")

# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("QDNAseq")
library(QDNAseq)

# Install QDNAseq.hg38 from GitHub
if (!requireNamespace("QDNAseq.hg38", quietly = TRUE))
  devtools::install_github("asntech/QDNAseq.hg38@main")

library(QDNAseq.hg38)


# List all BAM files
bamFiles <- list.files(path=".", pattern="\\.bam$", full.names=TRUE)

# Create bin annotations 
bins <- getBinAnnotations(binSize=50, genome="hg38")

# Process each BAM file
for (bamFile in bamFiles) {
  
  cat("Processing", bamFile, "\n")
  
  # Read counts
  readCounts <- binReadCounts(bins, bamfiles=bamFile, cache=TRUE) 
  
  # Filters
  readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
  
  # Corrections 
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers) 
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  
  # Segment
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  
  # Plot
  png(sub(".bam", ".png", bamFile))
  plot(copyNumbersSegmented, main=bamFile)
  dev.off()
  
}
