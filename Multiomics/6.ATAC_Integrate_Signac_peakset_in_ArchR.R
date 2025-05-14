#Integrate Signac peakset in ArchR
#Initial R Env with correct installation of packages - ArchR_Renv
renv::activate()

#Set Working Directory and Arch R parameters
setwd("~/Desktop/sn_clones/")

addArchRThreads(threads = 1) 
addArchRGenome("hg38")

# Load in Data Files
project <- readRDS("project_clones_working/Save-ArchR-Project.rds")
atac_metadata <- getCellColData(project)

saveArchRProject(ArchRProj = project, outputDirectory = "project_clones_working/", load = TRUE)

#Load in Signac peak set
peak_get_signac <- readRDS("datafiles/combined_peaks_reduced.RDS")

## Add Signac Peak Set to ArchR Project ----
project <- addPeakSet(
  ArchRProj = project,
  peakSet = peak_get_signac,
  force = TRUE)

#Generate histogram of peak width
pdf("plots/Peak_size_histogram.pdf")
hist(width(peak_get_signac), xlab="peak size", breaks=200)
abline(v=median(width(peak_get_signac)), col="black", lty=2)
text(x=median(width(peak_get_signac)), y=10000, paste0("median peak size = ", median(width(peak_get_signac)), "bp"), pos = 4)
abline(v=mean(width(peak_get_signac)), col="blue", lty=2)
text(x=mean(width(peak_get_signac)), y=8000, paste0("mean peak size = ", round(mean(width(peak_get_signac)), digits=1), "bp"), pos = 4, col="blue")
hist(width(peak_get_signac), xlab="peak size", breaks=200, main="zoom", xlim=c(0,1000))
dev.off()

## Add Peak Calling Matrix ----
#Had some issues with H5 files being unaccessible so set H5 locking to false - not always needed
#Or set threads to 1 if error is invovling nthreads 
HDF5_USE_FILE_LOCKING=FALSE
RHDF5_USE_FILE_LOCKING=FALSE
project <- addPeakMatrix(project)
getAvailableMatrices(project) # display which matrices are available: now there should be a "PeakMatrix"

#PeakMatrix - summarises the peaks that are highly specific to cell group
#GeneScoreMatrix - summarises the gene score information (inferred by accessibility data) 
#TileMatrix - summarises genomic regions that are highly specific to a certain cell group

#Save ArchR project
saveArchRProject(ArchRProj = project, outputDirectory = "project_clones_working/", load = TRUE)


