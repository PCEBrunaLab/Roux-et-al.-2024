# Set your working directory to the location of count files
setwd("/Users/asadr/Documents/BulkRNAseq/Results/filtered_counts")

# Get a list of all filtered count files
count_files <- list.files(pattern = "*.txt")

# Function to extract the desired part of each filename
get_basename <- function(filename) {
  basename <- sub(".txt", "", filename)
  return(basename)
}

# Read in all the tables into a list, selecting only the first and third columns
dfs <- lapply(count_files, function(x) {
  df <- read.table(x, header = FALSE, colClasses = c("character", "numeric"), 
                   col.names = c("Ensemble", get_basename(x)))
  return(df)
})

# Merge all dataframes by the 'Ensemble' column
merged_df <- Reduce(function(x, y) merge(x, y, by = 'Ensemble', all = TRUE), dfs)

# Replace NA values with 0
merged_df[is.na(merged_df)] <- 0

# Write the merged dataframe to a new csv file
write.csv(merged_df, file = "merged_counts.csv", row.names = FALSE)
