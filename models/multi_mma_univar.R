# Define the top-level datasets folder
base_dir <- "datasets"
out_dir  <- "estimates"

# Get all subdirectories starting with "data"
subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
subdirs <- subdirs[grepl("^.*/simulations", subdirs)]

# Loop over each subdirectory
for (d in subdirs) {
  cat("Entering directory:", d, "\n")

  # Get files in this subdir starting with "simulated"
  files <- list.files(d, pattern = "^simulated", full.names = TRUE)

  # Process each file one by one
  for (f in files) {
    cat("Processing file:", f, "\n")
    allVars   <- c('zyg','zyg2','sex1','sex2','p1_t1','p1_t2','p2_t1','p2_t2')
    dataTwo   <- read.csv(f, header=T, na.strings=".", col.names=allVars)
    source('models/mma_univar.R')
    file_num <- sub(".*_(\\d+)\\.csv$", "\\1", basename(f))
    subdir_name <- basename(d)
    out_file <- file.path(out_dir, paste0(subdir_name, "_", file_num, ".csv"))
    write.csv(mma$"Model-Average Estimates", out_file, row.names = FALSE)
    # Code here completes fully before moving to the next file
  }
}
