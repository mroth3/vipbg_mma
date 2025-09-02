# Define the top-level datasets folder
base_dir <- "datasets"
out_dir  <- "estimates"
# Ensure the output directory exists
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Get all subdirectories starting with "data"
subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
subdirs <- subdirs[grepl("^.*/simulations", subdirs)]

# Loop over each subdirectory
# Ensure the output directory exists
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Loop over each subdirectory
for (d in subdirs) {
  cat("Entering directory:", d, "\n")

  # Get all simulated dataset files
  files <- list.files(d, pattern = "^simulated.*\\.csv$", full.names = TRUE)

  for (f in files) {
    cat("Processing file:", f, "\n")
    allVars   <- c('zyg','zyg2','sex1','sex2','p1_t1','p1_t2','p2_t1','p2_t2')
    dataTwo   <- read.csv(f, header=TRUE, na.strings=".", col.names=allVars)

    source('models/mma_univar.R')

    # Keep the original dataset filename
    dataset_name <- basename(f)
    out_file     <- file.path(out_dir, dataset_name)

    write.csv(combined_results, out_file, row.names = FALSE)
  }
}
