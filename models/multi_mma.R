# ==============================================================================
# Parallel Multimodel Averaging For VIPBG Simulated Twins Data Study
# ==============================================================================

# Load required libraries
require(OpenMx)
require(psych)
require(parallel)
require(foreach)
require(doParallel)
require(flock)  # For file locking mechanism
source("models/miFunctions.R")
options(width=245)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Create a thread-safe queue management system
create_queue_system <- function(base_dir = "datasets") {

  # Find all data folders (those starting with "data")
  data_folders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  data_folders <- data_folders[grepl("data[a-z0-9_]+$", basename(data_folders))]

  all_files <- c()
  for (folder in data_folders) {
    # Find all simulated_data CSV files directly in each data folder
    csv_files <- list.files(folder, pattern = "^simulated_data_.*\\.csv$", full.names = TRUE)
    all_files <- c(all_files, csv_files)
  }

  # Create queue file
  queue_file <- "processing_queue.txt"
  writeLines(all_files, queue_file)

  # Create processed files log
  processed_file <- "processed_files.txt"
  if (!file.exists(processed_file)) {
    writeLines("", processed_file)
  }

  cat(sprintf("Found %d data folders with simulated datasets:\n", length(data_folders)))
  for (folder in data_folders) {
    folder_files <- list.files(folder, pattern = "^simulated_data_.*\\.csv$")
    cat(sprintf("  %s: %d files\n", basename(folder), length(folder_files)))
  }

  return(list(
    queue_file = queue_file,
    processed_file = processed_file,
    total_files = length(all_files)
  ))
}

#' Thread-safe function to get next file from queue
get_next_file <- function(queue_file, processed_file) {

  # Use file locking to ensure thread safety
  lock_file <- paste0(queue_file, ".lock")

  # Try to acquire lock with timeout
  max_attempts <- 50
  attempt <- 1

  while (attempt <= max_attempts) {
    if (!file.exists(lock_file)) {
      # Create lock file
      writeLines(paste("Locked by process", Sys.getpid(), "at", Sys.time()), lock_file)

      # Double check we got the lock
      Sys.sleep(0.01)  # Small delay
      lock_content <- readLines(lock_file, warn = FALSE)
      if (length(lock_content) > 0 && grepl(Sys.getpid(), lock_content[1])) {
        break
      }
    }

    # Wait and retry
    Sys.sleep(runif(1, 0.01, 0.1))  # Random wait to reduce collision
    attempt <- attempt + 1
  }

  if (attempt > max_attempts) {
    return(NULL)  # Could not acquire lock
  }

  # Now we have the lock, process the queue
  tryCatch({
    # Read current queue and processed files
    if (!file.exists(queue_file)) {
      return(NULL)
    }

    queue_lines <- readLines(queue_file, warn = FALSE)
    processed_lines <- readLines(processed_file, warn = FALSE)

    # Remove empty lines
    queue_lines <- queue_lines[queue_lines != ""]
    processed_lines <- processed_lines[processed_lines != ""]

    # Find unprocessed files
    remaining_files <- setdiff(queue_lines, processed_lines)

    if (length(remaining_files) == 0) {
      return(NULL)  # No more files to process
    }

    # Get the next file
    next_file <- remaining_files[1]

    # Add to processed list
    writeLines(c(processed_lines, next_file), processed_file)

    return(next_file)

  }, error = function(e) {
    cat("Error in get_next_file:", e$message, "\n")
    return(NULL)
  }, finally = {
    # Always remove lock
    if (file.exists(lock_file)) {
      file.remove(lock_file)
    }
  })
}

#' Process a single dataset file
process_single_dataset <- function(file_path, worker_id = NULL) {

  cat(sprintf("[Worker %s] Processing: %s\n",
              ifelse(is.null(worker_id), "?", worker_id),
              basename(file_path)))

  tryCatch({
    # Load and prepare data (adapted from your original code)
    allVars <- c('zyg','zyg2','sex1','sex2','p1_t1','p1_t2','p2_t1','p2_t2')
    dataTwo <- read.csv(file_path, header=T, na.strings=".", col.names=allVars)

    # Subset data by zygosity
    dataMZm <- subset(dataTwo, zyg==3,)
    dataMZf <- subset(dataTwo, zyg==1,)
    dataDZm <- subset(dataTwo, zyg==4,)
    dataDZf <- subset(dataTwo, zyg==2,)
    dataDZo <- subset(dataTwo, zyg==5,)

    # Reorder data to match model expectations
    dataMZm <- dataMZm[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]
    dataMZf <- dataMZf[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]
    dataDZm <- dataDZm[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]
    dataDZf <- dataDZf[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]
    dataDZo <- dataDZo[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]

    # Run the multimodel fitting (assuming this is defined in multimodel_fitFn.R)
    source('models/multimodel_fitFn.R')

    # Perform model averaging
    mma <- mxModelAverage(reference=c(colPra,"Prop"),
                         models=fitACE5,
                         include="onlyFree",
                         SE=NULL,
                         refAsBlock=FALSE,
                         type="AIC")

    # Create output directory structure based on input file path
    rel_path <- gsub("^datasets/", "", file_path)
    rel_path <- gsub("\\.csv$", "", rel_path)

    # Create directories for each component
    components <- c("Model-Average Estimates",
                   "Model-wise Estimates",
                   "Model-wise Sampling Variances",
                   "Akaike-Weights Table")

    for (component in components) {
      output_dir <- file.path(component, dirname(rel_path))
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      }

      # Save the component
      output_file <- file.path(output_dir, paste0(basename(rel_path), ".rds"))
      saveRDS(mma[[component]], output_file)
    }

    cat(sprintf("[Worker %s] ✅ Completed: %s\n",
                ifelse(is.null(worker_id), "?", worker_id),
                basename(file_path)))

    return(TRUE)

  }, error = function(e) {
    cat(sprintf("[Worker %s] ❌ Error processing %s: %s\n",
                ifelse(is.null(worker_id), "?", worker_id),
                basename(file_path),
                e$message))
    return(FALSE)
  })
}

#' Worker function for parallel processing
worker_function <- function(worker_id, queue_info) {

  processed_count <- 0

  while (TRUE) {
    # Get next file from queue
    next_file <- get_next_file(queue_info$queue_file, queue_info$processed_file)

    if (is.null(next_file)) {
      cat(sprintf("[Worker %d] No more files to process. Processed %d files.\n",
                  worker_id, processed_count))
      break
    }

    # Process the file
    success <- process_single_dataset(next_file, worker_id)
    if (success) {
      processed_count <- processed_count + 1
    }

    # Brief pause to prevent overwhelming the system
    Sys.sleep(0.01)
  }

  return(processed_count)
}

# ==============================================================================
# MAIN PARALLEL PROCESSING FUNCTION
# ==============================================================================

#' Main function to run parallel MMA analysis
run_parallel_mma <- function(n_cores = NULL, base_dir = "datasets") {

  # Determine number of cores
  if (is.null(n_cores)) {
    n_cores <- max(1, detectCores() - 1)  # Use n-1 cores
  }

  cat(sprintf("🚀 Starting parallel MMA analysis with %d cores\n", n_cores))

  # Create queue system
  cat("📋 Creating processing queue...\n")
  queue_info <- create_queue_system(base_dir)

  cat(sprintf("📊 Found %d total files to process\n", queue_info$total_files))

  # Start timing
  start_time <- Sys.time()

  # Set up parallel cluster
  cl <- makeCluster(n_cores)

  # Export necessary objects and functions to workers
  clusterEvalQ(cl, {
    require(OpenMx)
    require(psych)
    source("models/miFunctions.R")
    options(width=245)
  })

  clusterExport(cl, c("get_next_file", "process_single_dataset", "worker_function"))

  # Run workers in parallel
  results <- parLapply(cl, 1:n_cores, function(i) {
    worker_function(i, queue_info)
  })

  # Stop cluster
  stopCluster(cl)

  # Calculate summary statistics
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  total_processed <- sum(unlist(results))

  cat("\n" , rep("=", 50), "\n")
  cat("🎉 PARALLEL PROCESSING COMPLETE!\n")
  cat(sprintf("⏱️  Total time: %.2f minutes\n", total_time))
  cat(sprintf("📁 Files processed: %d\n", total_processed))
  cat(sprintf("⚡ Average rate: %.2f files/minute\n", total_processed / total_time))
  cat(sprintf("🖥️  Used %d cores\n", n_cores))
  cat(rep("=", 50), "\n")

  # Clean up temporary files
  if (file.exists(queue_info$queue_file)) {
    file.remove(queue_info$queue_file)
  }

  return(list(
    total_files = queue_info$total_files,
    processed_files = total_processed,
    processing_time = total_time,
    cores_used = n_cores
  ))
}

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# Run with default settings (n-1 cores)
# results <- run_parallel_mma()

# Run with specific number of cores
# results <- run_parallel_mma(n_cores = 95)

# Run with custom dataset directory
# results <- run_parallel_mma(n_cores = 50, base_dir = "my_datasets")

# ==============================================================================
# MONITORING FUNCTIONS
# ==============================================================================

#' Monitor progress during processing
monitor_progress <- function(queue_file = "processing_queue.txt",
                           processed_file = "processed_files.txt",
                           refresh_interval = 10) {

  if (!file.exists(queue_file) || !file.exists(processed_file)) {
    cat("Queue files not found. Make sure processing has started.\n")
    return()
  }

  total_files <- length(readLines(queue_file, warn = FALSE))

  cat("📊 Monitoring progress (press Ctrl+C to stop monitoring)...\n\n")

  repeat {
    tryCatch({
      processed_count <- length(readLines(processed_file, warn = FALSE)) - 1  # Subtract 1 for empty first line
      processed_count <- max(0, processed_count)

      progress <- (processed_count / total_files) * 100

      cat(sprintf("\r🔄 Progress: %d/%d (%.1f%%) completed",
                  processed_count, total_files, progress))
      flush.console()

      if (processed_count >= total_files) {
        cat("\n✅ All files processed!\n")
        break
      }

      Sys.sleep(refresh_interval)

    }, error = function(e) {
      cat("\nError monitoring progress:", e$message, "\n")
      break
    })
  }
}

# To monitor in a separate R session while processing:
# monitor_progress()
