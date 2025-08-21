#!/usr/bin/env Rscript

# ==============================================================================
# PARALLEL MMA WRAPPER - KEEPS WORKING CODE INTACT
# ==============================================================================
# Creates queue, assigns files to workers, lets them run working mma.R

library(parallel)

# Calculate cores (n-1)
n_cores <- 8
#n_cores <- parallel::detectCores() - 1
cat(sprintf("🚀 Using %d cores with file queue system\n", n_cores))

# ==============================================================================
# SIMPLE FILE QUEUE SYSTEM
# ==============================================================================

# Create queue file with all dataset paths
create_file_queue <- function() {

  # Find all datasets
  base_dir <- "datasets"
  condition_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  condition_dirs <- condition_dirs[grepl("^data", condition_dirs)]

  file_list <- list()

  for (condition_dir in condition_dirs) {
    cond_path <- file.path(base_dir, condition_dir)
    csv_files <- list.files(cond_path, pattern = "^simulated_data_\\d+\\.csv$")

    for (csv_file in csv_files) {
      file_path <- file.path(condition_dir, csv_file)
      dataset_num <- gsub("simulated_data_(\\d+)\\.csv", "\\1", csv_file)

      file_list[[length(file_list) + 1]] <- list(
        condition = condition_dir,
        dataset = as.numeric(dataset_num),
        file_path = file_path
      )
    }
  }

  # Save queue
  saveRDS(file_list, "file_queue.rds")
  cat(sprintf("📋 Created file queue with %d datasets\n", length(file_list)))

  return(length(file_list))
}

# Get next file from queue (thread-safe)
get_next_file <- function() {

  lock_file <- "file_queue.rds.lock"
  max_wait <- 30
  wait_time <- 0

  # Wait for lock
  while (file.exists(lock_file) && wait_time < max_wait) {
    Sys.sleep(0.1)
    wait_time <- wait_time + 0.1
  }

  # Create lock and get file
  tryCatch({
    writeLines("locked", lock_file)

    if (!file.exists("file_queue.rds")) {
      if (file.exists(lock_file)) file.remove(lock_file)
      return(NULL)
    }

    files <- readRDS("file_queue.rds")

    if (length(files) == 0) {
      if (file.exists(lock_file)) file.remove(lock_file)
      return(NULL)
    }

    # Get first file
    next_file <- files[[1]]

    # Remove from queue
    remaining_files <- files[-1]

    if (length(remaining_files) > 0) {
      saveRDS(remaining_files, "file_queue.rds")
    } else {
      if (file.exists("file_queue.rds")) file.remove("file_queue.rds")
    }

    # Release lock
    if (file.exists(lock_file)) file.remove(lock_file)

    return(next_file)

  }, error = function(e) {
    if (file.exists(lock_file)) file.remove(lock_file)
    return(NULL)
  })
}

# ==============================================================================
# WORKER FUNCTION
# ==============================================================================

worker_process_files <- function(worker_id) {

  results <- list()
  processed_count <- 0

  cat(sprintf("Worker %d started\n", worker_id))

  while (TRUE) {

    # Get next file from queue
    file_info <- get_next_file()

    if (is.null(file_info)) {
      cat(sprintf("Worker %d finished - processed %d files\n", worker_id, processed_count))
      break
    }

    # Process this file
    result <- process_single_file(file_info, worker_id)
    results[[length(results) + 1]] <- result
    processed_count <- processed_count + 1

    cat(sprintf("Worker %d completed: %s dataset %d\n",
                worker_id, file_info$condition, file_info$dataset))
  }

  return(results)
}

# Process a single file using the working mma.R
process_single_file <- function(file_info, worker_id) {

  tryCatch({

    # Set the file path for mma.R to use
    assign("current_file_path", file.path("simulator/binned_results", file_info$file_path), envir = .GlobalEnv)
    assign("current_condition", file_info$condition, envir = .GlobalEnv)
    assign("current_dataset", file_info$dataset, envir = .GlobalEnv)

    # Run the working mma.R script
    source("models/mma.R", local = FALSE)

    return(list(
      worker = worker_id,
      condition = file_info$condition,
      dataset = file_info$dataset,
      status = "SUCCESS",
      message = "Completed"
    ))

  }, error = function(e) {
    return(list(
      worker = worker_id,
      condition = file_info$condition,
      dataset = file_info$dataset,
      status = "ERROR",
      message = as.character(e)
    ))
  })
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

cat("📁 Scanning for datasets...\n")

# Create file queue
total_files <- create_file_queue()

if (total_files == 0) {
  stop("No datasets found!")
}

# Run parallel processing
if (n_cores > 1 && total_files > 1) {

  cat("🔄 Starting parallel processing...\n")
  start_time <- Sys.time()

  # Create cluster
  cl <- makeCluster(n_cores, type = "PSOCK")
  on.exit({
    tryCatch(stopCluster(cl), error = function(e) invisible())
    # Clean up queue files
    if (file.exists("file_queue.rds")) file.remove("file_queue.rds")
    if (file.exists("file_queue.rds.lock")) file.remove("file_queue.rds.lock")
  })

  # Export functions
  clusterExport(cl, c("worker_process_files", "get_next_file", "process_single_file"),
                envir = environment())

  # Start workers
  worker_results <- parLapply(cl, 1:n_cores, worker_process_files)

  # Flatten results
  all_results <- unlist(worker_results, recursive = FALSE)

} else {

  cat("🔄 Running sequentially...\n")
  start_time <- Sys.time()

  all_results <- worker_process_files(1)
}

end_time <- Sys.time()

# Process results
success_count <- sum(sapply(all_results, function(x) x$status == "SUCCESS"))
error_count <- sum(sapply(all_results, function(x) x$status == "ERROR"))

cat("\n📊 RESULTS SUMMARY\n")
cat("==================\n")
cat(sprintf("Total files: %d\n", total_files))
cat(sprintf("Successful: %d\n", success_count))
cat(sprintf("Errors: %d\n", error_count))
cat(sprintf("Processing time: %.2f minutes\n", as.numeric(difftime(end_time, start_time, units = "mins"))))

# Save summary
dir.create("estimates", showWarnings = FALSE)
summary_df <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    worker = x$worker,
    condition = x$condition,
    dataset = x$dataset,
    status = x$status,
    message = x$message,
    stringsAsFactors = FALSE
  )
}))

write.csv(summary_df, "estimates/processing_summary.csv", row.names = FALSE)

# Show errors if any
if (error_count > 0) {
  cat("\n❌ Sample errors:\n")
  error_results <- all_results[sapply(all_results, function(x) x$status == "ERROR")]
  for (err in head(error_results, 3)) {
    cat(sprintf("  Worker %d - %s dataset %d: %s\n", err$worker, err$condition, err$dataset, err$message))
  }
}

cat("\n✅ Complete! Results saved to estimates/ directory\n")
cat("🎉 Parallel processing finished!\n")
