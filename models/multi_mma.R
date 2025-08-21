# ==============================================================================
# Parallel Multimodel Averaging For VIPBG Simulated Twins Data Study
# ==============================================================================

# Load required libraries
require(OpenMx)
require(psych)
require(parallel)
options(width=245)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Quiet run function from original mma.R script
quiet_run <- function(file_path, log_file = "error_log.txt") {
  lines <- readLines(file_path, warn = FALSE)
  total_lines <- length(lines)
  con <- textConnection(lines)

  # Open log connection
  log_con <- file(log_file, open = "wt")

  on.exit({
    close(con)
    close(log_con)
  })

  current_line <- 0

  repeat {
    expr <- tryCatch(parse(con, n = 1), error = function(e) NULL)

    # Stop if there's nothing left to run
    if (is.null(expr) || length(expr) == 0) {
      break
    }

    # Track percent complete
    attr_lines <- attr(expr, "srcref")
    if (!is.null(attr_lines)) {
      last_line <- max(sapply(attr_lines, function(x) as.integer(x[3])))
      current_line <- last_line
    } else {
      current_line <- min(total_lines, current_line + 1)
    }

    # Execute expression
    withCallingHandlers(
      tryCatch(
        eval(expr, envir = .GlobalEnv),
        error = function(e) {
          writeLines(sprintf("[%s] ERROR at lines up to %d: %s",
                             format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                             current_line,
                             e$message), log_con)
        }
      ),
      warning = function(w) {
        writeLines(sprintf("[%s] WARNING at lines up to %d: %s",
                           format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                           current_line,
                           w$message), log_con)
        invokeRestart("muffleWarning")
      }
    )

    # Break if we've reached the last line
    if (current_line >= total_lines) break
  }
}

#' Create list of all files to process
create_file_list <- function(base_dir = "datasets") {
  
  # Find all data folders (those starting with "data")
  data_folders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  data_folders <- data_folders[grepl("data[a-z0-9_]+$", basename(data_folders))]
  
  all_files <- c()
  for (folder in data_folders) {
    # Find all simulated_data CSV files directly in each data folder
    csv_files <- list.files(folder, pattern = "^simulated_data_.*\\.csv$", full.names = TRUE)
    all_files <- c(all_files, csv_files)
  }
  
  cat(sprintf("Found %d data folders with simulated datasets:\n", length(data_folders)))
  for (folder in data_folders) {
    folder_files <- list.files(folder, pattern = "^simulated_data_.*\\.csv$")
    cat(sprintf("  %s: %d files\n", basename(folder), length(folder_files)))
  }
  
  return(all_files)
}

#' Process a single dataset file
process_single_dataset <- function(file_path) {
  
  cat(sprintf("Processing: %s\n", basename(file_path)))
  
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
    
    # Set up the data environment for multimodel_fitFn.R
    assign("dataMZm", dataMZm, envir = .GlobalEnv)
    assign("dataMZf", dataMZf, envir = .GlobalEnv)
    assign("dataDZm", dataDZm, envir = .GlobalEnv)
    assign("dataDZf", dataDZf, envir = .GlobalEnv)
    assign("dataDZo", dataDZo, envir = .GlobalEnv)
    
    # Use the quiet_run function like in your original script
    quiet_run('models/multimodel_fitFn.R')
    
    # Check if the required objects exist after running multimodel_fitFn.R
    if (!exists("fitACE5") || !exists("colPra")) {
      stop("Required objects (fitACE5, colPra) not created by multimodel_fitFn.R")
    }
    
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
    components <- names(mma)  # Use actual component names from mma object
    
    for (component in components) {
      # Clean component name for directory creation
      clean_component <- gsub("[^A-Za-z0-9_-]", "_", component)
      output_dir <- file.path(clean_component, dirname(rel_path))
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Save the component
      output_file <- file.path(output_dir, paste0(basename(rel_path), ".rds"))
      saveRDS(mma[[component]], output_file)
    }
    
    cat(sprintf("✅ Completed: %s\n", basename(file_path)))
    
    return(list(success = TRUE, file = file_path))
    
  }, error = function(e) {
    cat(sprintf("❌ Error processing %s: %s\n", basename(file_path), e$message))
    
    # Write detailed error to log file
    error_log <- "parallel_processing_errors.txt"
    write(paste(Sys.time(), "- Error in", basename(file_path), ":", e$message), 
          file = error_log, append = TRUE)
    
    return(list(success = FALSE, file = file_path, error = e$message))
  })
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
  
  # Get list of all files to process
  cat("📋 Creating file list...\n")
  all_files <- create_file_list(base_dir)
  
  cat(sprintf("📊 Found %d total files to process\n", length(all_files)))
  
  if (length(all_files) == 0) {
    cat("❌ No files found to process!\n")
    return(NULL)
  }
  
  # Start timing
  start_time <- Sys.time()
  
  # Create a safer processing function that includes all dependencies
  safe_process_file <- function(file_path) {
    
    # Load all required libraries within the worker
    require(OpenMx, quietly = TRUE)
    require(psych, quietly = TRUE)
    
    # Source required files within the worker
    tryCatch({
      source("models/miFunctions.R", local = TRUE)
    }, error = function(e) {
      return(list(success = FALSE, file = file_path, 
                 error = paste("Could not source miFunctions.R:", e$message)))
    })
    
    # Define quiet_run within the worker
    quiet_run <- function(file_path, log_file = "error_log.txt") {
      lines <- readLines(file_path, warn = FALSE)
      total_lines <- length(lines)
      con <- textConnection(lines)
      
      current_line <- 0
      
      repeat {
        expr <- tryCatch(parse(con, n = 1), error = function(e) NULL)
        
        if (is.null(expr) || length(expr) == 0) {
          break
        }
        
        attr_lines <- attr(expr, "srcref")
        if (!is.null(attr_lines)) {
          last_line <- max(sapply(attr_lines, function(x) as.integer(x[3])))
          current_line <- last_line
        } else {
          current_line <- min(total_lines, current_line + 1)
        }
        
        withCallingHandlers(
          tryCatch(
            eval(expr, envir = .GlobalEnv),
            error = function(e) {
              # Just continue on errors in quiet_run
            }
          ),
          warning = function(w) {
            invokeRestart("muffleWarning")
          }
        )
        
        if (current_line >= total_lines) break
      }
      
      close(con)
    }
    
    # Now do the actual processing
    tryCatch({
      cat(sprintf("Processing: %s\n", basename(file_path)))
      
      # Load and prepare data
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
      
      # Set up the data environment
      assign("dataMZm", dataMZm, envir = .GlobalEnv)
      assign("dataMZf", dataMZf, envir = .GlobalEnv)
      assign("dataDZm", dataDZm, envir = .GlobalEnv)
      assign("dataDZf", dataDZf, envir = .GlobalEnv)
      assign("dataDZo", dataDZo, envir = .GlobalEnv)
      
      # Check if multimodel_fitFn.R exists
      if (!file.exists('models/multimodel_fitFn.R')) {
        stop("models/multimodel_fitFn.R not found")
      }
      
      # Run the multimodel fitting
      quiet_run('models/multimodel_fitFn.R')
      
      # Check if required objects exist
      if (!exists("fitACE5") || !exists("colPra")) {
        stop("Required objects (fitACE5, colPra) not created by multimodel_fitFn.R")
      }
      
      # Perform model averaging
      mma <- mxModelAverage(reference=c(colPra,"Prop"), 
                           models=fitACE5, 
                           include="onlyFree", 
                           SE=NULL, 
                           refAsBlock=FALSE, 
                           type="AIC")
      
      # Create output directory structure
      rel_path <- gsub("^datasets/", "", file_path)
      rel_path <- gsub("\\.csv$", "", rel_path)
      
      # Create directories for each component
      components <- names(mma)
      
      for (component in components) {
        clean_component <- gsub("[^A-Za-z0-9_-]", "_", component)
        output_dir <- file.path(clean_component, dirname(rel_path))
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        output_file <- file.path(output_dir, paste0(basename(rel_path), ".rds"))
        saveRDS(mma[[component]], output_file)
      }
      
      cat(sprintf("✅ Completed: %s\n", basename(file_path)))
      
      return(list(success = TRUE, file = file_path))
      
    }, error = function(e) {
      cat(sprintf("❌ Error processing %s: %s\n", basename(file_path), e$message))
      
      return(list(success = FALSE, file = file_path, error = e$message))
    })
  }
  
  # Set up parallel cluster with better error handling
  cl <- makeCluster(n_cores, outfile = "")  # outfile = "" shows worker output
  
  # Set working directory on all workers
  clusterEvalQ(cl, {
    options(width=245)
  })
  
  # Process files in parallel with better error handling
  cat("🔄 Starting parallel processing...\n")
  
  results <- tryCatch({
    parLapply(cl, all_files, safe_process_file)
  }, error = function(e) {
    cat("Error in parallel processing:", e$message, "\n")
    return(NULL)
  })
  
  # Stop cluster
  stopCluster(cl)
  
  if (is.null(results)) {
    cat("❌ Parallel processing failed\n")
    return(NULL)
  }
  
  # Calculate summary statistics
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  # Count successful and failed processing
  successful <- sum(sapply(results, function(x) x$success))
  failed <- length(results) - successful
  
  cat("\n", rep("=", 50), "\n")
  cat("🎉 PARALLEL PROCESSING COMPLETE!\n")
  cat(sprintf("⏱️  Total time: %.2f minutes\n", total_time))
  cat(sprintf("📁 Total files: %d\n", length(all_files)))
  cat(sprintf("✅ Successfully processed: %d\n", successful))
  cat(sprintf("❌ Failed: %d\n", failed))
  cat(sprintf("⚡ Average rate: %.2f files/minute\n", successful / total_time))
  cat(sprintf("🖥️  Used %d cores\n", n_cores))
  cat(rep("=", 50), "\n")
  
  if (failed > 0) {
    cat("📝 Failed files:\n")
    failed_files <- sapply(results[!sapply(results, function(x) x$success)], function(x) x$file)
    for (f in failed_files[1:min(5, length(failed_files))]) {
      cat("  ", basename(f), "\n")
    }
    if (length(failed_files) > 5) {
      cat("  ... and", length(failed_files) - 5, "more\n")
    }
  }
  
  return(list(
    total_files = length(all_files),
    successful = successful,
    failed = failed,
    processing_time = total_time,
    cores_used = n_cores,
    results = results
  ))
}

# ==============================================================================
# TEST FUNCTION FOR SINGLE FILE
# ==============================================================================

#' Test processing a single file (for debugging)
test_single_file <- function(base_dir = "datasets") {
  all_files <- create_file_list(base_dir)
  if (length(all_files) > 0) {
    cat("Testing with first file:", all_files[1], "\n")
    result <- process_single_dataset(all_files[1])
    return(result)
  } else {
    cat("No files found to test\n")
    return(NULL)
  }
}

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# Test with a single file first:
# test_result <- test_single_file()

# Run with default settings (n-1 cores)
# results <- run_parallel_mma()

# Run with specific number of cores
# results <- run_parallel_mma(n_cores = 95)
