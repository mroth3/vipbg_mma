# ==============================================================================
# Iterative Multimodel Averaging For VIPBG Simulated Twins Data Study
# ==============================================================================

# Load Libraries & Options
require(OpenMx)
require(psych)

# Safe OpenMx settings to prevent segfaults
mxOption(key='Number of Threads', value=1)
Sys.setenv(OMP_NUM_THREADS=1)

source("models/miFunctions.R")
options(width=245)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Create list of all dataset files to process
get_all_datasets <- function(base_dir = "datasets") {
  
  # Find all data folders (those starting with "data")
  data_folders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  data_folders <- data_folders[grepl("data[a-z0-9_]+$", basename(data_folders))]
  
  all_files <- c()
  for (folder in data_folders) {
    # Find all simulated_data CSV files directly in each data folder
    csv_files <- list.files(folder, pattern = "^simulated_data_.*\\.csv$", full.names = TRUE)
    all_files <- c(all_files, csv_files)
  }
  
  cat(sprintf("📊 Found %d data folders with simulated datasets:\n", length(data_folders)))
  for (folder in data_folders) {
    folder_files <- list.files(folder, pattern = "^simulated_data_.*\\.csv$")
    cat(sprintf("  %s: %d files\n", basename(folder), length(folder_files)))
  }
  
  cat(sprintf("📁 Total files to process: %d\n\n", length(all_files)))
  
  return(all_files)
}

#' Enhanced quiet_run function with better error handling
quiet_run <- function(file_path, log_file = "error_log.txt") {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
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

    # Execute expression with error catching
    withCallingHandlers(
      tryCatch(
        eval(expr, envir = .GlobalEnv),
        error = function(e) {
          error_msg <- sprintf("[%s] ERROR at lines up to %d: %s",
                              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                              current_line,
                              e$message)
          writeLines(error_msg, log_con)
          
          # If it's a critical OpenMx error, stop
          if (grepl("mxRun|mxModel", deparse(expr))) {
            stop(e$message)
          }
        }
      ),
      warning = function(w) {
        warning_msg <- sprintf("[%s] WARNING at lines up to %d: %s",
                              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                              current_line,
                              w$message)
        writeLines(warning_msg, log_con)
        invokeRestart("muffleWarning")
      }
    )

    if (current_line >= total_lines) break
  }
}

#' Process a single dataset and save MMA results
process_single_dataset <- function(file_path, output_base_dir = "estimates") {
  
  cat(sprintf("🔄 Processing: %s\n", basename(file_path)))
  start_time <- Sys.time()
  
  tryCatch({
    # Clear previous results to avoid contamination
    if (exists("fitACE5")) rm(fitACE5, envir = .GlobalEnv)
    if (exists("colPra")) rm(colPra, envir = .GlobalEnv)
    if (exists("mma")) rm(mma, envir = .GlobalEnv)
    
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
    
    # Assign to global environment for multimodel_fitFn.R
    assign("dataMZm", dataMZm, envir = .GlobalEnv)
    assign("dataMZf", dataMZf, envir = .GlobalEnv)
    assign("dataDZm", dataDZm, envir = .GlobalEnv)
    assign("dataDZf", dataDZf, envir = .GlobalEnv)
    assign("dataDZo", dataDZo, envir = .GlobalEnv)
    
    # Run multimodel fitting
    quiet_run('models/multimodel_fitFn.R')
    
    # Check if required objects were created
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
    
    # Create estimates directory structure
    output_dir <- file.path(output_base_dir, dirname(rel_path))
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Save the complete MMA object
    output_file <- file.path(output_dir, paste0(basename(rel_path), "_mma.rds"))
    saveRDS(mma, output_file)
    
    # Also save individual components for easier access
    components <- names(mma)
    for (i in seq_along(components)) {
      component_name <- components[i]
      clean_component <- gsub("[^A-Za-z0-9_-]", "_", component_name)
      component_file <- file.path(output_dir, paste0(basename(rel_path), "_", clean_component, ".rds"))
      saveRDS(mma[[component_name]], component_file)
    }
    
    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    cat(sprintf("✅ Completed: %s (%.2f seconds)\n", basename(file_path), processing_time))
    
    return(list(
      success = TRUE, 
      file = file_path, 
      output_file = output_file,
      processing_time = processing_time,
      components = components
    ))
    
  }, error = function(e) {
    cat(sprintf("❌ Error processing %s: %s\n", basename(file_path), e$message))
    
    # Log error details
    error_log <- "failed_datasets.txt"
    write(paste(Sys.time(), "-", file_path, ":", e$message), 
          file = error_log, append = TRUE)
    
    return(list(
      success = FALSE, 
      file = file_path, 
      error = e$message
    ))
  })
}

# ==============================================================================
# MAIN PROCESSING FUNCTION
# ==============================================================================

#' Main function to process all datasets iteratively
run_iterative_mma <- function(base_dir = "datasets", output_dir = "estimates") {
  
  cat("🚀 Starting iterative MMA analysis\n")
  cat("📂 Input directory:", base_dir, "\n")
  cat("💾 Output directory:", output_dir, "\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("📁 Created output directory:", output_dir, "\n")
  }
  
  # Get list of all files to process
  all_files <- get_all_datasets(base_dir)
  
  if (length(all_files) == 0) {
    cat("❌ No files found to process!\n")
    return(NULL)
  }
  
  # Initialize tracking variables
  start_time <- Sys.time()
  results <- list()
  successful <- 0
  failed <- 0
  
  # Process each file
  for (i in seq_along(all_files)) {
    file_path <- all_files[i]
    
    cat(sprintf("\n[%d/%d] ", i, length(all_files)))
    
    # Process the file
    result <- process_single_dataset(file_path, output_dir)
    results[[i]] <- result
    
    if (result$success) {
      successful <- successful + 1
    } else {
      failed <- failed + 1
    }
    
    # Show progress every 10 files
    if (i %% 10 == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      rate <- i / elapsed
      remaining <- (length(all_files) - i) / rate
      cat(sprintf("\n📊 Progress: %d/%d (%.1f%%) | ✅ %d success | ❌ %d failed | ⏱️ %.1f min remaining\n", 
                  i, length(all_files), (i/length(all_files))*100, 
                  successful, failed, remaining))
    }
    
    # Force garbage collection every 50 files to prevent memory issues
    if (i %% 50 == 0) {
      gc()
    }
  }
  
  # Final summary
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("\n", rep("=", 60), "\n")
  cat("🎉 ITERATIVE PROCESSING COMPLETE!\n")
  cat(sprintf("⏱️  Total time: %.2f minutes\n", total_time))
  cat(sprintf("📁 Total files: %d\n", length(all_files)))
  cat(sprintf("✅ Successfully processed: %d\n", successful))
  cat(sprintf("❌ Failed: %d\n", failed))
  cat(sprintf("📈 Success rate: %.1f%%\n", (successful/length(all_files))*100))
  cat(sprintf("⚡ Average rate: %.2f files/minute\n", length(all_files) / total_time))
  cat(sprintf("💾 Results saved to: %s\n", output_dir))
  cat(rep("=", 60), "\n")
  
  if (failed > 0) {
    cat("📝 Check 'failed_datasets.txt' for error details\n")
    failed_files <- sapply(results[!sapply(results, function(x) x$success)], function(x) basename(x$file))
    cat("Failed files:\n")
    for (f in failed_files[1:min(5, length(failed_files))]) {
      cat("  ", f, "\n")
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
    success_rate = (successful/length(all_files))*100,
    results = results
  ))
}

# ==============================================================================
# USAGE
# ==============================================================================

# Run the iterative processing
# results <- run_iterative_mma()

# To run with custom directories:
# results <- run_iterative_mma(base_dir = "my_datasets", output_dir = "my_estimates")

# To check what files would be processed without running:
# files <- get_all_datasets()

cat("📋 Script loaded. Run with: results <- run_iterative_mma()\n")
