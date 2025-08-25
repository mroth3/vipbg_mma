#!/usr/bin/env Rscript

# ==============================================================================
# BATCH ACDE SIMULATION RUNNER
# ==============================================================================
# Calls engine.R to generate multiple datasets and organizes them into
# timestamped directories with comprehensive parameter logging.

library(lubridate)

# ==============================================================================
# COMMAND LINE ARGUMENT PARSING
# ==============================================================================

parse_batch_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  params <- list(
    n_datasets = 50,          # Number of datasets to generate
    n_pairs = 400,            # Twin pairs per dataset
    engine_path = "simulator/engine.R", # Path to the engine script
    base_dir = "datasets",    # Base directory for outputs
    seed_start = NULL,        # Starting seed (if NULL, uses random)
    max_correlation = 0.8,    # Maximum correlation magnitude
    stratified_prop = 0.6,    # Proportion using stratified sampling
    n_pool = 20,              # Trait pool size
    rebalance = TRUE,         # Enable correlation rebalancing
    verbose = TRUE,           # Show progress
    run_id = NULL,            # Custom run identifier (if NULL, uses timestamp)

    # NEW: ACDE component variance ranges
    a_range = c(0.01, 0.64),  # Default A component variance range (0.1^2 to 0.8^2)
    c_range = c(0.00, 0.49),  # Default C component variance range (0^2 to 0.7^2)
    d_range = c(0.00, 0.36),  # Default D component variance range (0^2 to 0.6^2)
    e_range = c(0.04, 0.81),  # Default E component variance range (0.2^2 to 0.9^2)

    # NEW: Variance and correlation constraints
    min_variance = 0.01,      # Minimum variance threshold for any component
    allow_negative_rho_D = FALSE  # Allow negative dominance correlations
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--n_datasets" && i < length(args)) {
      params$n_datasets <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--n_pairs" && i < length(args)) {
      params$n_pairs <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--engine_path" && i < length(args)) {
      params$engine_path <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--base_dir" && i < length(args)) {
      params$base_dir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--seed_start" && i < length(args)) {
      params$seed_start <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--max_corr" && i < length(args)) {
      params$max_correlation <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--stratified_prop" && i < length(args)) {
      params$stratified_prop <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--n_pool" && i < length(args)) {
      params$n_pool <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--run_id" && i < length(args)) {
      params$run_id <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--a_range" && i + 1 < length(args)) {
      params$a_range <- c(as.numeric(args[i + 1]), as.numeric(args[i + 2]))
      i <- i + 3
    } else if (args[i] == "--c_range" && i + 1 < length(args)) {
      params$c_range <- c(as.numeric(args[i + 1]), as.numeric(args[i + 2]))
      i <- i + 3
    } else if (args[i] == "--d_range" && i + 1 < length(args)) {
      params$d_range <- c(as.numeric(args[i + 1]), as.numeric(args[i + 2]))
      i <- i + 3
    } else if (args[i] == "--e_range" && i + 1 < length(args)) {
      params$e_range <- c(as.numeric(args[i + 1]), as.numeric(args[i + 2]))
      i <- i + 3
    } else if (args[i] == "--min_variance" && i < length(args)) {
      params$min_variance <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--allow_negative_rho_D" && i < length(args)) {
      params$allow_negative_rho_D <- as.logical(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--no_rebalance") {
      params$rebalance <- FALSE
      i <- i + 1
    } else if (args[i] == "--quiet") {
      params$verbose <- FALSE
      i <- i + 1
    } else if (args[i] == "--help") {
      cat("BATCH ACDE SIMULATION RUNNER\n")
      cat("============================\n")
      cat("Usage: Rscript batch_runner.R [options]\n\n")
      cat("Options:\n")
      cat("  --n_datasets N      Number of datasets to generate (default: 50)\n")
      cat("  --n_pairs N         Twin pairs per dataset (default: 400)\n")
      cat("  --engine_path PATH  Path to engine.R (default: engine.R)\n")
      cat("  --base_dir DIR      Base output directory (default: datasets)\n")
      cat("  --seed_start N      Starting seed for reproducibility (default: random)\n")
      cat("  --max_corr X        Maximum correlation magnitude (default: 0.8)\n")
      cat("  --stratified_prop X Proportion using stratified sampling (default: 0.6)\n")
      cat("  --n_pool N          Trait pool size (default: 20)\n")
      cat("  --run_id ID         Custom run identifier (default: timestamp)\n")
      cat("  --a_range MIN MAX   A component variance range (default: 0.01 0.64)\n")
      cat("  --c_range MIN MAX   C component variance range (default: 0.00 0.49)\n")
      cat("  --d_range MIN MAX   D component variance range (default: 0.00 0.36)\n")
      cat("  --e_range MIN MAX   E component variance range (default: 0.04 0.81)\n")
      cat("  --min_variance X    Minimum variance threshold (default: 0.01)\n")
      cat("  --allow_negative_rho_D BOOL  Allow negative dominance correlations (default: FALSE)\n")
      cat("  --no_rebalance      Disable correlation rebalancing\n")
      cat("  --quiet             Suppress progress output\n")
      cat("  --help              Show this help message\n\n")
      cat("Output Structure:\n")
      cat("  datasets/simulations_MMDDYYHHMM/\n")
      cat("    ├── simulated_data_1.csv\n")
      cat("    ├── simulated_data_2.csv\n")
      cat("    ├── ...\n")
      cat("    ├── simulated_data_N.csv\n")
      cat("    ├── params_log.csv\n")
      cat("    └── run_info.txt\n")
      quit(status = 0)
    } else {
      i <- i + 1
    }
  }

  return(params)
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

create_timestamp_id <- function() {
  # Create timestamp in MMDDYYHHMM format
  now <- Sys.time()
  formatted <- format(now, "%m%d%y%H%M")
  return(formatted)
}

create_output_directory <- function(base_dir, run_id = NULL) {
  # Create timestamped output directory

  if (is.null(run_id)) {
    timestamp_id <- create_timestamp_id()
    dir_name <- paste0("simulations_", timestamp_id)
  } else {
    dir_name <- paste0("simulations_", run_id)
  }

  output_dir <- file.path(base_dir, dir_name)

  # Create directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  } else {
    warning(sprintf("Directory %s already exists. Files may be overwritten.", output_dir))
  }

  return(list(
    dir = output_dir,
    timestamp = ifelse(is.null(run_id), timestamp_id, run_id)
  ))
}

source_engine <- function(engine_path) {
  # Source the engine script and verify required functions exist

  if (!file.exists(engine_path)) {
    stop(sprintf("Engine script not found: %s", engine_path))
  }

  tryCatch({
    source(engine_path)
  }, error = function(e) {
    stop(sprintf("Failed to source engine script %s: %s", engine_path, e$message))
  })

  # Verify required functions exist
  required_functions <- c("simulate_acde_dataset_enhanced", "get_truth", "get_params")
  missing_functions <- required_functions[!sapply(required_functions, exists)]

  if (length(missing_functions) > 0) {
    stop(sprintf("Missing required functions from engine: %s",
                 paste(missing_functions, collapse = ", ")))
  }

  return(TRUE)
}

# ==============================================================================
# BATCH SIMULATION FUNCTIONS
# ==============================================================================

generate_batch_datasets <- function(params, output_info, verbose = TRUE) {
  # Generate multiple datasets using the engine

  if (verbose) {
    cat("🚀 BATCH DATASET GENERATION\n")
    cat("===========================\n")
    cat(sprintf("Target datasets: %d\n", params$n_datasets))
    cat(sprintf("Twin pairs per dataset: %d\n", params$n_pairs))
    cat(sprintf("Output directory: %s\n", output_info$dir))
    cat(sprintf("Maximum correlation: ±%.2f\n", params$max_correlation))
    cat(sprintf("Stratified sampling: %.1f%%\n", params$stratified_prop * 100))
    cat(sprintf("Trait pool size: %d\n", params$n_pool))
    cat(sprintf("Rebalancing: %s\n", ifelse(params$rebalance, "enabled", "disabled")))
    if (!is.null(params$seed_start)) {
      cat(sprintf("Seed sequence: %d to %d\n", params$seed_start, params$seed_start + params$n_datasets - 1))
    }
    cat("\n")
  }

  # Initialize parameter log
  params_log <- data.frame()

  # Generate each dataset
  for (i in 1:params$n_datasets) {
    if (verbose && (i %% 10 == 0 || i == 1)) {
      cat(sprintf("Generating dataset %d/%d...\n", i, params$n_datasets))
    }

    # Set seed for reproducibility
    current_seed <- if (!is.null(params$seed_start)) params$seed_start + i - 1 else NULL
    if (!is.null(current_seed)) {
      set.seed(current_seed)
    }

    # Generate dataset
    tryCatch({
      dataset <- simulate_acde_dataset_enhanced(
        n_pairs = params$n_pairs,
        n_pool = params$n_pool,
        stratified_prop = params$stratified_prop,
        max_correlation = params$max_correlation,
        rebalance = params$rebalance,
        verbose = FALSE
      )

      # Save dataset
      data_filename <- sprintf("simulated_data_%d.csv", i)
      data_filepath <- file.path(output_info$dir, data_filename)
      write.csv(dataset, data_filepath, row.names = FALSE, quote = FALSE)

      # Extract truth values and metadata for parameter log
      truth <- get_truth(dataset)
      params_data <- get_params(dataset)

      # Create parameter log entry
      log_entry <- data.frame(
        replicate = i,
        seed = ifelse(is.null(current_seed), NA, current_seed),

        # Path coefficients (variance components)
        a11 = truth$a11,
        a22 = truth$a22,
        a12 = truth$a12,
        c11 = truth$c11,
        c22 = truth$c22,
        c12 = truth$c12,
        d11 = truth$d11,
        d22 = truth$d22,
        d12 = truth$d12,
        e11 = truth$e11,
        e22 = truth$e22,
        e12 = truth$e12,

        # Component correlations
        rho_A = params_data$correlations["rho_A"],
        rho_C = params_data$correlations["rho_C"],
        rho_D = params_data$correlations["rho_D"],
        rho_E = params_data$correlations["rho_E"],

        # Metadata
        implied_phenotypic_corr = ifelse(
          !is.null(params_data$meta$implied_phenotypic_corr),
          params_data$meta$implied_phenotypic_corr,
          NA
        ),
        sampling_strategy = ifelse(
          !is.null(params_data$meta$sampling_strategy),
          params_data$meta$sampling_strategy,
          "unknown"
        ),
        target_band = ifelse(
          !is.null(params_data$meta$target_band),
          params_data$meta$target_band,
          "none"
        ),
        retry_count = ifelse(
          !is.null(params_data$meta$retry_count),
          params_data$meta$retry_count,
          NA
        ),
        pool_size = ifelse(
          !is.null(params_data$meta$pool_size),
          params_data$meta$pool_size,
          params$n_pool
        ),

        stringsAsFactors = FALSE
      )

      # Add to parameter log
      params_log <- rbind(params_log, log_entry)

    }, error = function(e) {
      warning(sprintf("Failed to generate dataset %d: %s", i, e$message))

      # Create empty entry for failed generation
      empty_entry <- data.frame(
        replicate = i,
        seed = ifelse(is.null(current_seed), NA, current_seed),
        a11 = NA, a22 = NA, a12 = NA,
        c11 = NA, c22 = NA, c12 = NA,
        d11 = NA, d22 = NA, d12 = NA,
        e11 = NA, e22 = NA, e12 = NA,
        rho_A = NA, rho_C = NA, rho_D = NA, rho_E = NA,
        implied_phenotypic_corr = NA,
        sampling_strategy = "failed",
        target_band = "none",
        retry_count = NA,
        pool_size = NA,
        stringsAsFactors = FALSE
      )

      params_log <- rbind(params_log, empty_entry)
    })
  }

  return(params_log)
}

create_run_info <- function(params, output_info, params_log) {
  # Create comprehensive run information file

  info_file <- file.path(output_info$dir, "run_info.txt")

  # Calculate summary statistics
  successful_runs <- sum(!is.na(params_log$a11))
  failed_runs <- nrow(params_log) - successful_runs

  if (successful_runs > 0) {
    pheno_corrs <- params_log$implied_phenotypic_corr[!is.na(params_log$implied_phenotypic_corr)]
    corr_range <- if (length(pheno_corrs) > 0) range(pheno_corrs) else c(NA, NA)
    corr_mean <- if (length(pheno_corrs) > 0) mean(pheno_corrs) else NA
    corr_sd <- if (length(pheno_corrs) > 0) sd(pheno_corrs) else NA

    # Strategy counts
    strategy_counts <- table(params_log$sampling_strategy[!is.na(params_log$a11)])

    # Target band counts
    band_counts <- table(params_log$target_band[params_log$target_band != "none" &
                                               !is.na(params_log$a11)])
  }

  # Write run information
  cat("ACDE TWIN DATA SIMULATION RUN INFORMATION\n", file = info_file)
  cat("=========================================\n\n", file = info_file, append = TRUE)

  # Run parameters
  cat("RUN PARAMETERS:\n", file = info_file, append = TRUE)
  cat(sprintf("  Run ID: %s\n", output_info$timestamp), file = info_file, append = TRUE)
  cat(sprintf("  Generation Time: %s\n", Sys.time()), file = info_file, append = TRUE)
  cat(sprintf("  Engine Script: %s\n", params$engine_path), file = info_file, append = TRUE)
  cat(sprintf("  Target Datasets: %d\n", params$n_datasets), file = info_file, append = TRUE)
  cat(sprintf("  Twin Pairs per Dataset: %d\n", params$n_pairs), file = info_file, append = TRUE)
  cat(sprintf("  Maximum Correlation: ±%.2f\n", params$max_correlation), file = info_file, append = TRUE)
  cat(sprintf("  Stratified Sampling Proportion: %.2f\n", params$stratified_prop), file = info_file, append = TRUE)
  cat(sprintf("  Trait Pool Size: %d\n", params$n_pool), file = info_file, append = TRUE)
  cat(sprintf("  ACDE Variance Ranges:\n"), file = info_file, append = TRUE)
  cat(sprintf("    A: [%.4f, %.4f]\n", params$a_range[1], params$a_range[2]), file = info_file, append = TRUE)
  cat(sprintf("    C: [%.4f, %.4f]\n", params$c_range[1], params$c_range[2]), file = info_file, append = TRUE)
  cat(sprintf("    D: [%.4f, %.4f]\n", params$d_range[1], params$d_range[2]), file = info_file, append = TRUE)
  cat(sprintf("    E: [%.4f, %.4f]\n", params$e_range[1], params$e_range[2]), file = info_file, append = TRUE)
  cat(sprintf("  Minimum Variance Threshold: %.4f\n", params$min_variance), file = info_file, append = TRUE)
  cat(sprintf("  Allow Negative rho_D: %s\n", params$allow_negative_rho_D), file = info_file, append = TRUE)

  if (!is.null(params$seed_start)) {
    cat(sprintf("  Seed Range: %d to %d\n", params$seed_start, params$seed_start + params$n_datasets - 1), file = info_file, append = TRUE)
  }

  cat("\n", file = info_file, append = TRUE)

  # Generation results
  cat("GENERATION RESULTS:\n", file = info_file, append = TRUE)
  cat(sprintf("  Successful: %d/%d (%.1f%%)\n", successful_runs, params$n_datasets,
              100 * successful_runs / params$n_datasets), file = info_file, append = TRUE)
  if (failed_runs > 0) {
    cat(sprintf("  Failed: %d/%d (%.1f%%)\n", failed_runs, params$n_datasets,
                100 * failed_runs / params$n_datasets), file = info_file, append = TRUE)
  }

  if (successful_runs > 0) {
    cat("\n", file = info_file, append = TRUE)
    cat("CORRELATION ANALYSIS:\n", file = info_file, append = TRUE)
    cat(sprintf("  Phenotypic Correlation Range: [%.4f, %.4f]\n", corr_range[1], corr_range[2]), file = info_file, append = TRUE)
    cat(sprintf("  Mean Phenotypic Correlation: %.4f (±%.4f)\n", corr_mean, corr_sd), file = info_file, append = TRUE)

    cat("\n", file = info_file, append = TRUE)
    cat("SAMPLING STRATEGY DISTRIBUTION:\n", file = info_file, append = TRUE)
    for (strategy in names(strategy_counts)) {
      cat(sprintf("  %s: %d (%.1f%%)\n", strategy, strategy_counts[[strategy]],
                  100 * strategy_counts[[strategy]] / successful_runs), file = info_file, append = TRUE)
    }

    if (length(band_counts) > 0) {
      cat("\n", file = info_file, append = TRUE)
      cat("TARGET BAND DISTRIBUTION:\n", file = info_file, append = TRUE)
      for (band in names(band_counts)) {
        cat(sprintf("  %s: %d\n", band, band_counts[[band]]), file = info_file, append = TRUE)
      }
    }
  }

  cat("\n", file = info_file, append = TRUE)
  cat("OUTPUT FILES:\n", file = info_file, append = TRUE)
  cat("  Data files: simulated_data_1.csv through simulated_data_N.csv\n", file = info_file, append = TRUE)
  cat("  Parameter log: params_log.csv\n", file = info_file, append = TRUE)
  cat("  This file: run_info.txt\n", file = info_file, append = TRUE)

  return(invisible(info_file))
}

# ==============================================================================
# MAIN EXECUTION FUNCTION
# ==============================================================================

main_batch <- function() {
  # Parse arguments
  params <- parse_batch_args()

  if (params$verbose) {
    cat("📦 BATCH ACDE SIMULATION RUNNER\n")
    cat("================================\n\n")
  }

  # Check and source engine
  if (params$verbose) {
    cat(sprintf("🔧 Loading engine from: %s\n", params$engine_path))
  }
  source_engine(params$engine_path)

  # Create output directory
  output_info <- create_output_directory(params$base_dir, params$run_id)

  if (params$verbose) {
    cat(sprintf("📁 Created output directory: %s\n\n", output_info$dir))
  }

  # Generate datasets
  params_log <- generate_batch_datasets(params, output_info, params$verbose)

  # Save parameter log
  params_log_file <- file.path(output_info$dir, "params_log.csv")
  write.csv(params_log, params_log_file, row.names = FALSE, quote = FALSE)

  # Create run information
  create_run_info(params, output_info, params_log)

  if (params$verbose) {
    successful_runs <- sum(!is.na(params_log$a11))
    failed_runs <- nrow(params_log) - successful_runs

    cat("\n✅ BATCH GENERATION COMPLETE!\n")
    cat("==============================\n")
    cat(sprintf("Output directory: %s\n", output_info$dir))
    cat(sprintf("Successful datasets: %d/%d\n", successful_runs, params$n_datasets))
    if (failed_runs > 0) {
      cat(sprintf("Failed datasets: %d\n", failed_runs))
    }
    cat(sprintf("Parameter log: %s\n", params_log_file))

    if (successful_runs > 0) {
      pheno_corrs <- params_log$implied_phenotypic_corr[!is.na(params_log$implied_phenotypic_corr)]
      if (length(pheno_corrs) > 0) {
        cat(sprintf("Correlation range: [%.4f, %.4f]\n", min(pheno_corrs), max(pheno_corrs)))
      }
    }
  }

  return(invisible(list(
    output_dir = output_info$dir,
    params_log = params_log,
    successful_runs = sum(!is.na(params_log$a11)),
    total_runs = nrow(params_log)
  )))
}

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

# Function to quickly run a batch with common settings
quick_batch <- function(n_datasets = 50, n_pairs = 400, max_corr = 0.8,
                       seed_start = 12345, verbose = TRUE) {
  # Set up temporary parameters
  temp_args <- c(
    "--n_datasets", as.character(n_datasets),
    "--n_pairs", as.character(n_pairs),
    "--max_corr", as.character(max_corr),
    "--seed_start", as.character(seed_start)
  )

  if (!verbose) {
    temp_args <- c(temp_args, "--quiet")
  }

  # Override command line args temporarily
  old_args <- commandArgs(trailingOnly = TRUE)
  assignInNamespace("commandArgs", function(trailingOnly = FALSE) {
    if (trailingOnly) temp_args else c("Rscript", "batch_runner.R", temp_args)
  }, ns = "base")

  # Run batch
  result <- main_batch()

  # Restore original args
  assignInNamespace("commandArgs", function(trailingOnly = FALSE) {
    if (trailingOnly) old_args else c("Rscript", "batch_runner.R", old_args)
  }, ns = "base")

  return(result)
}

# ==============================================================================
# EXECUTE WHEN RUN AS SCRIPT
# ==============================================================================

if (!interactive()) {
  main_batch()
}

# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# COMMAND LINE USAGE:
#
# Basic usage (50 datasets, 400 pairs each):
# Rscript batch_runner.R
#
# Custom parameters:
# Rscript batch_runner.R --n_datasets 100 --n_pairs 500 --max_corr 0.7 --seed_start 12345
#
# High throughput run:
# Rscript batch_runner.R --n_datasets 200 --n_pairs 1000 --n_pool 30 --no_rebalance
#
# FUNCTION USAGE:
#
# Quick batch generation:
# result <- quick_batch(n_datasets = 100, n_pairs = 400, seed_start = 12345)
#
# Custom batch:
# source("batch_runner.R")
# result <- main_batch()  # Uses command line args or interactive prompts
