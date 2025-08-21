#!/usr/bin/env Rscript

# ==============================================================================
# BINNED TWIN DATA SIMULATION - GENERATE GOOD DATA THEN SORT BY EFFECT SIZES
# ==============================================================================
# Generates realistic ACE data, then bins by actual component magnitudes

library(MASS)

# ==============================================================================
# COMMAND LINE ARGUMENT PARSING
# ==============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  params <- list(
    n_pairs = 400,
    n_reps_per_bin = 100,
    outdir = "datasets",
    n_attempts = 5000  # Generate this many datasets to fill bins
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--n_pairs" && i < length(args)) {
      params$n_pairs <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--n_reps_per_bin" && i < length(args)) {
      params$n_reps_per_bin <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--outdir" && i < length(args)) {
      params$outdir <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--n_attempts" && i < length(args)) {
      params$n_attempts <- as.integer(args[i + 1])
      i <- i + 2
    } else {
      i <- i + 1
    }
  }

  return(params)
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Generate realistic ACE parameters with broader ranges
generate_realistic_ace_params <- function() {
  # Use broader uniform distributions to ensure we hit all effect sizes
  # Generate proportions that will naturally span small, medium, large effects

  # Generate raw proportions with wider ranges
  A_t1 <- runif(1, 0.10, 0.80)  # Allow wide range for A
  C_t1 <- runif(1, 0.00, 0.60)  # Allow wide range for C
  E_t1 <- runif(1, 0.15, 0.85)  # Allow wide range for E

  # Normalize trait 1
  total_t1 <- A_t1 + C_t1 + E_t1
  A_t1 <- A_t1 / total_t1
  C_t1 <- C_t1 / total_t1
  E_t1 <- E_t1 / total_t1

  # Generate trait 2 with some correlation to trait 1 but allow variation
  correlation_strength <- 0.3  # Moderate correlation between traits

  A_t2 <- A_t1 + runif(1, -0.3, 0.3) * (1 - correlation_strength)
  C_t2 <- C_t1 + runif(1, -0.3, 0.3) * (1 - correlation_strength)
  E_t2 <- E_t1 + runif(1, -0.3, 0.3) * (1 - correlation_strength)

  # Ensure bounds
  A_t2 <- pmax(0.05, pmin(0.85, A_t2))
  C_t2 <- pmax(0.00, pmin(0.65, C_t2))
  E_t2 <- pmax(0.10, pmin(0.90, E_t2))

  # Normalize trait 2
  total_t2 <- A_t2 + C_t2 + E_t2
  A_t2 <- A_t2 / total_t2
  C_t2 <- C_t2 / total_t2
  E_t2 <- E_t2 / total_t2

  # Generate cross-trait correlations
  rho_A <- runif(1, -0.30, 0.30)
  rho_C <- runif(1, -0.20, 0.20)
  rho_E <- runif(1, -0.15, 0.15)

  return(list(
    A_diag = c(A_t1, A_t2),
    C_diag = c(C_t1, C_t2),
    E_diag = c(E_t1, E_t2),
    rho_A = rho_A,
    rho_C = rho_C,
    rho_E = rho_E
  ))
}

# Build component matrix safely
build_component_matrix <- function(diag_vals, rho) {
  # Ensure correlation is safe for PD
  max_safe_rho <- 0.95 / sqrt(max(diag_vals) / min(diag_vals))
  rho_safe <- sign(rho) * min(abs(rho), max_safe_rho)

  off_diag <- rho_safe * sqrt(diag_vals[1] * diag_vals[2])

  matrix(c(diag_vals[1], off_diag, off_diag, diag_vals[2]), nrow = 2, ncol = 2)
}

# Generate twin data
simulate_twin_group <- function(n_pairs, A, C, E, r_A, r_C, group_info) {
  # Build covariance matrix
  Sigma_within <- A + C + E
  Sigma_between <- r_A * A + r_C * C

  Sigma <- rbind(
    cbind(Sigma_within, Sigma_between),
    cbind(Sigma_between, Sigma_within)
  )

  # Check PD and add small regularization if needed
  eigenvals <- eigen(Sigma, only.values = TRUE)$values
  if (any(eigenvals <= 1e-8)) {
    Sigma <- Sigma + diag(1e-6, 4)
  }

  # Generate data
  data <- mvrnorm(n = n_pairs, mu = rep(0, 4), Sigma = Sigma)

  # Format with correct column order for model
  df <- data.frame(
    zyg = group_info$zyg,
    zyg2 = group_info$zyg2,
    sex1 = group_info$sex1,
    sex2 = group_info$sex2,
    p1_t1 = data[, 1],  # Twin 1, Trait 1
    p2_t1 = data[, 3],  # Twin 2, Trait 1
    p1_t2 = data[, 2],  # Twin 1, Trait 2
    p2_t2 = data[, 4]   # Twin 2, Trait 2
  )

  return(df)
}

# ==============================================================================
# EFFECT SIZE CLASSIFICATION
# ==============================================================================

classify_effect_size <- function(value, thresholds = c(sm = 0.25, lg = 0.45)) {
  if (value < thresholds["sm"]) return("sm")
  else if (value < thresholds["lg"]) return("md")
  else return("lg")
}

determine_bin <- function(params) {
  # Simple classification: just look at individual components
  # This is much more permissive than the previous complex logic

  A_t1_class <- classify_effect_size(params$A_diag[1])
  A_t2_class <- classify_effect_size(params$A_diag[2])
  C_t1_class <- classify_effect_size(params$C_diag[1])
  C_t2_class <- classify_effect_size(params$C_diag[2])
  E_t1_class <- classify_effect_size(params$E_diag[1])
  E_t2_class <- classify_effect_size(params$E_diag[2])

  # Return all possible bins this dataset could fill
  bins <- c(
    paste0("a11_", A_t1_class),
    paste0("a22_", A_t2_class),
    paste0("c11_", C_t1_class),
    paste0("c22_", C_t2_class),
    paste0("e11_", E_t1_class),
    paste0("e22_", E_t2_class)
  )

  return(bins)
}

# ==============================================================================
# MAIN SIMULATION FUNCTION
# ==============================================================================

main <- function() {
  params <- parse_args()

  cat("🎯 BINNED TWIN DATA SIMULATION\n")
  cat("===============================\n")
  cat(sprintf("Target: %d datasets per bin\n", params$n_reps_per_bin))
  cat(sprintf("Sample size: %d pairs per group\n", params$n_pairs))
  cat(sprintf("Max attempts: %d\n", params$n_attempts))
  cat(sprintf("Output directory: %s\n\n", params$outdir))

  # Define target bins
  target_bins <- c(
    "a11_sm", "a11_md", "a11_lg",
    "a22_sm", "a22_md", "a22_lg",
    "c11_sm", "c11_md", "c11_lg",
    "c22_sm", "c22_md", "c22_lg",
    "e11_sm", "e11_md", "e11_lg",
    "e22_sm", "e22_md", "e22_lg"
  )

  # Initialize bin storage
  bin_storage <- list()
  bin_counts <- list()
  for (bin in target_bins) {
    bin_storage[[bin]] <- list()
    bin_counts[[bin]] <- 0
  }

  # Twin group definitions
  groups <- list(
    MZf = list(zyg = 1, zyg2 = 1, sex1 = 2, sex2 = 2, r_A = 1.0, r_C = 1.0),
    DZf = list(zyg = 2, zyg2 = 2, sex1 = 2, sex2 = 2, r_A = 0.5, r_C = 1.0),
    MZm = list(zyg = 3, zyg2 = 1, sex1 = 1, sex2 = 1, r_A = 1.0, r_C = 1.0),
    DZm = list(zyg = 4, zyg2 = 2, sex1 = 1, sex2 = 1, r_A = 0.5, r_C = 1.0),
    DZo = list(zyg = 5, zyg2 = 2, sex1 = 1, sex2 = 2, r_A = 0.5, r_C = 1.0)
  )

  # Generate datasets and bin them
  attempt <- 1
  while (attempt <= params$n_attempts && any(sapply(bin_counts, function(x) x < params$n_reps_per_bin))) {

    # Progress update
    if (attempt %% 100 == 0) {
      filled_bins <- sum(sapply(bin_counts, function(x) x >= params$n_reps_per_bin))
      cat(sprintf("Attempt %d: %d/%d bins filled\n", attempt, filled_bins, length(target_bins)))
    }

    # Generate realistic parameters
    ace_params <- generate_realistic_ace_params()

    # Build component matrices
    A <- build_component_matrix(ace_params$A_diag, ace_params$rho_A)
    C <- build_component_matrix(ace_params$C_diag, ace_params$rho_C)
    E <- build_component_matrix(ace_params$E_diag, ace_params$rho_E)

    # Check which bins this dataset could fill
    applicable_bins <- determine_bin(ace_params)
    available_bins <- intersect(applicable_bins, target_bins)
    available_bins <- available_bins[sapply(available_bins, function(b) bin_counts[[b]] < params$n_reps_per_bin)]

    if (length(available_bins) > 0) {
      # Generate the twin data
      tryCatch({
        all_data <- list()
        for (group_name in names(groups)) {
          group_info <- groups[[group_name]]
          group_data <- simulate_twin_group(params$n_pairs, A, C, E,
                                           group_info$r_A, group_info$r_C, group_info)
          all_data[[group_name]] <- group_data
        }

        df <- do.call(rbind, all_data)

        # Store in the first available bin (could be modified to store in multiple)
        chosen_bin <- available_bins[1]
        bin_counts[[chosen_bin]] <- bin_counts[[chosen_bin]] + 1

        # Store dataset info
        dataset_info <- list(
          data = df,
          params = ace_params,
          bin = chosen_bin,
          replicate = bin_counts[[chosen_bin]]
        )

        bin_storage[[chosen_bin]][[bin_counts[[chosen_bin]]]] <- dataset_info

      }, error = function(e) {
        # Skip this attempt if data generation fails
      })
    }

    attempt <- attempt + 1
  }

  # Final status
  cat("\n📊 FINAL BIN STATUS:\n")
  cat("====================\n")
  for (bin in target_bins) {
    status <- if (bin_counts[[bin]] >= params$n_reps_per_bin) "✅ COMPLETE" else sprintf("❌ %d/%d", bin_counts[[bin]], params$n_reps_per_bin)
    cat(sprintf("%-8s: %s\n", bin, status))
  }

  # Save completed bins
  cat("\n💾 SAVING DATASETS...\n")
  saved_bins <- 0
  for (bin_name in names(bin_storage)) {
    if (bin_counts[[bin_name]] >= params$n_reps_per_bin) {

      # Create directory
      bin_dir <- file.path(params$outdir, paste0("data", bin_name))
      dir.create(bin_dir, recursive = TRUE, showWarnings = FALSE)

      # Save datasets and parameters
      params_log <- data.frame()

      for (i in 1:params$n_reps_per_bin) {
        dataset_info <- bin_storage[[bin_name]][[i]]

        # Save data file
        data_file <- file.path(bin_dir, sprintf("simulated_data_%d.csv", i))
        write.csv(dataset_info$data, data_file, row.names = FALSE, quote = FALSE)

        # Prepare parameter log entry
        param_row <- data.frame(
          comp = strsplit(bin_name, "_")[[1]][1],
          trait = paste0(strsplit(bin_name, "_")[[1]][2:3], collapse = "_"),
          eff_level = strsplit(bin_name, "_")[[1]][length(strsplit(bin_name, "_")[[1]])],
          rep = i,
          A_diag_t1 = dataset_info$params$A_diag[1],
          C_diag_t1 = dataset_info$params$C_diag[1],
          E_diag_t1 = dataset_info$params$E_diag[1],
          A_diag_t2 = dataset_info$params$A_diag[2],
          C_diag_t2 = dataset_info$params$C_diag[2],
          E_diag_t2 = dataset_info$params$E_diag[2],
          rho_A = dataset_info$params$rho_A,
          rho_C = dataset_info$params$rho_C,
          rho_E = dataset_info$params$rho_E
        )

        params_log <- rbind(params_log, param_row)
      }

      # Save parameter log
      params_file <- file.path(bin_dir, "params_log.csv")
      write.csv(params_log, params_file, row.names = FALSE, quote = FALSE)

      saved_bins <- saved_bins + 1
      cat(sprintf("✅ Saved bin: %s (%d datasets)\n", bin_name, params$n_reps_per_bin))
    }
  }

  cat(sprintf("\n🎉 SIMULATION COMPLETE! Saved %d/%d bins\n", saved_bins, length(target_bins)))
  cat(sprintf("Results saved to: %s\n", params$outdir))

  # Return summary for inspection
  return(invisible(list(
    bin_counts = bin_counts,
    target_bins = target_bins,
    saved_bins = saved_bins
  )))
}

# ==============================================================================
# EXECUTE
# ==============================================================================

if (!interactive()) {
  main()
}
