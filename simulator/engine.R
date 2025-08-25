#!/usr/bin/env Rscript

# ==============================================================================
# ENHANCED ACDE TWIN DATA SIMULATION ENGINE
# ==============================================================================
# Generates realistic twin datasets with stratified sampling across the full
# cross-trait correlation space. Supports both random and structured coverage
# with automatic rebalancing for under-sampled correlation regions.

library(MASS)
library(Matrix)

# ==============================================================================
# COMMAND LINE ARGUMENT PARSING
# ==============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  params <- list(
    n_pairs = 400,              # Total pairs (distributed across 5 groups)
    outfile = NULL,             # Output CSV file (if NULL, returns data frame)
    seed = NULL,                # Random seed for reproducibility
    verbose = TRUE,             # Print progress information
    n_pool = 20,                # Size of trait pool for pairing
    stratified_prop = 0.6,      # Proportion using stratified sampling
    target_bands = list(        # Target correlation bands for stratification
      high = c(0.5, 0.8),
      moderate = c(0.2, 0.5),
      low = c(0.0, 0.2)
    ),
    rebalance = TRUE,           # Enable rebalancing for under-sampled regions
    max_correlation = 0.8       # Maximum absolute correlation allowed
  )

  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--n_pairs" && i < length(args)) {
      params$n_pairs <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--outfile" && i < length(args)) {
      params$outfile <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--seed" && i < length(args)) {
      params$seed <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--n_pool" && i < length(args)) {
      params$n_pool <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--max_corr" && i < length(args)) {
      params$max_correlation <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--no_stratify") {
      params$stratified_prop <- 0.0
      i <- i + 1
    } else if (args[i] == "--no_rebalance") {
      params$rebalance <- FALSE
      i <- i + 1
    } else if (args[i] == "--quiet") {
      params$verbose <- FALSE
      i <- i + 1
    } else if (args[i] == "--help") {
      cat("Usage: Rscript engine.R [options]\n")
      cat("Options:\n")
      cat("  --n_pairs N       Total number of twin pairs (default: 400)\n")
      cat("  --outfile F       Output CSV file (default: return data frame)\n")
      cat("  --seed S          Random seed (default: random)\n")
      cat("  --n_pool N        Size of trait pool (default: 20)\n")
      cat("  --max_corr X      Maximum correlation magnitude (default: 0.8)\n")
      cat("  --no_stratify     Disable stratified sampling\n")
      cat("  --no_rebalance    Disable correlation rebalancing\n")
      cat("  --quiet           Suppress progress output\n")
      cat("  --help            Show this help message\n")
      quit(status = 0)
    } else {
      i <- i + 1
    }
  }

  return(params)
}

# ==============================================================================
# PARAMETER GENERATION FUNCTIONS
# ==============================================================================

generate_single_trait_params <- function(a_range = c(0.01, 0.64),
                                         c_range = c(0.00, 0.49),
                                         d_range = c(0.00, 0.36),
                                         e_range = c(0.04, 0.81),
                                         min_variance = 0.01) {
  # Generate path coefficients for a single trait with specified variance ranges
  # Note: Input ranges are for VARIANCES, we generate PATHS and then square them

  # Convert variance ranges to path ranges (square root)
  a_path_range <- sqrt(a_range)
  c_path_range <- sqrt(c_range)
  d_path_range <- sqrt(d_range)
  e_path_range <- sqrt(e_range)

  # Generate paths from specified ranges
  a <- runif(1, a_path_range[1], a_path_range[2])
  c <- runif(1, c_path_range[1], c_path_range[2])
  d <- runif(1, d_path_range[1], d_path_range[2])
  e <- runif(1, e_path_range[1], e_path_range[2])

  # Normalize to ensure total variance = 1
  total_var <- a^2 + c^2 + d^2 + e^2
  scale_factor <- sqrt(1 / total_var)

  a_scaled <- a * scale_factor
  c_scaled <- c * scale_factor
  d_scaled <- d * scale_factor
  e_scaled <- e * scale_factor

  # Check minimum variance constraint
  variances <- c(a_scaled^2, c_scaled^2, d_scaled^2, e_scaled^2)
  if (any(variances < min_variance)) {
    return(NULL)  # Signal retry needed
  }

  return(list(
    a = a_scaled,
    c = c_scaled,
    d = d_scaled,
    e = e_scaled
  ))
}

simulate_trait_pool <- function(n_pool = 20, a_range = c(0.01, 0.64),
                               c_range = c(0.00, 0.49), d_range = c(0.00, 0.36),
                               e_range = c(0.04, 0.81), min_variance = 0.01,
                               max_attempts = 1000) {
  # Generate a pool of independent trait structures with constraints

  pool <- list()
  attempts <- 0

  while (length(pool) < n_pool && attempts < max_attempts) {
    attempts <- attempts + 1

    trait_candidate <- generate_single_trait_params(
      a_range = a_range, c_range = c_range,
      d_range = d_range, e_range = e_range,
      min_variance = min_variance
    )

    if (!is.null(trait_candidate)) {
      pool[[length(pool) + 1]] <- trait_candidate
    }
  }

  if (length(pool) < n_pool) {
    warning(sprintf("Could only generate %d/%d valid traits after %d attempts",
                    length(pool), n_pool, max_attempts))
  }

  return(pool)
}

calculate_implied_correlations <- function(trait1, trait2, rho_A, rho_C, rho_D, rho_E) {
  # Calculate what the actual cross-trait correlations would be
  # given the component correlations and trait structures

  # Phenotypic correlation is sum of component contributions
  rho_pheno <- (trait1$a * trait2$a * rho_A +
                trait1$c * trait2$c * rho_C +
                trait1$d * trait2$d * rho_D +
                trait1$e * trait2$e * rho_E)

  return(rho_pheno)
}

generate_stratified_correlations <- function(target_band, max_correlation = 0.8,
                                             allow_negative_rho_D = FALSE) {
  # Generate correlations targeting specific phenotypic correlation bands

  # Expand ranges based on scientific plausibility
  # Genetic correlations can be quite high
  rho_A_range <- c(-max_correlation, max_correlation)

  # Environmental correlations typically more moderate
  rho_C_range <- c(-min(0.5, max_correlation), min(0.5, max_correlation))
  rho_E_range <- c(-min(0.4, max_correlation), min(0.4, max_correlation))

  # Dominance correlation range - constrain based on allow_negative_rho_D
  if (allow_negative_rho_D) {
    rho_D_range <- c(-min(0.5, max_correlation), min(0.5, max_correlation))
  } else {
    rho_D_range <- c(0.0, min(0.5, max_correlation))
  }

  # Sample correlations with bias toward target band
  if (!is.null(target_band)) {
    # For structured sampling, bias toward combinations likely to yield target range
    target_magnitude <- mean(target_band)

    # Bias genetic correlation toward target
    rho_A <- runif(1,
                   max(rho_A_range[1], -target_magnitude - 0.2),
                   min(rho_A_range[2], target_magnitude + 0.2))

    # Other correlations sampled more broadly but with some constraint
    rho_C <- runif(1, rho_C_range[1], rho_C_range[2])
    rho_D <- runif(1, rho_D_range[1], rho_D_range[2])
    rho_E <- runif(1, rho_E_range[1], rho_E_range[2])

  } else {
    # Pure random sampling across full ranges
    rho_A <- runif(1, rho_A_range[1], rho_A_range[2])
    rho_C <- runif(1, rho_C_range[1], rho_C_range[2])
    rho_D <- runif(1, rho_D_range[1], rho_D_range[2])
    rho_E <- runif(1, rho_E_range[1], rho_E_range[2])
  }

  return(list(
    rho_A = rho_A,
    rho_C = rho_C,
    rho_D = rho_D,
    rho_E = rho_E
  ))
}

find_correlated_trait_stratified <- function(trait1, pool, target_band = NULL,
                                           max_attempts = 200, max_correlation = 0.8,
                                           allow_negative_rho_D = FALSE) {
  # Enhanced trait pairing with stratified correlation sampling and constraints

  for (attempt in 1:max_attempts) {
    # Randomly select a trait from pool
    trait2_candidate <- pool[[sample(length(pool), 1)]]

    # Generate correlations (stratified or random) with constraints
    correlations <- generate_stratified_correlations(target_band, max_correlation,
                                                   allow_negative_rho_D)

    # Calculate implied phenotypic correlation
    implied_corr <- calculate_implied_correlations(
      trait1, trait2_candidate,
      correlations$rho_A, correlations$rho_C,
      correlations$rho_D, correlations$rho_E
    )

    # Check if correlation is in acceptable range
    acceptable <- if (!is.null(target_band)) {
      abs(implied_corr) >= target_band[1] && abs(implied_corr) <= target_band[2]
    } else {
      abs(implied_corr) <= max_correlation  # Any correlation up to max
    }

    if (acceptable) {
      return(list(
        trait2 = trait2_candidate,
        rho_A = correlations$rho_A,
        rho_C = correlations$rho_C,
        rho_D = correlations$rho_D,
        rho_E = correlations$rho_E,
        implied_phenotypic_corr = implied_corr,
        target_band = target_band
      ))
    }
  }

  # Return best attempt even if not perfect
  return(list(
    trait2 = trait2_candidate,
    rho_A = correlations$rho_A,
    rho_C = correlations$rho_C,
    rho_D = correlations$rho_D,
    rho_E = correlations$rho_E,
    implied_phenotypic_corr = implied_corr,
    target_band = target_band,
    warning = "Target band not achieved"
  ))
}

select_sampling_strategy <- function(iteration, total_iterations, stratified_prop = 0.6,
                                   target_bands = list(high = c(0.5, 0.8),
                                                      moderate = c(0.2, 0.5),
                                                      low = c(0.0, 0.2))) {
  # Determine sampling strategy for current iteration

  if (runif(1) < stratified_prop) {
    # Stratified sampling - cycle through target bands
    band_names <- names(target_bands)
    band_index <- ((iteration - 1) %% length(band_names)) + 1
    return(list(
      strategy = "stratified",
      target_band = target_bands[[band_index]],
      band_name = band_names[band_index]
    ))
  } else {
    # Random sampling
    return(list(
      strategy = "random",
      target_band = NULL,
      band_name = "random"
    ))
  }
}

generate_acde_params_enhanced <- function(iteration = 1, total_iterations = 1,
                                        n_pool = 20, stratified_prop = 0.6,
                                        target_bands = list(high = c(0.5, 0.8),
                                                           moderate = c(0.2, 0.5),
                                                           low = c(0.0, 0.2)),
                                        max_correlation = 0.8, max_retries = 15,

                                        # NEW: ACDE component constraints
                                        a_range = c(0.01, 0.64),
                                        c_range = c(0.00, 0.49),
                                        d_range = c(0.00, 0.36),
                                        e_range = c(0.04, 0.81),
                                        min_variance = 0.01,
                                        allow_negative_rho_D = FALSE) {
  # Enhanced parameter generation with mathematical constraints

  # Select sampling strategy
  strategy <- select_sampling_strategy(iteration, total_iterations,
                                     stratified_prop, target_bands)

  for (retry in 1:max_retries) {
    # Step 1: Generate trait 1 with constraints
    trait1 <- NULL
    trait1_attempts <- 0
    max_trait_attempts <- 100

    while (is.null(trait1) && trait1_attempts < max_trait_attempts) {
      trait1_attempts <- trait1_attempts + 1
      trait1 <- generate_single_trait_params(
        a_range = a_range, c_range = c_range,
        d_range = d_range, e_range = e_range,
        min_variance = min_variance
      )
    }

    if (is.null(trait1)) {
      next  # Skip this retry if can't generate valid trait1
    }

    # Step 2: Generate pool of candidate trait 2 structures with constraints
    trait_pool <- simulate_trait_pool(
      n_pool = n_pool,
      a_range = a_range, c_range = c_range,
      d_range = d_range, e_range = e_range,
      min_variance = min_variance
    )

    if (length(trait_pool) == 0) {
      next  # Skip this retry if no valid traits in pool
    }

    # Step 3: Find compatible trait 2 from pool with correlation constraints
    match_result <- find_correlated_trait_stratified(
      trait1, trait_pool, strategy$target_band,
      max_correlation = max_correlation,
      allow_negative_rho_D = allow_negative_rho_D
    )

    if (!is.null(match_result) && is.null(match_result$warning)) {
      # Success! Construct final parameter set
      trait2 <- match_result$trait2

      return(list(
        paths = c(a1 = trait1$a, c1 = trait1$c, d1 = trait1$d, e1 = trait1$e,
                  a2 = trait2$a, c2 = trait2$c, d2 = trait2$d, e2 = trait2$e),
        correlations = c(rho_A = match_result$rho_A,
                        rho_C = match_result$rho_C,
                        rho_D = match_result$rho_D,
                        rho_E = match_result$rho_E),
        meta = list(
          implied_phenotypic_corr = match_result$implied_phenotypic_corr,
          retry_count = retry,
          pool_size = length(trait_pool),
          sampling_strategy = strategy$strategy,
          target_band = strategy$band_name,
          iteration = iteration,
          constraint_info = list(
            a_range = a_range, c_range = c_range,
            d_range = d_range, e_range = e_range,
            min_variance = min_variance,
            allow_negative_rho_D = allow_negative_rho_D
          )
        )
      ))
    }
  }

  # Use the last attempt even if not perfect
  if (!is.null(match_result)) {
    trait2 <- match_result$trait2

    return(list(
      paths = c(a1 = trait1$a, c1 = trait1$c, d1 = trait1$d, e1 = trait1$e,
                a2 = trait2$a, c2 = trait2$c, d2 = trait2$d, e2 = trait2$e),
      correlations = c(rho_A = match_result$rho_A,
                      rho_C = match_result$rho_C,
                      rho_D = match_result$rho_D,
                      rho_E = match_result$rho_E),
      meta = list(
        implied_phenotypic_corr = match_result$implied_phenotypic_corr,
        retry_count = max_retries,
        pool_size = length(trait_pool),
        sampling_strategy = strategy$strategy,
        target_band = strategy$band_name,
        iteration = iteration,
        warning = "Used imperfect match",
        constraint_info = list(
          a_range = a_range, c_range = c_range,
          d_range = d_range, e_range = e_range,
          min_variance = min_variance,
          allow_negative_rho_D = allow_negative_rho_D
        )
      )
    ))
  }

  # Final fallback with constraints
  return(generate_acde_params_fallback(
    a_range = a_range, c_range = c_range,
    d_range = d_range, e_range = e_range,
    min_variance = min_variance,
    allow_negative_rho_D = allow_negative_rho_D
  ))
}

generate_acde_params_fallback <- function(a_range = c(0.01, 0.64),
                                         c_range = c(0.00, 0.49),
                                         d_range = c(0.00, 0.36),
                                         e_range = c(0.04, 0.81),
                                         min_variance = 0.01,
                                         allow_negative_rho_D = FALSE) {
  # Simplified fallback method with constraints

  # Generate traits with constraints
  trait1 <- NULL
  trait2 <- NULL
  attempts <- 0
  max_attempts <- 200

  while ((is.null(trait1) || is.null(trait2)) && attempts < max_attempts) {
    attempts <- attempts + 1

    if (is.null(trait1)) {
      trait1 <- generate_single_trait_params(
        a_range = a_range, c_range = c_range,
        d_range = d_range, e_range = e_range,
        min_variance = min_variance
      )
    }

    if (is.null(trait2)) {
      trait2 <- generate_single_trait_params(
        a_range = a_range, c_range = c_range,
        d_range = d_range, e_range = e_range,
        min_variance = min_variance
      )
    }
  }

  # If still null, use very basic generation
  if (is.null(trait1)) {
    trait1 <- list(a = 0.4, c = 0.3, d = 0.2, e = 0.5)
    total1 <- sqrt(trait1$a^2 + trait1$c^2 + trait1$d^2 + trait1$e^2)
    trait1 <- lapply(trait1, function(x) x / total1)
  }

  if (is.null(trait2)) {
    trait2 <- list(a = 0.5, c = 0.2, d = 0.3, e = 0.4)
    total2 <- sqrt(trait2$a^2 + trait2$c^2 + trait2$d^2 + trait2$e^2)
    trait2 <- lapply(trait2, function(x) x / total2)
  }

  # Generate correlations with constraints
  rho_A <- runif(1, -0.5, 0.5)
  rho_C <- runif(1, -0.3, 0.3)
  rho_E <- runif(1, -0.25, 0.25)

  # Handle dominance correlation constraint
  if (allow_negative_rho_D) {
    rho_D <- runif(1, -0.3, 0.3)
  } else {
    rho_D <- runif(1, 0.0, 0.3)
  }

  return(list(
    paths = c(a1 = trait1$a, c1 = trait1$c, d1 = trait1$d, e1 = trait1$e,
              a2 = trait2$a, c2 = trait2$c, d2 = trait2$d, e2 = trait2$e),
    correlations = c(rho_A = rho_A, rho_C = rho_C, rho_D = rho_D, rho_E = rho_E),
    meta = list(
      implied_phenotypic_corr = NA,
      sampling_strategy = "fallback",
      constraint_info = list(
        a_range = a_range, c_range = c_range,
        d_range = d_range, e_range = e_range,
        min_variance = min_variance,
        allow_negative_rho_D = allow_negative_rho_D
      )
    )
  ))
}

# ==============================================================================
# CORRELATION SPACE ANALYSIS FUNCTIONS
# ==============================================================================

analyze_correlation_distribution <- function(phenotypic_corrs, target_bands = list(
                                            high = c(0.5, 0.8),
                                            moderate = c(0.2, 0.5),
                                            low = c(0.0, 0.2))) {
  # Analyze empirical distribution of phenotypic correlations

  abs_corrs <- abs(phenotypic_corrs)

  # Count samples in each band
  band_counts <- list()
  for (band_name in names(target_bands)) {
    band <- target_bands[[band_name]]
    in_band <- abs_corrs >= band[1] & abs_corrs <= band[2]
    band_counts[[band_name]] <- sum(in_band, na.rm = TRUE)
  }

  # Calculate total and proportions
  total_samples <- length(phenotypic_corrs[!is.na(phenotypic_corrs)])
  band_proportions <- lapply(band_counts, function(x) x / total_samples)

  return(list(
    total_samples = total_samples,
    band_counts = band_counts,
    band_proportions = band_proportions,
    empirical_range = range(phenotypic_corrs, na.rm = TRUE),
    empirical_mean = mean(phenotypic_corrs, na.rm = TRUE),
    empirical_sd = sd(phenotypic_corrs, na.rm = TRUE)
  ))
}

identify_undersample_regions <- function(correlation_analysis, min_proportion = 0.15) {
  # Identify correlation bands that are under-sampled

  undersample_bands <- c()

  for (band_name in names(correlation_analysis$band_proportions)) {
    if (correlation_analysis$band_proportions[[band_name]] < min_proportion) {
      undersample_bands <- c(undersample_bands, band_name)
    }
  }

  return(undersample_bands)
}

# ==============================================================================
# MATRIX CONSTRUCTION FUNCTIONS
# ==============================================================================

build_component_matrix <- function(path1, path2, rho) {
  # Build 2x2 covariance matrix from paths and correlation

  # Variance components (squared paths)
  var1 <- path1^2
  var2 <- path2^2

  # Ensure correlation is feasible for positive definiteness
  max_corr <- 0.95  # Conservative upper bound
  rho_safe <- sign(rho) * min(abs(rho), max_corr)

  # Covariance (correlation * geometric mean of standard deviations)
  cov12 <- rho_safe * path1 * path2

  # Build matrix
  matrix(c(var1, cov12, cov12, var2), nrow = 2, ncol = 2)
}

ensure_positive_definite <- function(mat, min_eigenval = 1e-6) {
  # Ensure matrix is positive definite using eigenvalue flooring

  eig <- eigen(mat, symmetric = TRUE)

  # Floor eigenvalues
  eig$values[eig$values < min_eigenval] <- min_eigenval

  # Reconstruct matrix
  reconstructed <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

  # Ensure symmetry (numerical precision)
  (reconstructed + t(reconstructed)) / 2
}

build_acde_covariances <- function(params) {
  # Extract parameters
  paths <- params$paths
  cors <- params$correlations

  # Build component matrices
  A <- build_component_matrix(paths["a1"], paths["a2"], cors["rho_A"])
  C <- build_component_matrix(paths["c1"], paths["c2"], cors["rho_C"])
  D <- build_component_matrix(paths["d1"], paths["d2"], cors["rho_D"])
  E <- build_component_matrix(paths["e1"], paths["e2"], cors["rho_E"])

  # Ensure all matrices are positive definite
  A <- ensure_positive_definite(A)
  C <- ensure_positive_definite(C)
  D <- ensure_positive_definite(D)
  E <- ensure_positive_definite(E)

  # Build full 4x4 twin covariance matrices

  # MZ twins: full sharing for A, C, D
  Sigma_MZ <- rbind(
    cbind(A + C + D + E, A + C + D),      # Twin1 with Twin1, Twin1 with Twin2
    cbind(A + C + D, A + C + D + E)       # Twin2 with Twin1, Twin2 with Twin2
  )

  # DZ twins: 0.5 for A, 1.0 for C, 0.25 for D
  Sigma_DZ <- rbind(
    cbind(A + C + D + E, 0.5 * A + C + 0.25 * D),
    cbind(0.5 * A + C + 0.25 * D, A + C + D + E)
  )

  # Ensure positive definiteness of full matrices
  Sigma_MZ <- ensure_positive_definite(Sigma_MZ)
  Sigma_DZ <- ensure_positive_definite(Sigma_DZ)

  return(list(
    A = A, C = C, D = D, E = E,
    Sigma_MZ = Sigma_MZ,
    Sigma_DZ = Sigma_DZ
  ))
}

# ==============================================================================
# DATA SIMULATION FUNCTIONS
# ==============================================================================

simulate_twin_group <- function(n_pairs, Sigma, group_info) {
  # Generate multivariate normal data for twin pairs

  # Generate data: 4 columns (trait1_twin1, trait2_twin1, trait1_twin2, trait2_twin2)
  raw_data <- mvrnorm(n = n_pairs, mu = rep(0, 4), Sigma = Sigma)

  # Format with consistent column names and metadata
  df <- data.frame(
    zyg = rep(group_info$zyg, n_pairs),
    zyg2 = rep(group_info$zyg2, n_pairs),
    sex1 = rep(group_info$sex1, n_pairs),
    sex2 = rep(group_info$sex2, n_pairs),
    p1_t1 = raw_data[, 1],  # Twin 1, Trait 1
    p1_t2 = raw_data[, 2],  # Twin 1, Trait 2
    p2_t1 = raw_data[, 3],  # Twin 2, Trait 1
    p2_t2 = raw_data[, 4]   # Twin 2, Trait 2
  )

  return(df)
}

simulate_acde_dataset_enhanced <- function(n_pairs = 400, n_pool = 20,
                                         stratified_prop = 0.6,
                                         target_bands = list(high = c(0.5, 0.8),
                                                            moderate = c(0.2, 0.5),
                                                            low = c(0.0, 0.2)),
                                         max_correlation = 0.8,
                                         rebalance = TRUE,
                                         verbose = FALSE,

                                         # NEW: ACDE component constraints
                                         a_range = c(0.01, 0.64),
                                         c_range = c(0.00, 0.49),
                                         d_range = c(0.00, 0.36),
                                         e_range = c(0.04, 0.81),
                                         min_variance = 0.01,
                                         allow_negative_rho_D = FALSE) {
  # Enhanced simulation function with mathematical constraints

  if (verbose) {
    cat("🧬 Generating enhanced ACDE parameters...\n")
    cat(sprintf("   Stratified sampling: %.1f%% of datasets\n", stratified_prop * 100))
    cat(sprintf("   Maximum correlation: ±%.2f\n", max_correlation))
    cat(sprintf("   Trait pool size: %d\n", n_pool))
    cat("   Component variance ranges:\n")
    cat(sprintf("     A: [%.4f, %.4f]\n", a_range[1], a_range[2]))
    cat(sprintf("     C: [%.4f, %.4f]\n", c_range[1], c_range[2]))
    cat(sprintf("     D: [%.4f, %.4f]\n", d_range[1], d_range[2]))
    cat(sprintf("     E: [%.4f, %.4f]\n", e_range[1], e_range[2]))
    cat(sprintf("   Minimum variance threshold: %.4f\n", min_variance))
    cat(sprintf("   Allow negative rho_D: %s\n", allow_negative_rho_D))
  }

  # Generate random ACDE parameters with enhanced strategy and constraints
  acde_params <- generate_acde_params_enhanced(
    iteration = 1, total_iterations = 1,
    n_pool = n_pool, stratified_prop = stratified_prop,
    target_bands = target_bands, max_correlation = max_correlation,

    # Pass constraint parameters
    a_range = a_range, c_range = c_range,
    d_range = d_range, e_range = e_range,
    min_variance = min_variance,
    allow_negative_rho_D = allow_negative_rho_D
  )

  if (verbose) {
    cat("📊 Building covariance matrices...\n")
  }

  # Build covariance matrices
  covariances <- build_acde_covariances(acde_params)

  # Define twin groups with metadata
  groups <- list(
    MZf = list(zyg = 1, zyg2 = 1, sex1 = 2, sex2 = 2, sigma = covariances$Sigma_MZ, n_frac = 0.2),
    DZf = list(zyg = 2, zyg2 = 2, sex1 = 2, sex2 = 2, sigma = covariances$Sigma_DZ, n_frac = 0.2),
    MZm = list(zyg = 3, zyg2 = 1, sex1 = 1, sex2 = 1, sigma = covariances$Sigma_MZ, n_frac = 0.2),
    DZm = list(zyg = 4, zyg2 = 2, sex1 = 1, sex2 = 1, sigma = covariances$Sigma_DZ, n_frac = 0.2),
    DZo = list(zyg = 5, zyg2 = 2, sex1 = 1, sex2 = 2, sigma = covariances$Sigma_DZ, n_frac = 0.2)
  )

  if (verbose) {
    cat("👥 Simulating twin groups...\n")
  }

  # Simulate each group
  all_data <- list()
  for (group_name in names(groups)) {
    group_info <- groups[[group_name]]
    group_n <- round(n_pairs * group_info$n_frac)

    if (verbose) {
      cat(sprintf("   %s: %d pairs\n", group_name, group_n))
    }

    group_data <- simulate_twin_group(group_n, group_info$sigma, group_info)
    all_data[[group_name]] <- group_data
  }

  # Combine all groups
  final_data <- do.call(rbind, all_data)
  rownames(final_data) <- NULL

  # Attach truth values as attributes
  paths <- acde_params$paths
  cors <- acde_params$correlations

  # Truth values in expected format
  truth <- list(
    a11 = paths["a1"]^2,  # Convert paths to variances
    a22 = paths["a2"]^2,
    a12 = cors["rho_A"] * paths["a1"] * paths["a2"],  # Cross-trait covariance
    c11 = paths["c1"]^2,
    c22 = paths["c2"]^2,
    c12 = cors["rho_C"] * paths["c1"] * paths["c2"],
    d11 = paths["d1"]^2,
    d22 = paths["d2"]^2,
    d12 = cors["rho_D"] * paths["d1"] * paths["d2"],
    e11 = paths["e1"]^2,
    e22 = paths["e2"]^2,
    e12 = cors["rho_E"] * paths["e1"] * paths["e2"]
  )

  attr(final_data, ".truth") <- truth
  attr(final_data, ".covariances") <- covariances
  attr(final_data, ".params") <- acde_params

  if (verbose) {
    cat("\n✅ Enhanced simulation complete!\n")
    cat(sprintf("   Dataset: %d rows × %d columns\n", nrow(final_data), ncol(final_data)))
    cat(sprintf("   A components: %.3f, %.3f (r = %.3f)\n",
                truth$a11, truth$a22, cors["rho_A"]))
    cat(sprintf("   C components: %.3f, %.3f (r = %.3f)\n",
                truth$c11, truth$c22, cors["rho_C"]))
    cat(sprintf("   D components: %.3f, %.3f (r = %.3f)\n",
                truth$d11, truth$d22, cors["rho_D"]))
    cat(sprintf("   E components: %.3f, %.3f (r = %.3f)\n",
                truth$e11, truth$e22, cors["rho_E"]))

    # Show trait pairing metadata if available
    if (!is.null(acde_params$meta)) {
      cat(sprintf("   Implied phenotypic r: %.3f\n", acde_params$meta$implied_phenotypic_corr))
      cat(sprintf("   Sampling strategy: %s\n", acde_params$meta$sampling_strategy))
      if (!is.null(acde_params$meta$target_band) && acde_params$meta$target_band != "random") {
        cat(sprintf("   Target band: %s\n", acde_params$meta$target_band))
      }

      # Show constraint compliance
      if (!is.null(acde_params$meta$constraint_info)) {
        cat("   Constraint compliance: ✓\n")
        if (!allow_negative_rho_D && cors["rho_D"] < 0) {
          cat("   Warning: Negative rho_D generated despite constraint\n")
        }
      }
    }
  }

  return(final_data)
}

# ==============================================================================
# BATCH SIMULATION WITH REBALANCING
# ==============================================================================

simulate_acde_batch <- function(n_datasets = 100, n_pairs = 400,
                               rebalance = TRUE, verbose = TRUE, ...) {
  # Simulate multiple datasets with optional rebalancing

  if (verbose) {
    cat("🎯 BATCH ACDE SIMULATION WITH REBALANCING\n")
    cat("==========================================\n")
    cat(sprintf("Target: %d datasets\n", n_datasets))
    cat(sprintf("Rebalancing: %s\n\n", ifelse(rebalance, "enabled", "disabled")))
  }

  datasets <- list()
  phenotypic_corrs <- numeric(n_datasets)

  # Phase 1: Initial generation
  if (verbose) cat("Phase 1: Initial dataset generation\n")

  for (i in 1:n_datasets) {
    if (verbose && i %% 10 == 0) {
      cat(sprintf("   Generated %d/%d datasets\n", i, n_datasets))
    }

    # Generate dataset with iteration-aware strategy
    data <- simulate_acde_dataset_enhanced(
      n_pairs = n_pairs,
      verbose = FALSE,
      ...
    )

    datasets[[i]] <- data

    # Extract phenotypic correlation
    params <- attr(data, ".params")
    if (!is.null(params$meta$implied_phenotypic_corr)) {
      phenotypic_corrs[i] <- params$meta$implied_phenotypic_corr
    }
  }

  if (!rebalance) {
    if (verbose) cat("\n✅ Batch simulation complete (no rebalancing)\n")
    return(list(datasets = datasets, phenotypic_correlations = phenotypic_corrs))
  }

  # Phase 2: Analyze correlation distribution
  if (verbose) cat("\nPhase 2: Analyzing correlation distribution\n")

  correlation_analysis <- analyze_correlation_distribution(phenotypic_corrs)

  if (verbose) {
    cat(sprintf("   Empirical range: [%.3f, %.3f]\n",
                correlation_analysis$empirical_range[1],
                correlation_analysis$empirical_range[2]))
    cat(sprintf("   Mean correlation: %.3f (SD: %.3f)\n",
                correlation_analysis$empirical_mean,
                correlation_analysis$empirical_sd))
    cat("   Band coverage:\n")
    for (band_name in names(correlation_analysis$band_proportions)) {
      prop <- correlation_analysis$band_proportions[[band_name]]
      count <- correlation_analysis$band_counts[[band_name]]
      cat(sprintf("     %s: %.1f%% (%d datasets)\n",
                  band_name, prop * 100, count))
    }
  }

  # Phase 3: Identify and resample under-represented regions
  undersample_bands <- identify_undersample_regions(correlation_analysis)

  if (length(undersample_bands) > 0 && verbose) {
    cat("\nPhase 3: Rebalancing under-sampled regions\n")
    cat(sprintf("   Under-sampled bands: %s\n", paste(undersample_bands, collapse = ", ")))
  }

  additional_datasets <- 0
  for (band_name in undersample_bands) {
    current_count <- correlation_analysis$band_counts[[band_name]]
    target_count <- max(20, round(n_datasets * 0.2))  # At least 20% representation
    needed <- target_count - current_count

    if (needed > 0) {
      if (verbose) {
        cat(sprintf("   Generating %d additional datasets for %s band\n",
                    needed, band_name))
      }

      # Get target band definition
      target_bands <- list(high = c(0.5, 0.8), moderate = c(0.2, 0.5), low = c(0.0, 0.2))
      target_band <- target_bands[[band_name]]

      for (j in 1:needed) {
        # Generate dataset specifically targeting this band
        data <- simulate_acde_dataset_enhanced(
          n_pairs = n_pairs,
          stratified_prop = 1.0,  # Force stratified sampling
          target_bands = setNames(list(target_band), band_name),
          verbose = FALSE,
          ...
        )

        datasets[[length(datasets) + 1]] <- data

        # Update correlation tracking
        params <- attr(data, ".params")
        if (!is.null(params$meta$implied_phenotypic_corr)) {
          phenotypic_corrs <- c(phenotypic_corrs, params$meta$implied_phenotypic_corr)
        }

        additional_datasets <- additional_datasets + 1
      }
    }
  }

  if (verbose) {
    if (additional_datasets > 0) {
      cat(sprintf("\n✅ Rebalancing complete! Added %d datasets\n", additional_datasets))
      cat(sprintf("   Total datasets: %d\n", length(datasets)))
    } else {
      cat("\n✅ No rebalancing needed - good coverage achieved\n")
    }
  }

  return(list(
    datasets = datasets,
    phenotypic_correlations = phenotypic_corrs,
    correlation_analysis = correlation_analysis,
    additional_datasets = additional_datasets
  ))
}

# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

# Function to extract truth values from simulated data
get_truth <- function(data) {
  return(attr(data, ".truth"))
}

# Function to extract covariance matrices
get_covariances <- function(data) {
  return(attr(data, ".covariances"))
}

# Function to extract original parameters
get_params <- function(data) {
  return(attr(data, ".params"))
}

# Function to verify data structure matches expected format
verify_data_structure <- function(data) {
  expected_vars <- c('zyg', 'zyg2', 'sex1', 'sex2', 'p1_t1', 'p1_t2', 'p2_t1', 'p2_t2')

  missing_vars <- setdiff(expected_vars, names(data))
  extra_vars <- setdiff(names(data), expected_vars)

  if (length(missing_vars) > 0) {
    warning("Missing expected variables: ", paste(missing_vars, collapse = ", "))
  }

  if (length(extra_vars) > 0) {
    warning("Extra variables found: ", paste(extra_vars, collapse = ", "))
  }

  # Check zygosity groups
  zyg_counts <- table(data$zyg)
  expected_zyg <- 1:5

  if (!all(expected_zyg %in% names(zyg_counts))) {
    missing_zyg <- setdiff(expected_zyg, names(zyg_counts))
    warning("Missing zygosity groups: ", paste(missing_zyg, collapse = ", "))
  }

  return(invisible(TRUE))
}

# ==============================================================================
# MAIN EXECUTION FUNCTION
# ==============================================================================

main <- function() {
  # Parse command line arguments
  params <- parse_args()

  # Set seed if provided
  if (!is.null(params$seed)) {
    set.seed(params$seed)
  }

  if (params$verbose) {
    cat("🧬 ENHANCED ACDE TWIN DATA SIMULATION\n")
    cat("=====================================\n")
    cat(sprintf("Sample size: %d twin pairs\n", params$n_pairs))
    cat(sprintf("Trait pool size: %d\n", params$n_pool))
    cat(sprintf("Max correlation: ±%.2f\n", params$max_correlation))
    cat(sprintf("Stratified sampling: %.1f%%\n", params$stratified_prop * 100))
    cat(sprintf("Rebalancing: %s\n", ifelse(params$rebalance, "enabled", "disabled")))
    if (!is.null(params$seed)) {
      cat(sprintf("Random seed: %d\n", params$seed))
    }
    if (!is.null(params$outfile)) {
      cat(sprintf("Output file: %s\n", params$outfile))
    }
    cat("\n")
  }

  # Generate data
  simulated_data <- simulate_acde_dataset_enhanced(
    n_pairs = params$n_pairs,
    n_pool = params$n_pool,
    stratified_prop = params$stratified_prop,
    target_bands = params$target_bands,
    max_correlation = params$max_correlation,
    rebalance = params$rebalance,
    verbose = params$verbose
  )

  # Verify structure
  if (params$verbose) {
    cat("\n🔍 Verifying data structure...\n")
    verify_data_structure(simulated_data)
  }

  # Save or return data
  if (!is.null(params$outfile)) {
    # Create output directory if needed
    outdir <- dirname(params$outfile)
    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

    # Write CSV
    write.csv(simulated_data, params$outfile, row.names = FALSE, quote = FALSE)

    if (params$verbose) {
      cat(sprintf("\n💾 Data saved to: %s\n", params$outfile))
    }

    return(invisible(params$outfile))
  } else {
    return(simulated_data)
  }
}

# ==============================================================================
# EXAMPLE USAGE AND TESTING
# ==============================================================================

# Function to demonstrate enhanced usage
demo_enhanced_usage <- function() {
  cat("🎯 ENHANCED ACDE SIMULATION DEMO\n")
  cat("=================================\n\n")

  # Generate example dataset with enhanced features
  data <- simulate_acde_dataset_enhanced(
    n_pairs = 200,
    stratified_prop = 0.7,
    max_correlation = 0.8,
    verbose = TRUE
  )

  # Show structure
  cat("\n📋 Data structure:\n")
  print(str(data))

  cat("\n📊 First few rows:\n")
  print(head(data))

  # Show enhanced truth values
  cat("\n🎯 True parameter values:\n")
  truth <- get_truth(data)
  for (name in names(truth)) {
    cat(sprintf("   %s: %.4f\n", name, truth[[name]]))
  }

  # Show enhanced metadata
  cat("\n🔬 Enhanced metadata:\n")
  params <- get_params(data)
  if (!is.null(params$meta)) {
    meta <- params$meta
    cat(sprintf("   Implied phenotypic correlation: %.4f\n", meta$implied_phenotypic_corr))
    cat(sprintf("   Sampling strategy: %s\n", meta$sampling_strategy))
    cat(sprintf("   Target band: %s\n", ifelse(is.null(meta$target_band), "random", meta$target_band)))
    cat(sprintf("   Pool search attempts: %d\n", meta$retry_count))
    cat(sprintf("   Trait pool size: %d\n", meta$pool_size))
  }

  # Test batch simulation with rebalancing
  cat("\n🔄 Testing batch simulation (small example):\n")
  batch_result <- simulate_acde_batch(
    n_datasets = 20,
    n_pairs = 100,
    rebalance = TRUE,
    verbose = TRUE,
    max_correlation = 0.7
  )

  cat("\n📈 Batch results:\n")
  cat(sprintf("   Generated datasets: %d\n", length(batch_result$datasets)))
  cat(sprintf("   Additional datasets from rebalancing: %d\n", batch_result$additional_datasets))
  cat(sprintf("   Phenotypic correlation range: [%.3f, %.3f]\n",
              min(batch_result$phenotypic_correlations, na.rm = TRUE),
              max(batch_result$phenotypic_correlations, na.rm = TRUE)))

  return(invisible(list(single = data, batch = batch_result)))
}

# Function to analyze correlation coverage across multiple datasets
analyze_simulation_coverage <- function(n_datasets = 100, n_pairs = 400, ...) {
  cat("📊 CORRELATION COVERAGE ANALYSIS\n")
  cat("=================================\n")

  batch_result <- simulate_acde_batch(
    n_datasets = n_datasets,
    n_pairs = n_pairs,
    rebalance = TRUE,
    verbose = TRUE,
    ...
  )

  # Extract all correlations
  pheno_corrs <- batch_result$phenotypic_correlations

  # Create detailed analysis
  cat("\n📈 Detailed Coverage Report:\n")
  cat(sprintf("Total datasets: %d\n", length(batch_result$datasets)))
  cat(sprintf("Correlation range: [%.3f, %.3f]\n", min(pheno_corrs, na.rm = TRUE), max(pheno_corrs, na.rm = TRUE)))
  cat(sprintf("Mean correlation: %.3f (±%.3f)\n", mean(pheno_corrs, na.rm = TRUE), sd(pheno_corrs, na.rm = TRUE)))

  # Histogram-like analysis
  breaks <- seq(-1, 1, by = 0.2)
  hist_counts <- hist(pheno_corrs, breaks = breaks, plot = FALSE)$counts
  cat("\nCorrelation distribution:\n")
  for (i in 1:length(hist_counts)) {
    range_str <- sprintf("[%.1f, %.1f)", breaks[i], breaks[i+1])
    bar <- paste(rep("█", round(hist_counts[i] / max(hist_counts) * 20)), collapse = "")
    cat(sprintf("   %s: %3d %s\n", range_str, hist_counts[i], bar))
  }

  return(invisible(batch_result))
}

# ==============================================================================
# EXECUTE WHEN RUN AS SCRIPT
# ==============================================================================

if (!interactive()) {
  result <- main()

  # If running from command line and no output file specified,
  # show a sample of the data
  params <- parse_args()
  if (is.null(params$outfile) && params$verbose) {
    cat("\n📋 Sample of generated data:\n")
    if (is.data.frame(result)) {
      print(head(result))
      cat(sprintf("\nTotal rows: %d\n", nrow(result)))

      # Show correlation info
      params_data <- get_params(result)
      if (!is.null(params_data$meta$implied_phenotypic_corr)) {
        cat(sprintf("Phenotypic correlation: %.4f\n", params_data$meta$implied_phenotypic_corr))
      }
    }
  }
}

# ==============================================================================
# NOTES FOR USERS
# ==============================================================================

# COMMAND LINE USAGE:
# Rscript engine.R --n_pairs 500 --outfile data.csv --max_corr 0.8 --seed 12345
# Rscript engine.R --n_pairs 1000 --n_pool 30 --no_rebalance --quiet
#
# FUNCTION USAGE:
# # Single dataset
# data <- simulate_acde_dataset_enhanced(n_pairs = 400, max_correlation = 0.8)
#
# # Batch with rebalancing
# batch <- simulate_acde_batch(n_datasets = 100, rebalance = TRUE)
#
# # Coverage analysis
# coverage <- analyze_simulation_coverage(n_datasets = 200)
#
# ENHANCED FEATURES:
# - Stratified sampling across correlation bands (high/moderate/low)
# - Expanded correlation ranges (up to ±0.8 for genetic, ±0.5 for environmental)
# - Automatic rebalancing for under-sampled correlation regions
# - Enhanced metadata tracking sampling strategy and correlation targets
# - Batch simulation with coverage analysis
#
# CORRELATION RANGES:
# - rho_A (genetic): ±0.8
# - rho_C (shared environment): ±0.5
# - rho_D (dominance): ±0.5
# - rho_E (unique environment): ±0.4
#
# VARIABLE NAMES (unchanged for compatibility):
# allVars <- c('zyg','zyg2','sex1','sex2','p1_t1','p1_t2','p2_t1','p2_t2')
