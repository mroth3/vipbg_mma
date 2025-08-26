# ==============================================================================
# UNIVARIATE ACDE MULTIMODEL FITTING - OpenMx Implementation (MZ/DZ Only)
# ==============================================================================
# Fits univariate twin models under ACDE framework with mxModelAverage
# Compatible with bivariate simulator output: p1_t1,p2_t1,p1_t2,p2_t2

library(OpenMx)

# ==============================================================================
# USER-TUNABLE CONTROLS
# ==============================================================================

target_trait <- 1L        # 1 or 2; which trait to analyze (p1_t1/p2_t1 or p1_t2/p2_t2)
useTryHard   <- TRUE      # Use mxTryHard instead of mxRun for robustness
optimizer    <- "NPSOL"   # Optimizer: "NPSOL", "SLSQP", or "CSOLNP"
startVals    <- c(A = 0.3, C = 0.3, E = 0.3)  # Starting values for variance components
includeOSDZ  <- TRUE     # Include opposite-sex DZ twins in DZ pool
includeACEc  <- FALSE     # Include constrained ACE model (e fixed)
zscore       <- FALSE     # Z-score data within groups before fitting
smartStarts  <- FALSE     # Use observed correlations for starting values
logResults   <- TRUE      # Log fit statistics to CSV file

# Lower bound variant controls
lb_threshold <- 0.05     # Only generate variants if component estimated below this
lb_start     <- 0.01     # Always start lower bounds at this value
lb_upper     <- 0.3      # Upper limit for lower bound variants
lb_step      <- 0.02      # Step size for lower bound variants

# Set optimizer with fallback
availOpt <- tryCatch(mxOption(NULL, "Default optimizer", optimizer), error = function(e) e)
if (inherits(availOpt, "error")) {
  mxOption(NULL, "Default optimizer", "SLSQP")
  cat("Requested optimizer", optimizer, "not available; falling back to SLSQP\n")
}

# Print configuration settings
cat("=== CONFIGURATION SETTINGS ===\n")
cat("Target trait:", target_trait, "\n")
cat("Use TryHard:", useTryHard, "\n")
cat("Optimizer:", optimizer, "\n")
cat("Starting values - A:", startVals[["A"]], "C:", startVals[["C"]], "E:", startVals[["E"]], "\n")
cat("Include OSDZ:", includeOSDZ, "\n")
cat("Include ACEc:", includeACEc, "\n")
cat("Z-score data:", zscore, "\n")
cat("Smart starts:", smartStarts, "\n")
cat("Log results:", logResults, "\n")
cat("Lower bound threshold:", lb_threshold, "\n")
cat("Lower bound start:", lb_start, "\n")
cat("Lower bound upper limit:", lb_upper, "\n")
cat("Lower bound step size:", lb_step, "\n")

# ==============================================================================
# DATA PREPARATION
# ==============================================================================
data <- dataTwo
# Check for required columns
required_cols <- c("p1_t1", "p2_t1", "p1_t2", "p2_t2", "zyg")
missing_cols <- setdiff(required_cols, names(data))
if(length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Select trait-specific columns
if(target_trait == 1L) {
  selVars <- c("p1_t1", "p2_t1")
} else if(target_trait == 2L) {
  selVars <- c("p1_t2", "p2_t2")
} else {
  stop("target_trait must be 1 or 2")
}

# Extract twin data subsets by zygosity (pooled across sex, no sex differences)
dataMZ <- subset(data, zyg %in% c(1,3), select = c(selVars, "zyg"))  # MZ twins
dataDZ <- subset(data, zyg %in% c(2,4), select = c(selVars, "zyg"))  # DZ twins

# Optionally include opposite-sex DZ twins
if(includeOSDZ && 5 %in% data$zyg) {
  dataDZos <- subset(data, zyg == 5, select = c(selVars, "zyg"))
  dataDZ <- rbind(dataDZ, dataDZos)
  cat("Including opposite-sex DZ twins in DZ pool (n =", nrow(dataDZos), ")\n")
}

cat("Using trait variables:", paste(selVars, collapse = ", "), "\n")
cat("MZ sample size:", nrow(dataMZ), "\n")
cat("DZ sample size:", nrow(dataDZ), "\n")

# Remove rows with missing data and check for empty groups
dataMZ <- dataMZ[complete.cases(dataMZ[, selVars]), ]
dataDZ <- dataDZ[complete.cases(dataDZ[, selVars]), ]

# Guard against empty groups after cleaning
if(nrow(dataMZ) == 0L) stop("No MZ rows after filtering/cleaning for selected trait.")
if(nrow(dataDZ) == 0L) stop("No DZ rows after filtering/cleaning for selected trait.")

cat("After removing missing data - MZ:", nrow(dataMZ), "DZ:", nrow(dataDZ), "\n")

# Optional z-scoring within groups
if(zscore) {
  cat("Applying z-score standardization within groups\n")
  for(v in selVars) {
    dataMZ[[v]] <- scale(dataMZ[[v]])[,1]
    dataDZ[[v]] <- scale(dataDZ[[v]])[,1]
  }
}

# Calculate observed correlations for reference and smart starts
mz_cor <- cor(dataMZ[, selVars], use = "complete.obs")[1,2]
dz_cor <- cor(dataDZ[, selVars], use = "complete.obs")[1,2]

# Guard smart starts against non-finite correlations
if(smartStarts && (!is.finite(mz_cor) || !is.finite(dz_cor))) {
  warning("smartStarts disabled: non-finite MZ/DZ correlations.")
  smartStarts <- FALSE
}

cat("\n=== OBSERVED TWIN CORRELATIONS ===\n")
cat("MZ correlation:", round(mz_cor, 3), "\n")
cat("DZ correlation:", round(dz_cor, 3), "\n")

# Data objects for model fitting
datMZ <- mxData(observed = dataMZ[, selVars], type = "raw")
datDZ <- mxData(observed = dataDZ[, selVars], type = "raw")

# ==============================================================================
# COMMON MODEL COMPONENTS (no name/label collisions)
# ==============================================================================

# Mean structure (equal across the two twins by using the same label 'm')
Mean <- mxMatrix(
  type = "Full", nrow = 1, ncol = 2,
  free = TRUE, values = 0,
  labels = c("m","m"),              # labels tie the two means; name is 'Mean'
  name = "Mean"
)

# Expected covariance algebras (added into each submodel, so A/C/D/E resolve locally)
mzCov <- mxAlgebra(
  rbind(
    cbind(A + C + D + E, A + C + D),
    cbind(A + C + D,     A + C + D + E)
  ),
  name = "mzCov"
)

dzCov <- mxAlgebra(
  rbind(
    cbind(A + C + D + E, 0.5*A + C + 0.25*D),
    cbind(0.5*A + C + 0.25*D, A + C + D + E)
  ),
  name = "dzCov"
)

# Expectations and fit functions
expMZ <- mxExpectationNormal(covariance = "mzCov", means = "Mean", dimnames = selVars)
expDZ <- mxExpectationNormal(covariance = "dzCov", means = "Mean", dimnames = selVars)
fitML <- mxFitFunctionML()

# Multigroup fit function
multiTwin <- mxFitFunctionMultigroup(c("MZ","DZ"))
# ==============================================================================
# HELPER FUNCTION TO BUILD MODELS WITH SMART STARTS
# ==============================================================================

buildTwinModel <- function(modelName,
                           a_free = FALSE, c_free = FALSE, d_free = FALSE, e_free = TRUE,
                           a_lbound = NULL, c_lbound = NULL, d_lbound = NULL, e_lbound = NULL) {

  # Set starting values from user config or smart starts
  a_val <- startVals[["A"]]
  c_val <- startVals[["C"]]
  d_val <- startVals[["A"]]  # start D like A unless user sets otherwise
  e_val <- startVals[["E"]]

  # Optional smart starts
  if (smartStarts) {
    if (modelName == "AE" && a_free && e_free && !c_free && !d_free) {
      smart_a <- sqrt(pmax(mz_cor - dz_cor, 0.01))
      smart_e <- sqrt(pmax(1 - mz_cor, 0.01))
      if (all(is.finite(c(smart_a, smart_e)))) { a_val <- smart_a; e_val <- smart_e }
    } else if (modelName == "ACE" && a_free && c_free && e_free && !d_free) {
      smart_c <- sqrt(pmax(2*dz_cor - mz_cor, 0))
      smart_a <- sqrt(pmax(mz_cor - smart_c^2, 0))
      smart_e <- sqrt(pmax(1 - smart_a^2 - smart_c^2, 0.01))
      if (all(is.finite(c(smart_a, smart_c, smart_e)))) { a_val <- smart_a; c_val <- smart_c; e_val <- smart_e }
    } else if (modelName == "ADE" && a_free && d_free && e_free && !c_free) {
      smart_d <- sqrt(pmax(2*mz_cor - 2*dz_cor, 0))
      smart_a <- sqrt(pmax(2*dz_cor - smart_d^2, 0))
      smart_e <- sqrt(pmax(1 - smart_a^2 - smart_d^2, 0.01))
      if (all(is.finite(c(smart_a, smart_d, smart_e)))) { a_val <- smart_a; d_val <- smart_d; e_val <- smart_e }
    }
  }

  # Set lower bounds (use provided or default)
  a_lb <- if(a_free) (if(is.null(a_lbound)) 0 else a_lbound) else NA
  c_lb <- if(c_free) (if(is.null(c_lbound)) 0 else c_lbound) else NA
  d_lb <- if(d_free) (if(is.null(d_lbound)) 0 else d_lbound) else NA
  e_lb <- if(e_free) (if(is.null(e_lbound)) 1e-6 else e_lbound) else NA

  # ---- Path coefficient matrices (labels MUST NOT equal object names) ----
    a_mat <- mxMatrix(type="Lower", nrow=1, ncol=1, free=a_free,
                  values = ifelse(a_free, a_val, 0),
                  labels = "a11", lbound = a_lb, name="a")

    c_mat <- mxMatrix(type="Lower", nrow=1, ncol=1, free=c_free,
                  values = ifelse(c_free, c_val, 0),
                  labels = "c11", lbound = c_lb, name="c")

    d_mat <- mxMatrix(type="Lower", nrow=1, ncol=1, free=d_free,
                  values = ifelse(d_free, d_val, 0),
                  labels = "d11", lbound = d_lb, name="d")

    e_mat <- mxMatrix(type="Lower", nrow=1, ncol=1, free=e_free,
                  values = ifelse(e_free, max(e_val, 0.1), 0),
                  labels = "e11", lbound = e_lb, name="e")
  # Variance component algebras
  A_alg <- mxAlgebra(a %*% t(a), name="A")
  C_alg <- mxAlgebra(c %*% t(c), name="C")
  D_alg <- mxAlgebra(d %*% t(d), name="D")
  E_alg <- mxAlgebra(e %*% t(e), name="E")

  # Raw component scalars (for logging/model averaging)
  a2_alg <- mxAlgebra(A[1,1], name="a2")
  c2_alg <- mxAlgebra(C[1,1], name="c2")
  d2_alg <- mxAlgebra(D[1,1], name="d2")
  e2_alg <- mxAlgebra(E[1,1], name="e2")

  # Submodels
    modelMZ <- mxModel("MZ", a_mat, c_mat, d_mat, e_mat,
                   A_alg, C_alg, D_alg, E_alg,
                   a2_alg, c2_alg, d2_alg, e2_alg,
                   Mean, mzCov, datMZ, expMZ, fitML)

    modelDZ <- mxModel("DZ", a_mat, c_mat, d_mat, e_mat,
                   A_alg, C_alg, D_alg, E_alg,
                   a2_alg, c2_alg, d2_alg, e2_alg,
                   Mean, dzCov, datDZ, expDZ, fitML)

  # Top-level: standardized proportions (use namespaced MZ.*)
  Tot  <- mxAlgebra(MZ.a2 + MZ.c2 + MZ.d2 + MZ.e2, name="Tot")
  Prop <- mxAlgebra(cbind(MZ.a2/Tot, MZ.c2/Tot, MZ.d2/Tot, MZ.e2/Tot),
                    name="Prop",
                    dimnames=list(NULL, c("a2","c2","d2","e2")))

  mxModel(modelName, modelMZ, modelDZ, Tot, Prop, multiTwin)
}

# ==============================================================================
# FUNCTION TO LOG MODEL RESULTS
# ==============================================================================

logModelResults <- function(fit, results_list) {
  # If fit is not a valid MxModel, log a failure row and return
  if (!inherits(fit, "MxModel") || is.null(fit$output)) {
    result_row <- data.frame(
      Model      = if (!is.null(fit) && !is.null(fit$name)) fit$name else NA_character_,
      MZ_N       = if (exists("dataMZ")) nrow(dataMZ) else NA_integer_,
      DZ_N       = if (exists("dataDZ")) nrow(dataDZ) else NA_integer_,
      MZ_cor     = if (exists("mz_cor")) round(mz_cor, 4) else NA_real_,
      DZ_cor     = if (exists("dz_cor")) round(dz_cor, 4) else NA_real_,
      minus2LL   = NA_real_,
      AIC        = NA_real_,
      BIC        = NA_real_,
      df         = NA_integer_,
      status     = NA_integer_,
      a_param    = NA_real_, c_param = NA_real_, d_param = NA_real_, e_param = NA_real_,
      A_var      = NA_real_, C_var = NA_real_, D_var = NA_real_, E_var = NA_real_,
      Total_var  = NA_real_,
      stringsAsFactors = FALSE
    )
    return(rbind(results_list, result_row))
  }

  # Safe extractors ------------------------------------------------------------
  safe_eval <- function(expr) {
    out <- try(expr, silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
  }

  # Fit stats
  minus2LL <- if (!is.null(fit$output$Minus2LogLikelihood)) fit$output$Minus2LogLikelihood else NA_real_
  AIC      <- if (!is.null(fit$output$AIC))                fit$output$AIC                else NA_real_
  BIC      <- if (!is.null(fit$output$BIC))                fit$output$BIC                else NA_real_
  df       <- if (!is.null(fit$output$degreesOfFreedom))   fit$output$degreesOfFreedom   else NA_integer_
  status   <- if (!is.null(fit$output$status$code))        fit$output$status$code        else NA_integer_

  # Parameter values (only if those matrices exist)
  a_param <- if ("a" %in% names(fit$MZ)) safe_eval(mxEval(a, fit)) else NA_real_
  c_param <- if ("c" %in% names(fit$MZ)) safe_eval(mxEval(c, fit)) else NA_real_
  d_param <- if ("d" %in% names(fit$MZ)) safe_eval(mxEval(d, fit)) else NA_real_
  e_param <- if ("e" %in% names(fit$MZ)) safe_eval(mxEval(e, fit)) else NA_real_

  # Variance components (as defined in your top model via child MZ algebras)
  a2_val <- safe_eval(mxEval(MZ.a2, fit))
  c2_val <- safe_eval(mxEval(MZ.c2, fit))
  d2_val <- safe_eval(mxEval(MZ.d2, fit))
  e2_val <- safe_eval(mxEval(MZ.e2, fit))
  total_var <- sum(c(a2_val, c2_val, d2_val, e2_val), na.rm = TRUE)

  # Assemble result row
  result_row <- data.frame(
    Model     = fit$name,
    MZ_N      = if (exists("dataMZ")) nrow(dataMZ) else NA_integer_,
    DZ_N      = if (exists("dataDZ")) nrow(dataDZ) else NA_integer_,
    MZ_cor    = if (exists("mz_cor")) round(mz_cor, 4) else NA_real_,
    DZ_cor    = if (exists("dz_cor")) round(dz_cor, 4) else NA_real_,
    minus2LL  = round(minus2LL, 4),
    AIC       = round(AIC, 4),
    BIC       = round(BIC, 4),
    df        = df,
    status    = status,
    a_param   = if (!is.na(a_param)) round(a_param, 4) else NA_real_,
    c_param   = if (!is.na(c_param)) round(c_param, 4) else NA_real_,
    d_param   = if (!is.na(d_param)) round(d_param, 4) else NA_real_,
    e_param   = if (!is.na(e_param)) round(e_param, 4) else NA_real_,
    A_var     = if (!is.na(a2_val))  round(a2_val, 4)  else NA_real_,
    C_var     = if (!is.na(c2_val))  round(c2_val, 4)  else NA_real_,
    D_var     = if (!is.na(d2_val))  round(d2_val, 4)  else NA_real_,
    E_var     = if (!is.na(e2_val))  round(e2_val, 4)  else NA_real_,
    Total_var = if (is.finite(total_var)) round(total_var, 4) else NA_real_,
    stringsAsFactors = FALSE
  )

  rbind(results_list, result_row)
}
# ==============================================================================
# SANITY CHECK FUNCTION
# ==============================================================================

checkVarianceComponents <- function(fit, verbose = TRUE, tol = 1e-8) {
  # If fit failed, return NAs
  if (!inherits(fit, "MxModel") || is.null(fit$output)) {
    if (verbose) cat("WARNING: fit is not a valid MxModel; returning NA proportions.\n")
    return(setNames(rep(NA_real_, 4), c("A","C","D","E")))
  }

  safe_eval <- function(expr) {
    out <- try(expr, silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
  }

  a2 <- safe_eval(mxEval(MZ.a2, fit))
  c2 <- safe_eval(mxEval(MZ.c2, fit))
  d2 <- safe_eval(mxEval(MZ.d2, fit))
  e2 <- safe_eval(mxEval(MZ.e2, fit))

  comps <- c(A = a2, C = c2, D = d2, E = e2)

  if (all(is.na(comps))) {
    if (verbose) cat("WARNING: Could not evaluate variance components; returning NA.\n")
    return(comps)
  }

  total <- sum(comps, na.rm = TRUE)
  if (!is.finite(total) || total < tol) {
    if (verbose) cat("WARNING: Total variance is zero/NA/singular in model", fit$name, "\n")
    return(setNames(rep(NA_real_, 4), names(comps)))
  }

  props <- comps / total

  # Soft range checks
  bad <- which(props < -tol | props > 1 + tol)
  if (length(bad)) {
    if (verbose) {
      cat("WARNING: Variance proportions outside [0,1] in model", fit$name, "\n")
      print(round(props, 4))
    }
  }

  if (verbose && (total < 0.1 || total > 10)) {
    cat("WARNING: Total variance unusual in model", fit$name, ":", round(total, 4), "\n")
  }

  # Clamp tiny numerical drift
  props[pmax(props, 0)] -> props  # ensure non-negative
  props[props < 0 & props > -tol] <- 0
  props[props > 1 & props < 1 + tol] <- 1
  props
}

# ==============================================================================
# LOWER BOUND VARIANT GENERATION FUNCTIONS
# ==============================================================================

# Function to extract variance component estimates from a fitted model
extract_variance_estimates <- function(fit) {
  if (!inherits(fit, "MxModel") || is.null(fit$output)) {
    return(c(A = NA_real_, C = NA_real_, D = NA_real_))
  }

  safe_eval <- function(expr) {
    out <- try(expr, silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else as.numeric(out)
  }

  a2 <- safe_eval(mxEval(MZ.a2, fit))
  c2 <- safe_eval(mxEval(MZ.c2, fit))
  d2 <- safe_eval(mxEval(MZ.d2, fit))

  c(A = a2, C = c2, D = d2)
}

# Function to check which components are free in a model
check_free_components <- function(model_name) {
  base_specs <- list(
    "E"   = c(a_free = FALSE, c_free = FALSE, d_free = FALSE),
    "AE"  = c(a_free = TRUE,  c_free = FALSE, d_free = FALSE),
    "CE"  = c(a_free = FALSE, c_free = TRUE,  d_free = FALSE),
    "ACE" = c(a_free = TRUE,  c_free = TRUE,  d_free = FALSE),
    "DE"  = c(a_free = FALSE, c_free = FALSE, d_free = TRUE),
    "ADE" = c(a_free = TRUE,  c_free = FALSE, d_free = TRUE)
  )

  # Extract base model name (remove lb suffix if present)
  base_name <- gsub("lb\\d+$", "", model_name)

  if (base_name %in% names(base_specs)) {
    return(base_specs[[base_name]])
  } else {
    # Default assumption - no components free
    return(c(a_free = FALSE, c_free = FALSE, d_free = FALSE))
  }
}

# Function to generate lower bound variants for a single model
generate_lb_variants <- function(base_fit, base_model_specs, threshold = 0.05, lb_start = 0.01, lb_upper = 0.3, step = 0.1) {
  if (is.null(base_fit) || !inherits(base_fit, "MxModel") || is.null(base_fit$output)) {
    return(list())
  }

  # Extract variance component estimates
  var_estimates <- extract_variance_estimates(base_fit)

  # Check which components are free using the base_model_specs
  a_free <- base_model_specs$a_free
  c_free <- base_model_specs$c_free
  d_free <- base_model_specs$d_free

  # Find components that are both free and below threshold (indicating potential underestimation)
  components_to_vary <- c()
  if (a_free && is.finite(var_estimates["A"]) && var_estimates["A"] < threshold) {
    components_to_vary <- c(components_to_vary, "A")
  }
  if (c_free && is.finite(var_estimates["C"]) && var_estimates["C"] < threshold) {
    components_to_vary <- c(components_to_vary, "C")
  }
  if (d_free && is.finite(var_estimates["D"]) && var_estimates["D"] < threshold) {
    components_to_vary <- c(components_to_vary, "D")
  }

  if (length(components_to_vary) == 0) {
    return(list())
  }

  cat("Generating lb variants for", base_fit$name, "- components below threshold:",
      paste(components_to_vary, "=", round(var_estimates[components_to_vary], 3), collapse = ", "), "\n")

  # Generate lower bound variants - always start at lb_start, go up to lb_upper
  lb_variants <- list()

  for (comp in components_to_vary) {
    # Generate sequence of lower bound values from lb_start to lb_upper
    lb_values <- seq(lb_start, lb_upper, by = step)
    lb_values <- round(lb_values, 3)
    lb_values <- lb_values[lb_values <= lb_upper]  # Ensure we don't exceed upper limit

    for (lb_val in lb_values) {
      # Create variant name - convert to integer for sprintf
      variant_name <- paste0(base_fit$name, "lb", sprintf("%02d", as.integer(round(lb_val * 100))))

      # Set up lower bounds
      a_lb <- if (comp == "A") lb_val else NULL
      c_lb <- if (comp == "C") lb_val else NULL
      d_lb <- if (comp == "D") lb_val else NULL

      # Build the variant model
      variant_model <- buildTwinModel(
        modelName = variant_name,
        a_free = base_model_specs$a_free,
        c_free = base_model_specs$c_free,
        d_free = base_model_specs$d_free,
        e_free = TRUE,  # E is always free
        a_lbound = a_lb,
        c_lbound = c_lb,
        d_lbound = d_lb
      )

      lb_variants[[variant_name]] <- variant_model
    }
  }

  # If multiple components are below threshold, generate one combined variant
  if (length(components_to_vary) > 1) {
    # Generate one variant where all low components get lb_start as lower bound
    cross_variant_name <- paste0(base_fit$name, "lbAll", sprintf("%02d", as.integer(round(lb_start * 100))))

    a_lb <- if ("A" %in% components_to_vary) lb_start else NULL
    c_lb <- if ("C" %in% components_to_vary) lb_start else NULL
    d_lb <- if ("D" %in% components_to_vary) lb_start else NULL

    cross_variant_model <- buildTwinModel(
      modelName = cross_variant_name,
      a_free = base_model_specs$a_free,
      c_free = base_model_specs$c_free,
      d_free = base_model_specs$d_free,
      e_free = TRUE,
      a_lbound = a_lb,
      c_lbound = c_lb,
      d_lbound = d_lb
    )

    lb_variants[[cross_variant_name]] <- cross_variant_model
  }

  return(lb_variants)
}

# ==============================================================================
# BUILD AND FIT MODELS
# ==============================================================================

# Safe runner + helpers
# -----------------------------
runner <- if (exists("useTryHard") && isTRUE(useTryHard)) {
  # let user set tryHardTries elsewhere; default to 10
  function(m) mxTryHard(m, extraTries = if (exists("tryHardTries")) tryHardTries else 10)
} else {
  function(m) mxRun(m)
}

fit_ok <- function(x) inherits(x, "MxModel") && !is.null(x$output) && is.null(x$output$status$code)

safe_run <- function(model, label) {
  out <- try(runner(model), silent = TRUE)
  if (inherits(out, "try-error") || !inherits(out, "MxModel")) {
    message(sprintf("⚠ %s model failed to return a fitted MxModel. Skipping.", label))
    return(NULL)
  }
  if (!is.null(out$output$status$code) && out$output$status$code != 0) {
    message(sprintf("⚠ %s model status code = %s (non-zero). Proceeding with caution.",
                    label, out$output$status$code))
  }
  out
}

# -----------------------------
# Initialize results logging
# -----------------------------
if (exists("logResults") && isTRUE(logResults)) {
  if (!exists("results_log")) results_log <- data.frame()
}

cat("\n=== FITTING BASE MODELS (univariate twins) ===\n")

fits <- list()  # collect only successful fits

# Define model specifications for reuse
model_specs <- list(
  "E"   = list(a_free = FALSE, c_free = FALSE, d_free = FALSE, e_free = TRUE),
  "AE"  = list(a_free = TRUE,  c_free = FALSE, d_free = FALSE, e_free = TRUE),
  "CE"  = list(a_free = FALSE, c_free = TRUE,  d_free = FALSE, e_free = TRUE),
  "ACE" = list(a_free = TRUE,  c_free = TRUE,  d_free = FALSE, e_free = TRUE),
  "DE"  = list(a_free = FALSE, c_free = FALSE, d_free = TRUE,  e_free = TRUE),
  "ADE" = list(a_free = TRUE,  c_free = FALSE, d_free = TRUE,  e_free = TRUE)
)

# MODEL 1: E (A=0, C=0, D=0)
cat("Fitting E model...\n")
modelE <- buildTwinModel("E", a_free=FALSE, c_free=FALSE, d_free=FALSE, e_free=TRUE)
fitE   <- safe_run(modelE, "E")
if (!is.null(fitE)) {
  checkVarianceComponents(fitE)
  if (exists("logResults") && isTRUE(logResults)) results_log <- logModelResults(fitE, results_log)
  fits$E <- fitE
}

# MODEL 2: AE (C=0, D=0)
cat("Fitting AE model...\n")
modelAE <- buildTwinModel("AE", a_free = TRUE, c_free = FALSE, d_free = FALSE, e_free = TRUE)
fitAE <- safe_run(modelAE, "AE")
if (!is.null(fitAE)) {
  checkVarianceComponents(fitAE)
  if (exists("logResults") && isTRUE(logResults)) results_log <- logModelResults(fitAE, results_log)
  fits$AE <- fitAE
}

# MODEL 3: CE (A=0, D=0)
cat("Fitting CE model...\n")
modelCE <- buildTwinModel("CE", a_free = FALSE, c_free = TRUE, d_free = FALSE, e_free = TRUE)
fitCE <- safe_run(modelCE, "CE")
if (!is.null(fitCE)) {
  checkVarianceComponents(fitCE)
  if (exists("logResults") && isTRUE(logResults)) results_log <- logModelResults(fitCE, results_log)
  fits$CE <- fitCE
}

# MODEL 4: ACE (D=0)
cat("Fitting ACE model...\n")
modelACE <- buildTwinModel("ACE", a_free = TRUE, c_free = TRUE, d_free = FALSE, e_free = TRUE)
fitACE <- safe_run(modelACE, "ACE")
if (!is.null(fitACE)) {
  checkVarianceComponents(fitACE)
  if (exists("logResults") && isTRUE(logResults)) results_log <- logModelResults(fitACE, results_log)
  fits$ACE <- fitACE
}

# MODEL 5: DE (A=0, C=0)
cat("Fitting DE model...\n")
modelDE <- buildTwinModel("DE", a_free = FALSE, c_free = FALSE, d_free = TRUE, e_free = TRUE)
fitDE <- safe_run(modelDE, "DE")
if (!is.null(fitDE)) {
  checkVarianceComponents(fitDE)
  if (exists("logResults") && isTRUE(logResults)) results_log <- logModelResults(fitDE, results_log)
  fits$DE <- fitDE
}

# MODEL 6: ADE (C=0)
cat("Fitting ADE model...\n")
modelADE <- buildTwinModel("ADE", a_free = TRUE, c_free = FALSE, d_free = TRUE, e_free = TRUE)
fitADE <- safe_run(modelADE, "ADE")
if (!is.null(fitADE)) {
  checkVarianceComponents(fitADE)
  if (exists("logResults") && isTRUE(logResults)) results_log <- logModelResults(fitADE, results_log)
  fits$ADE <- fitADE
}

# MODEL 7: ACEc MODEL (optional constrained version)
if (includeACEc) {
  cat("Fitting constrained ACE model (ACEc)...\n")

  ## --- parameter matrices (labels distinct from object names) ---
  a_ACEc <- mxMatrix(type="Lower", nrow=1, ncol=1, free=TRUE,
                     values=startVals[["A"]], labels="par_a", lbound=0, name="aM")
  c_ACEc <- mxMatrix(type="Lower", nrow=1, ncol=1, free=TRUE,
                     values=startVals[["C"]], labels="par_c", lbound=0, name="cM")
  d_ACEc <- mxMatrix(type="Lower", nrow=1, ncol=1, free=FALSE,
                     values=0,                 labels=NA,     name="dM")
  # fixed e for the constrained model
  e_ACEc <- mxMatrix(type="Lower", nrow=1, ncol=1, free=FALSE,
                     values=sqrt(0.4),         labels=NA,     name="eM")

  ## --- variance algebras (within-group namespace) ---
  A_ACEc <- mxAlgebra(aM %*% t(aM), name="A")
  C_ACEc <- mxAlgebra(cM %*% t(cM), name="C")
  D_ACEc <- mxAlgebra(dM %*% t(dM), name="D")
  E_ACEc <- mxAlgebra(eM %*% t(eM), name="E")

  ## --- scalar extracts for logging/standardization ---
  a2_ACEc <- mxAlgebra(A[1,1], name="a2")
  c2_ACEc <- mxAlgebra(C[1,1], name="c2")
  d2_ACEc <- mxAlgebra(D[1,1], name="d2")
  e2_ACEc <- mxAlgebra(E[1,1], name="e2")

  ## --- group-specific mean and covariance ---
  # use distinct matrix name and labels to avoid collision
  expMean_ACEc <- mxMatrix(type="Full", nrow=1, ncol=2, free=TRUE, values=0,
                           labels=c("par_mean","par_mean"), name="expMean")

  # Build expCov inside each group so A/C/D/E are in scope
  mzCov_ACEc <- mxAlgebra(
    rbind(
      cbind(A + C + D + E, A + C + D),
      cbind(A + C + D,     A + C + D + E)
    ),
    name="expCov"
  )
  dzCov_ACEc <- mxAlgebra(
    rbind(
      cbind(A + C + D + E, 0.5 %x% A + C + 0.25 %x% D),
      cbind(0.5 %x% A + C + 0.25 %x% D, A + C + D + E)
    ),
    name="expCov"
  )

  expMZ_ACEc <- mxExpectationNormal(covariance="expCov", means="expMean", dimnames=selVars)
  expDZ_ACEc <- mxExpectationNormal(covariance="expCov", means="expMean", dimnames=selVars)
  fitML      <- mxFitFunctionML()

  ## --- assemble group models ---
  modelMZ_ACEc <- mxModel("MZ",
    a_ACEc, c_ACEc, d_ACEc, e_ACEc,
    A_ACEc, C_ACEc, D_ACEc, E_ACEc,
    a2_ACEc, c2_ACEc, d2_ACEc, e2_ACEc,
    expMean_ACEc, mzCov_ACEc, datMZ, expMZ_ACEc, fitML
  )
  modelDZ_ACEc <- mxModel("DZ",
    a_ACEc, c_ACEc, d_ACEc, e_ACEc,
    A_ACEc, C_ACEc, D_ACEc, E_ACEc,
    a2_ACEc, c2_ACEc, d2_ACEc, e2_ACEc,
    expMean_ACEc, dzCov_ACEc, datDZ, expDZ_ACEc, fitML
  )

  ## --- top-level standardization & multigroup fit ---
  Tot_ACEc  <- mxAlgebra(MZ.a2 + MZ.c2 + MZ.d2 + MZ.e2, name="Tot")
  Prop_ACEc <- mxAlgebra(
    cbind(MZ.a2/Tot, MZ.c2/Tot, MZ.d2/Tot, MZ.e2/Tot),
    name="Prop",
    dimnames=list(NULL, c("a2","c2","d2","e2"))
  )
  multiTwin <- mxFitFunctionMultigroup(c("MZ","DZ"))

  fitACEc <- mxModel("ACEc", modelMZ_ACEc, modelDZ_ACEc, Tot_ACEc, Prop_ACEc, multiTwin)
  fitACEc <- safe_run(fitACEc, "ACEc")

  if (!is.null(fitACEc)) {
    checkVarianceComponents(fitACEc)
    if (logResults) results_log <- logModelResults(fitACEc, results_log)
    fits$ACEc <- fitACEc
  }
}

# ==============================================================================
# GENERATE AND FIT LOWER BOUND VARIANTS
# ==============================================================================

cat("\n=== GENERATING LOWER BOUND VARIANTS ===\n")

# Store base models to check for lb variants
base_models <- c("AE", "CE", "ACE", "DE", "ADE")  # E model has no free variance components

for (model_name in base_models) {
  if (model_name %in% names(fits)) {
    base_fit <- fits[[model_name]]
    base_specs <- model_specs[[model_name]]

    # Generate lb variants for this model
    lb_variants <- generate_lb_variants(base_fit, base_specs, lb_threshold, lb_start, lb_upper, lb_step)

    if (length(lb_variants) > 0) {
      cat("Fitting", length(lb_variants), "lb variants for", model_name, "model...\n")

      # Fit each variant
      for (variant_name in names(lb_variants)) {
        variant_model <- lb_variants[[variant_name]]
        variant_fit <- safe_run(variant_model, variant_name)

        if (!is.null(variant_fit)) {
          checkVarianceComponents(variant_fit, verbose = FALSE)  # Less verbose for variants
          if (exists("logResults") && isTRUE(logResults)) {
            results_log <- logModelResults(variant_fit, results_log)
          }
          fits[[variant_name]] <- variant_fit
        }
      }
    }
  }
}

# ==============================================================================
# SAVE RESULTS LOG
# ==============================================================================

if(logResults) {
  log_filename <- paste0("univariate_fit_results_trait", target_trait, "_",
                        format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  write.csv(results_log, log_filename, row.names = FALSE)
  cat("Results logged to:", log_filename, "\n")
}

# ==============================================================================
# HELPER FUNCTIONS AND REPORTING
# ==============================================================================

# Helper function to report model fits with both raw and standardized components
report_fits <- function(model_list) {
  cat("\n=== MODEL FIT SUMMARY ===\n")
  for(i in seq_along(model_list)) {
    model <- model_list[[i]]
    name <- names(model_list)[i]

    cat(sprintf("%-10s: -2LL = %8.2f, df = %3d, AIC = %8.2f\n",
                name,
                model$output$Minus2LogLikelihood,
                model$output$degreesOfFreedom,
                model$output$AIC))

    # Print raw variance components (namespaced)
    rv <- model$output$algebras
    if(all(c("MZ.a2", "MZ.c2", "MZ.d2", "MZ.e2") %in% names(rv))) {
      raw_vals <- c(rv$`MZ.a2`, rv$`MZ.c2`, rv$`MZ.d2`, rv$`MZ.e2`)
      cat(sprintf("           Raw:  a2=%6.4f, c2=%6.4f, d2=%6.4f, e2=%6.4f (tot=%6.4f)\n",
                  raw_vals[1], raw_vals[2], raw_vals[3], raw_vals[4], sum(raw_vals)))
    }

    # Print standardized proportions
    if(!is.null(rv$Prop)) {
      pp <- as.numeric(rv$Prop)
      cat(sprintf("           Prop: a2=%.3f c2=%.3f d2=%.3f e2=%.3f (sum=%.3f)\n",
                  pp[1], pp[2], pp[3], pp[4], sum(pp)))
    }
  }
}

# Report all model fits
## ---------- Robust helpers for AIC weights & model averaging ----------

# Is a fitted object usable?
is_ok_fit <- function(m) {
  inherits(m, "MxModel") &&
    !is.null(m$output) &&
    is.finite(m$output$Minus2LogLikelihood)
}

# Number of free parameters (k)
n_params <- function(m) {
  if (!is_ok_fit(m)) return(NA_integer_)
  k <- length(m$output$estimate)
  if (is.null(k) || !is.finite(k)) k <- NA_integer_
  k
}

# -2 log-likelihood
get_m2ll <- function(m) if (is_ok_fit(m)) as.numeric(m$output$Minus2LogLikelihood) else NA_real_

# AIC = -2lnL + 2k
get_aic <- function(m) {
  m2 <- get_m2ll(m); k <- n_params(m)
  if (is.finite(m2) && is.finite(k)) m2 + 2*k else NA_real_
}

# Vector of AIC weights (named like fits)
aic_weights <- function(fits) {
  aics <- sapply(fits, get_aic)
  ok   <- is.finite(aics)
  if (!any(ok)) {
    warning("aic_weights: no valid AICs found in 'fits'. Returning equal weights.")
    w <- rep(1/length(aics), length(aics))
  } else {
    dAIC <- aics - min(aics[ok])
    w    <- exp(-0.5*dAIC)
    w[!ok] <- 0
    w    <- w / sum(w)
  }
  names(w) <- names(fits)
  w
}

# Does a model actually estimate a component? (checks 'MZ' submodel free flags)
comp_is_free <- function(m, comp = c("a2","c2","d2","e2")) {
  comp <- match.arg(comp)
  nm <- switch(comp, a2="a", c2="c", d2="d", e2="e")
  if (!inherits(m, "MxModel") || is.null(m$MZ) || is.null(m$MZ[[nm]])) return(FALSE)
  fr <- m$MZ[[nm]]$free
  isTRUE(any(fr))
}

# Extract raw variance component from output algebras (prefer what's already computed)
extract_comp <- function(m, comp = c("a2","c2","d2","e2")) {
  comp <- match.arg(comp)
  if (!is_ok_fit(m)) return(NA_real_)
  alg <- switch(comp, a2="MZ.a2", c2="MZ.c2", d2="MZ.d2", e2="MZ.e2")
  if (!is.null(m$output$algebras) && alg %in% names(m$output$algebras)) {
    as.numeric(m$output$algebras[[alg]])
  } else {
    # Fallback (rare): try to evaluate by name
    val <- tryCatch(as.numeric(mxEvalByName(alg, m, compute=TRUE)),
                    error = function(e) NA_real_)
    val
  }
}

# Conditional average (weights renormalized over models that estimate the component)
conditional_avg <- function(fits, comp) {
  w  <- aic_weights(fits)
  ok <- sapply(fits, comp_is_free, comp = comp)
  if (!any(ok)) return(NA_real_)
  ww <- w; ww[!ok] <- 0; ww <- ww / sum(ww)
  vals <- sapply(fits, extract_comp, comp = comp)
  sum(ww * vals, na.rm = TRUE)
}

# Unconditional average (global weights; models without the comp contribute 0)
unconditional_avg <- function(fits, comp) {
  w    <- aic_weights(fits)
  vals <- sapply(fits, extract_comp, comp = comp)
  vals[!sapply(fits, comp_is_free, comp = comp)] <- 0
  sum(w * vals, na.rm = TRUE)
}

# Average *proportions* across models (choose "conditional" or "unconditional")
avg_props <- function(fits, mode = c("conditional","unconditional")) {
  mode <- match.arg(mode)
  fun  <- if (mode == "conditional") conditional_avg else unconditional_avg
  raw  <- c(A = fun(fits, "a2"),
            C = fun(fits, "c2"),
            D = fun(fits, "d2"),
            E = fun(fits, "e2"))
  if (all(!is.finite(raw))) return(setNames(rep(NA_real_,4), c("A","C","D","E")))
  raw[is.na(raw)] <- 0
  s <- sum(raw)
  if (!is.finite(s) || s <= 0) return(setNames(rep(NA_real_,4), c("A","C","D","E")))
  raw / s
}

# Pretty AIC table
print_aic_table <- function(fits) {
  aics <- sapply(fits, get_aic)
  ok   <- is.finite(aics)
  dlt  <- if (any(ok)) aics - min(aics[ok]) else rep(NA_real_, length(aics))
  wt   <- aic_weights(fits)
  out  <- data.frame(model   = names(fits),
                     AIC     = aics,
                     delta   = dlt,
                     weight  = wt,
                     hasAIC  = ok,
                     row.names = NULL, check.names = FALSE)
  out[order(out$AIC, na.last = TRUE), ]
}
## ---------------------------------------------------------------------

# After you successfully fit all models (including lb variants):
cat("\n=== MODEL COMPARISON TABLE ===\n")
print(print_aic_table(fits))

# Raw component averages (variance scale)
A_cond <- conditional_avg(fits, "a2")
C_cond <- conditional_avg(fits, "c2")
D_cond <- conditional_avg(fits, "d2")
E_cond <- conditional_avg(fits, "e2")

cat("\n=== MODEL AVERAGED RESULTS ===\n")
cat("Conditional averages (variance scale):\n")
cat("A:", round(A_cond, 4), "C:", round(C_cond, 4), "D:", round(D_cond, 4), "E:", round(E_cond, 4), "\n")

# Proportion scale
cond_props <- avg_props(fits, "conditional")
uncond_props <- avg_props(fits, "unconditional")

cat("Conditional proportions:\n")
cat("A:", round(cond_props["A"], 3), "C:", round(cond_props["C"], 3),
    "D:", round(cond_props["D"], 3), "E:", round(cond_props["E"], 3), "\n")

cat("Unconditional proportions:\n")
cat("A:", round(uncond_props["A"], 3), "C:", round(uncond_props["C"], 3),
    "D:", round(uncond_props["D"], 3), "E:", round(uncond_props["E"], 3), "\n")

# Define parameter labels for model averaging
colPra <- c("a11", "c11", "d11", "e11")
# Test mxModelAverage(), with free parameters only (which is faster), for all permutations of arguments 'include', 'SE', 'refAsBlock', & 'type'
mma <- mxModelAverage(reference=c(colPra,'Prop'), models=fits, include='onlyFree', SE=NULL, refAsBlock=FALSE, type='AIC')
