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
startVals    <- c(A = 0.3, C = 0.3, D = 0.3, E = 0.3)  # Starting values for variance components
includeOSDZ  <- TRUE      # Include opposite-sex DZ twins in DZ pool

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
cat("Starting values - A:", startVals[["A"]], "C:", startVals[["C"]],
    "D:", startVals[["D"]], "E:", startVals[["E"]], "\n")
cat("Include OSDZ:", includeOSDZ, "\n")

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

# Calculate observed correlations for reference
mz_cor <- cor(dataMZ[, selVars], use = "complete.obs")[1,2]
dz_cor <- cor(dataDZ[, selVars], use = "complete.obs")[1,2]

cat("\n=== OBSERVED TWIN CORRELATIONS ===\n")
cat("MZ correlation:", round(mz_cor, 3), "\n")
cat("DZ correlation:", round(dz_cor, 3), "\n")

# Determine best model based on correlations
if (mz_cor > 2 * dz_cor) {
  cat("MZ correlation > 2 * DZ correlation suggests dominance effects (ADE model)\n")
  correlation_based_model <- "ADE"
} else {
  cat("MZ correlation < 2 * DZ correlation suggests shared environment effects (ACE model)\n")
  correlation_based_model <- "ACE"
}

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

# Expectations and fit functions
expMZ <- mxExpectationNormal(covariance = "mzCov", means = "Mean", dimnames = selVars)
expDZ <- mxExpectationNormal(covariance = "dzCov", means = "Mean", dimnames = selVars)
fitML <- mxFitFunctionML()

# Multigroup fit function
multiTwin <- mxFitFunctionMultigroup(c("MZ","DZ"))

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Function to translate model name to "free flags" for both builders
freeSpec <- function(name) {
  switch(name,
         "E"   = c(A=FALSE, C=FALSE, D=FALSE, E=TRUE),
         "AE"  = c(A=TRUE,  C=FALSE, D=FALSE, E=TRUE),
         "CE"  = c(A=FALSE, C=TRUE,  D=FALSE, E=TRUE),
         "ACE" = c(A=TRUE,  C=TRUE,  D=FALSE, E=TRUE),
         "DE"  = c(A=FALSE, C=FALSE, D=TRUE,  E=TRUE),
         "ADE" = c(A=TRUE,  C=FALSE, D=TRUE,  E=TRUE),
         stop("Unknown model: ", name)
  )
}

# Path-coefficient builder (traditional approach)
buildTwinModel_path <- function(modelName, a_free = FALSE, c_free = FALSE, d_free = FALSE, e_free = TRUE) {

  # Set starting values
  a_val <- if(a_free) startVals[["A"]] else 0
  c_val <- if(c_free) startVals[["C"]] else 0
  d_val <- if(d_free) startVals[["D"]] else 0
  e_val <- if(e_free) max(startVals[["E"]], 0.1) else 0

  # Set lower bounds
  a_lb <- if(a_free) 0 else NA
  c_lb <- if(c_free) 0 else NA
  d_lb <- if(d_free) 0 else NA
  e_lb <- if(e_free) 1e-6 else NA

  # Path coefficient matrices (labels MUST NOT equal object names)
  a_mat <- mxMatrix(type="Lower", nrow=1, ncol=1, free=a_free,
                    values = a_val, labels = "a11", lbound = a_lb, name="a")

  c_mat <- mxMatrix(type="Lower", nrow=1, ncol=1, free=c_free,
                    values = c_val, labels = "c11", lbound = c_lb, name="c")

  d_mat <- mxMatrix(type="Lower", nrow=1, ncol=1, free=d_free,
                    values = d_val, labels = "d11", lbound = d_lb, name="d")

  e_mat <- mxMatrix(type="Lower", nrow=1, ncol=1, free=e_free,
                    values = e_val, labels = "e11", lbound = e_lb, name="e")

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

  # Expected covariance algebras
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

# Direct-variance builder (new approach)
buildTwinModel_var <- function(modelName, A_free = FALSE, C_free = FALSE, D_free = FALSE, E_free = TRUE) {

  # Set starting values
  A_val <- if(A_free) startVals[["A"]] else 0
  C_val <- if(C_free) startVals[["C"]] else 0
  D_val <- if(D_free) startVals[["D"]] else 0
  E_val <- if(E_free) max(startVals[["E"]], 0.1) else 0

  # Set lower bounds
  A_lb <- if(A_free) 0 else NA
  C_lb <- if(C_free) 0 else NA
  D_lb <- if(D_free) 0 else NA
  E_lb <- if(E_free) 1e-6 else NA

  # Direct variance matrices
  A <- mxMatrix(type="Symm", nrow=1, ncol=1, free=A_free,
                values = A_val, labels = "A11", lbound = A_lb, name="A")

  C <- mxMatrix(type="Symm", nrow=1, ncol=1, free=C_free,
                values = C_val, labels = "C11", lbound = C_lb, name="C")

  D <- mxMatrix(type="Symm", nrow=1, ncol=1, free=D_free,
                values = D_val, labels = "D11", lbound = D_lb, name="D")

  E <- mxMatrix(type="Symm", nrow=1, ncol=1, free=E_free,
                values = E_val, labels = "E11", lbound = E_lb, name="E")

  # Raw component scalars (for logging/model averaging)
  a2_alg <- mxAlgebra(A[1,1], name="a2")
  c2_alg <- mxAlgebra(C[1,1], name="c2")
  d2_alg <- mxAlgebra(D[1,1], name="d2")
  e2_alg <- mxAlgebra(E[1,1], name="e2")

  # Expected covariance algebras
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

  # Submodels
  modelMZ <- mxModel("MZ", A, C, D, E,
                     a2_alg, c2_alg, d2_alg, e2_alg,
                     Mean, mzCov, datMZ, expMZ, fitML)

  modelDZ <- mxModel("DZ", A, C, D, E,
                     a2_alg, c2_alg, d2_alg, e2_alg,
                     Mean, dzCov, datDZ, expDZ, fitML)

  # Top-level: standardized proportions (use namespaced MZ.*)
  Tot  <- mxAlgebra(MZ.a2 + MZ.c2 + MZ.d2 + MZ.e2, name="Tot")
  Prop <- mxAlgebra(cbind(MZ.a2/Tot, MZ.c2/Tot, MZ.d2/Tot, MZ.e2/Tot),
                    name="Prop",
                    dimnames=list(NULL, c("a2","c2","d2","e2")))

  mxModel(modelName, modelMZ, modelDZ, Tot, Prop, multiTwin)
}

# Safe runner + helpers
runner <- if (exists("useTryHard") && isTRUE(useTryHard)) {
  # let user set tryHardTries elsewhere; default to 10
  function(m) mxTryHard(m, extraTries = if (exists("tryHardTries")) tryHardTries else 10)
} else {
  function(m) mxRun(m)
}

is_ok_fit <- function(m) {
  inherits(m, "MxModel") &&
    !is.null(m$output) &&
    is.finite(m$output$Minus2LogLikelihood)
}

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

# Helper functions for model comparison and averaging
n_params <- function(m) {
  if (!is_ok_fit(m)) return(NA_integer_)
  k <- length(m$output$estimate)
  if (is.null(k) || !is.finite(k)) k <- NA_integer_
  k
}

get_m2ll <- function(m) if (is_ok_fit(m)) as.numeric(m$output$Minus2LogLikelihood) else NA_real_

get_aic <- function(m) {
  m2 <- get_m2ll(m); k <- n_params(m)
  if (is.finite(m2) && is.finite(k)) m2 + 2*k else NA_real_
}

# Robust component extraction that works for both model types
extract_comp <- function(m, comp = c("a2","c2","d2","e2")) {
  comp <- match.arg(comp)
  if (!is_ok_fit(m)) return(NA_real_)

  # Try to get from output algebras first
  alg <- paste0("MZ.", comp)
  if (!is.null(m$output$algebras) && alg %in% names(m$output$algebras)) {
    return(as.numeric(m$output$algebras[[alg]]))
  }

  # Fallback: try to evaluate by name
  val <- tryCatch(as.numeric(mxEvalByName(alg, m, compute=TRUE)),
                  error = function(e) NA_real_)
  return(val)
}

# Check if a component is free in a model (works for both model types)
comp_is_free <- function(m, comp = c("a2","c2","d2","e2")) {
  comp <- match.arg(comp)

  # Map component to matrix name
  path_matrix <- switch(comp, a2="a", c2="c", d2="d", e2="e")
  var_matrix <- switch(comp, a2="A", c2="C", d2="D", e2="E")

  # Check if model has MZ submodel
  if (!inherits(m, "MxModel") || is.null(m$MZ)) return(FALSE)

  # Try path approach first
  if (!is.null(m$MZ[[path_matrix]])) {
    return(isTRUE(any(m$MZ[[path_matrix]]$free)))
  }

  # Try variance approach
  if (!is.null(m$MZ[[var_matrix]])) {
    return(isTRUE(any(m$MZ[[var_matrix]]$free)))
  }

  # Component not found or not free
  return(FALSE)
}

# Function to create AIC table
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

# Function to create per-model parameter table with SEs
create_parameter_table <- function(fits, model_type) {
  # Initialize empty data frame
  result_df <- data.frame()

  # Check if fits is empty
  if (length(fits) == 0) {
    return(result_df)  # Return empty data frame
  }

  # Process each model
  for (model_name in names(fits)) {
    fit <- fits[[model_name]]

    # Skip if fit is not valid (check for proper output structure)
    if (!inherits(fit, "MxModel") || is.null(fit$output) ||
        is.null(fit$output$Minus2LogLikelihood) ||
        !is.finite(fit$output$Minus2LogLikelihood)) {
      cat("Warning: Model", model_name, "does not have valid fitted results. Skipping.\n")
      next
    }

    # Extract basic fit statistics
    m2ll <- as.numeric(fit$output$Minus2LogLikelihood)

    # Calculate degrees of freedom
    n_params <- length(fit$output$estimate)
    total_obs <- sum(sapply(fit$output$data, function(d) d$numObs))
    n_vars <- 2  # Twin 1 and Twin 2
    df <- total_obs * n_vars - n_params

    # Calculate AIC
    aic <- m2ll + 2 * n_params

    # Extract status code
    status <- if (!is.null(fit$output$status$code)) fit$output$status$code else NA_integer_

    # Extract variance components from algebras
    a2 <- if (!is.null(fit$output$algebras$MZ.a2)) as.numeric(fit$output$algebras$MZ.a2) else NA_real_
    c2 <- if (!is.null(fit$output$algebras$MZ.c2)) as.numeric(fit$output$algebras$MZ.c2) else NA_real_
    d2 <- if (!is.null(fit$output$algebras$MZ.d2)) as.numeric(fit$output$algebras$MZ.d2) else NA_real_
    e2 <- if (!is.null(fit$output$algebras$MZ.e2)) as.numeric(fit$output$algebras$MZ.e2) else NA_real_

    # Calculate standard errors with improved error handling
    a2_se <- tryCatch({
      if (comp_is_free(fit, "a2")) {
        mxSE(MZ.a2, fit)
      } else {
        0  # Fixed parameter has SE = 0
      }
    }, error = function(e) NA_real_)

    c2_se <- tryCatch({
      if (comp_is_free(fit, "c2")) {
        mxSE(MZ.c2, fit)
      } else {
        0  # Fixed parameter has SE = 0
      }
    }, error = function(e) NA_real_)

    d2_se <- tryCatch({
      if (comp_is_free(fit, "d2")) {
        mxSE(MZ.d2, fit)
      } else {
        0  # Fixed parameter has SE = 0
      }
    }, error = function(e) NA_real_)

    e2_se <- tryCatch({
      if (comp_is_free(fit, "e2")) {
        mxSE(MZ.e2, fit)
      } else {
        0  # Fixed parameter has SE = 0
      }
    }, error = function(e) NA_real_)

    # Calculate proportions
    total <- sum(c(a2, c2, d2, e2), na.rm = TRUE)
    a2_prop <- if (is.finite(a2) && is.finite(total) && total > 0) a2/total else NA_real_
    c2_prop <- if (is.finite(c2) && is.finite(total) && total > 0) c2/total else NA_real_
    d2_prop <- if (is.finite(d2) && is.finite(total) && total > 0) d2/total else NA_real_
    e2_prop <- if (is.finite(e2) && is.finite(total) && total > 0) e2/total else NA_real_

    # Create row for this model
    model_row <- data.frame(
      Model = model_name,
      Type = model_type,
      minus2LL = m2ll,
      AIC = aic,
      df = df,
      status = status,
      A_var = a2,
      A_SE = a2_se,
      C_var = c2,
      C_SE = c2_se,
      D_var = d2,
      D_SE = d2_se,
      E_var = e2,
      E_SE = e2_se,
      A_prop = a2_prop,
      C_prop = c2_prop,
      D_prop = d2_prop,
      E_prop = e2_prop,
      stringsAsFactors = FALSE
    )

    # Add to result data frame
    result_df <- rbind(result_df, model_row)
  }

  return(result_df)
}

# Function to find best model based on AIC
find_best_model <- function(fits) {
  if (length(fits) == 0) {
    return(NULL)
  }

  aics <- sapply(fits, get_aic)
  ok <- is.finite(aics)

  if (!any(ok)) {
    return(NULL)
  }

  best_model_name <- names(fits)[which.min(aics)]
  return(list(name = best_model_name, model = fits[[best_model_name]]))
}

# ==============================================================================
# BUILD AND FIT MODELS
# ==============================================================================

cat("\n=== FITTING MODELS (PATH COEFFICIENT APPROACH) ===\n")
fits_path <- list()

# Build and fit the six standard models using path approach
for (model_name in c("E", "AE", "CE", "ACE", "DE", "ADE")) {
  cat("Fitting", model_name, "model (path approach)...\n")
  sp <- freeSpec(model_name)
  model <- buildTwinModel_path(model_name,
                              a_free=sp["A"],
                              c_free=sp["C"],
                              d_free=sp["D"],
                              e_free=sp["E"])
  fit <- safe_run(model, model_name)
  if (is_ok_fit(fit)) {
    fits_path[[model_name]] <- fit
  }
}

cat("\n=== FITTING MODELS (DIRECT VARIANCE APPROACH) ===\n")
fits_var <- list()

# Build and fit the six standard models using variance approach
for (model_name in c("E", "AE", "CE", "ACE", "DE", "ADE")) {
  cat("Fitting", model_name, "model (variance approach)...\n")
  sp <- freeSpec(model_name)
  model <- buildTwinModel_var(model_name,
                             A_free=sp["A"],
                             C_free=sp["C"],
                             D_free=sp["D"],
                             E_free=sp["E"])
  fit <- safe_run(model, model_name)
  if (is_ok_fit(fit)) {
    fits_var[[model_name]] <- fit
  }
}

# ==============================================================================
# MODEL COMPARISON AND AVERAGING
# ==============================================================================

# Define parameter labels for model averaging
colPra_path <- c("a11", "c11", "d11", "e11")
colPra_var <- c("A11", "C11", "D11", "E11")

# Print AIC tables
cat("\n=== MODEL COMPARISON TABLE (PATH APPROACH) ===\n")
if (length(fits_path) > 0) {
  aic_table_path <- print_aic_table(fits_path)
  print(aic_table_path)
} else {
  cat("No valid path models to compare.\n")
}

cat("\n=== MODEL COMPARISON TABLE (VARIANCE APPROACH) ===\n")
if (length(fits_var) > 0) {
  aic_table_var <- print_aic_table(fits_var)
  print(aic_table_var)
} else {
  cat("No valid variance models to compare.\n")
}

# Find best models based on AIC
best_path_model <- find_best_model(fits_path)
best_var_model <- find_best_model(fits_var)

# Print best models
cat("\n=== BEST MODELS BASED ON AIC ===\n")
if (!is.null(best_path_model)) {
  cat("Best path model:", best_path_model$name, "\n")
} else {
  cat("No valid path model found.\n")
}

if (!is.null(best_var_model)) {
  cat("Best variance model:", best_var_model$name, "\n")
} else {
  cat("No valid variance model found.\n")
}

cat("\n=== BEST MODEL BASED ON CORRELATIONS ===\n")
cat("Correlation-based model selection:", correlation_based_model, "\n")

# Run model averaging for path approach
cat("\n=== MODEL AVERAGING (PATH APPROACH) ===\n")
mma_path <- NULL
if (length(fits_path) > 0) {
  mma_path <- tryCatch(
    mxModelAverage(reference=c(colPra_path,"Prop"),
                  models=fits_path,
                  include='onlyFree',
                  SE=NULL,
                  refAsBlock=FALSE,
                  type='AIC'),
    error = function(e) {
      cat("Error in model averaging (path approach):", e$message, "\n")
      NULL
    }
  )
} else {
  cat("No valid path models for model averaging.\n")
}

# Run model averaging for variance approach
cat("\n=== MODEL AVERAGING (VARIANCE APPROACH) ===\n")
mma_var <- NULL
if (length(fits_var) > 0) {
  mma_var <- tryCatch(
    mxModelAverage(reference=c(colPra_var,"Prop"),
                  models=fits_var,
                  include='onlyFree',
                  SE=NULL,
                  refAsBlock=FALSE,
                  type='AIC'),
    error = function(e) {
      cat("Error in model averaging (variance approach):", e$message, "\n")
      NULL
    }
  )
} else {
  cat("No valid variance models for model averaging.\n")
}

# ==============================================================================
# COMBINE RESULTS INTO A SINGLE OUTPUT
# ==============================================================================

# Create timestamp for filenames
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create combined results dataframe
combined_results <- data.frame(
  Parameter = c("a11", "c11", "d11", "e11", "prop1", "prop2", "prop3", "prop4"),
  stringsAsFactors = FALSE
)

# Add multimodel path results
if (!is.null(mma_path) && !is.null(mma_path$"Model-Average Estimates")) {
  # Extract model-averaged estimates and convert to data frame
  path_estimates <- as.data.frame(mma_path$"Model-Average Estimates")
  path_estimates$name <- rownames(path_estimates)

  # Extract path coefficients
  path_a11 <- path_estimates$Estimate[path_estimates$name == "a11"]
  path_c11 <- path_estimates$Estimate[path_estimates$name == "c11"]
  path_d11 <- path_estimates$Estimate[path_estimates$name == "d11"]
  path_e11 <- path_estimates$Estimate[path_estimates$name == "e11"]

  # Extract standard errors
  path_a11_se <- path_estimates$SE[path_estimates$name == "a11"]
  path_c11_se <- path_estimates$SE[path_estimates$name == "c11"]
  path_d11_se <- path_estimates$SE[path_estimates$name == "d11"]
  path_e11_se <- path_estimates$SE[path_estimates$name == "e11"]

  # Extract proportions
  path_prop <- path_estimates[grep("Prop", path_estimates$name), ]
  path_prop1 <- path_prop$Estimate[1]
  path_prop2 <- path_prop$Estimate[2]
  path_prop3 <- path_prop$Estimate[3]
  path_prop4 <- path_prop$Estimate[4]

  # Extract proportion SEs
  path_prop1_se <- path_prop$SE[1]
  path_prop2_se <- path_prop$SE[2]
  path_prop3_se <- path_prop$SE[3]
  path_prop4_se <- path_prop$SE[4]

  # Add to combined results
  combined_results$MM_Path_Est <- c(path_a11, path_c11, path_d11, path_e11,
                                   path_prop1, path_prop2, path_prop3, path_prop4)
  combined_results$MM_Path_SE <- c(path_a11_se, path_c11_se, path_d11_se, path_e11_se,
                                  path_prop1_se, path_prop2_se, path_prop3_se, path_prop4_se)
} else {
  combined_results$MM_Path_Est <- NA
  combined_results$MM_Path_SE <- NA
}

# Add multimodel variance results
if (!is.null(mma_var) && !is.null(mma_var$"Model-Average Estimates")) {
  # Extract model-averaged estimates and convert to data frame
  var_estimates <- as.data.frame(mma_var$"Model-Average Estimates")
  var_estimates$name <- rownames(var_estimates)

  # Extract variance components
  var_A11 <- var_estimates$Estimate[var_estimates$name == "A11"]
  var_C11 <- var_estimates$Estimate[var_estimates$name == "C11"]
  var_D11 <- var_estimates$Estimate[var_estimates$name == "D11"]
  var_E11 <- var_estimates$Estimate[var_estimates$name == "E11"]

  # Extract standard errors
  var_A11_se <- var_estimates$SE[var_estimates$name == "A11"]
  var_C11_se <- var_estimates$SE[var_estimates$name == "C11"]
  var_D11_se <- var_estimates$SE[var_estimates$name == "D11"]
  var_E11_se <- var_estimates$SE[var_estimates$name == "E11"]

  # Extract proportions
  var_prop <- var_estimates[grep("Prop", var_estimates$name), ]
  var_prop1 <- var_prop$Estimate[1]
  var_prop2 <- var_prop$Estimate[2]
  var_prop3 <- var_prop$Estimate[3]
  var_prop4 <- var_prop$Estimate[4]

  # Extract proportion SEs
  var_prop1_se <- var_prop$SE[1]
  var_prop2_se <- var_prop$SE[2]
  var_prop3_se <- var_prop$SE[3]
  var_prop4_se <- var_prop$SE[4]

  # Add to combined results
  combined_results$MM_Var_Est <- c(var_A11, var_C11, var_D11, var_E11,
                                  var_prop1, var_prop2, var_prop3, var_prop4)
  combined_results$MM_Var_SE <- c(var_A11_se, var_C11_se, var_D11_se, var_E11_se,
                                 var_prop1_se, var_prop2_se, var_prop3_se, var_prop4_se)
} else {
  combined_results$MM_Var_Est <- NA
  combined_results$MM_Var_SE <- NA
}

# Add best single model path results
if (!is.null(best_path_model)) {
  best_path <- best_path_model$model

  # Extract path coefficients
  best_path_a11 <- if (comp_is_free(best_path, "a2")) sqrt(extract_comp(best_path, "a2")) else 0
  best_path_c11 <- if (comp_is_free(best_path, "c2")) sqrt(extract_comp(best_path, "c2")) else 0
  best_path_d11 <- if (comp_is_free(best_path, "d2")) sqrt(extract_comp(best_path, "d2")) else 0
  best_path_e11 <- if (comp_is_free(best_path, "e2")) sqrt(extract_comp(best_path, "e2")) else 0

  # Calculate standard errors
  best_path_a11_se <- tryCatch({
    if (comp_is_free(best_path, "a2")) {
      sqrt(mxSE(MZ.a2, best_path) / (4 * extract_comp(best_path, "a2")))
    } else {
      0
    }
  }, error = function(e) NA_real_)

  best_path_c11_se <- tryCatch({
    if (comp_is_free(best_path, "c2")) {
      sqrt(mxSE(MZ.c2, best_path) / (4 * extract_comp(best_path, "c2")))
    } else {
      0
    }
  }, error = function(e) NA_real_)

  best_path_d11_se <- tryCatch({
    if (comp_is_free(best_path, "d2")) {
      sqrt(mxSE(MZ.d2, best_path) / (4 * extract_comp(best_path, "d2")))
    } else {
      0
    }
  }, error = function(e) NA_real_)

  best_path_e11_se <- tryCatch({
    if (comp_is_free(best_path, "e2")) {
      sqrt(mxSE(MZ.e2, best_path) / (4 * extract_comp(best_path, "e2")))
    } else {
      0
    }
  }, error = function(e) NA_real_)

  # Extract proportions
  best_path_a2 <- extract_comp(best_path, "a2")
  best_path_c2 <- extract_comp(best_path, "c2")
  best_path_d2 <- extract_comp(best_path, "d2")
  best_path_e2 <- extract_comp(best_path, "e2")

  best_path_total <- sum(c(best_path_a2, best_path_c2, best_path_d2, best_path_e2), na.rm = TRUE)

  best_path_prop1 <- best_path_a2 / best_path_total
  best_path_prop2 <- best_path_c2 / best_path_total
  best_path_prop3 <- best_path_d2 / best_path_total
  best_path_prop4 <- best_path_e2 / best_path_total

  # Calculate proportion SEs (simplified approximation)
  best_path_prop1_se <- NA
  best_path_prop2_se <- NA
  best_path_prop3_se <- NA
  best_path_prop4_se <- NA

  # Add to combined results
  combined_results$SM_Path_Est <- c(best_path_a11, best_path_c11, best_path_d11, best_path_e11,
                                   best_path_prop1, best_path_prop2, best_path_prop3, best_path_prop4)
  combined_results$SM_Path_SE <- c(best_path_a11_se, best_path_c11_se, best_path_d11_se, best_path_e11_se,
                                  best_path_prop1_se, best_path_prop2_se, best_path_prop3_se, best_path_prop4_se)
} else {
  combined_results$SM_Path_Est <- NA
  combined_results$SM_Path_SE <- NA
}

# Add best single model variance results
if (!is.null(best_var_model)) {
  best_var <- best_var_model$model

  # Extract variance components
  best_var_A11 <- extract_comp(best_var, "a2")
  best_var_C11 <- extract_comp(best_var, "c2")
  best_var_D11 <- extract_comp(best_var, "d2")
  best_var_E11 <- extract_comp(best_var, "e2")

  # Calculate standard errors
  best_var_A11_se <- tryCatch({
    if (comp_is_free(best_var, "a2")) {
      mxSE(MZ.a2, best_var)
    } else {
      0
    }
  }, error = function(e) NA_real_)

  best_var_C11_se <- tryCatch({
    if (comp_is_free(best_var, "c2")) {
      mxSE(MZ.c2, best_var)
    } else {
      0
    }
  }, error = function(e) NA_real_)

  best_var_D11_se <- tryCatch({
    if (comp_is_free(best_var, "d2")) {
      mxSE(MZ.d2, best_var)
    } else {
      0
    }
  }, error = function(e) NA_real_)

  best_var_E11_se <- tryCatch({
    if (comp_is_free(best_var, "e2")) {
      mxSE(MZ.e2, best_var)
    } else {
      0
    }
  }, error = function(e) NA_real_)

  # Extract proportions
  best_var_total <- sum(c(best_var_A11, best_var_C11, best_var_D11, best_var_E11), na.rm = TRUE)

  best_var_prop1 <- best_var_A11 / best_var_total
  best_var_prop2 <- best_var_C11 / best_var_total
  best_var_prop3 <- best_var_D11 / best_var_total
  best_var_prop4 <- best_var_E11 / best_var_total

  # Calculate proportion SEs (simplified approximation)
  best_var_prop1_se <- NA
  best_var_prop2_se <- NA
  best_var_prop3_se <- NA
  best_var_prop4_se <- NA

  # Add to combined results
  combined_results$SM_Var_Est <- c(best_var_A11, best_var_C11, best_var_D11, best_var_E11,
                                  best_var_prop1, best_var_prop2, best_var_prop3, best_var_prop4)
  combined_results$SM_Var_SE <- c(best_var_A11_se, best_var_C11_se, best_var_D11_se, best_var_E11_se,
                                 best_var_prop1_se, best_var_prop2_se, best_var_prop3_se, best_var_prop4_se)
} else {
  combined_results$SM_Var_Est <- NA
  combined_results$SM_Var_SE <- NA
}

# Add correlation-based model results
if (correlation_based_model == "ACE") {
  # Get ACE model results
  if ("ACE" %in% names(fits_path)) {
    corr_model_path <- fits_path[["ACE"]]

    # Extract path coefficients
    corr_path_a11 <- if (comp_is_free(corr_model_path, "a2")) sqrt(extract_comp(corr_model_path, "a2")) else 0
    corr_path_c11 <- if (comp_is_free(corr_model_path, "c2")) sqrt(extract_comp(corr_model_path, "c2")) else 0
    corr_path_d11 <- if (comp_is_free(corr_model_path, "d2")) sqrt(extract_comp(corr_model_path, "d2")) else 0
    corr_path_e11 <- if (comp_is_free(corr_model_path, "e2")) sqrt(extract_comp(corr_model_path, "e2")) else 0

    # Extract proportions
    corr_path_a2 <- extract_comp(corr_model_path, "a2")
    corr_path_c2 <- extract_comp(corr_model_path, "c2")
    corr_path_d2 <- extract_comp(corr_model_path, "d2")
    corr_path_e2 <- extract_comp(corr_model_path, "e2")

    corr_path_total <- sum(c(corr_path_a2, corr_path_c2, corr_path_d2, corr_path_e2), na.rm = TRUE)

    corr_path_prop1 <- corr_path_a2 / corr_path_total
    corr_path_prop2 <- corr_path_c2 / corr_path_total
    corr_path_prop3 <- corr_path_d2 / corr_path_total
    corr_path_prop4 <- corr_path_e2 / corr_path_total

    # Add to combined results
    combined_results$Corr_Path_Est <- c(corr_path_a11, corr_path_c11, corr_path_d11, corr_path_e11,
                                       corr_path_prop1, corr_path_prop2, corr_path_prop3, corr_path_prop4)
  } else {
    combined_results$Corr_Path_Est <- NA
  }

  if ("ACE" %in% names(fits_var)) {
    corr_model_var <- fits_var[["ACE"]]

    # Extract variance components
    corr_var_A11 <- extract_comp(corr_model_var, "a2")
    corr_var_C11 <- extract_comp(corr_model_var, "c2")
    corr_var_D11 <- extract_comp(corr_model_var, "d2")
    corr_var_E11 <- extract_comp(corr_model_var, "e2")

    # Extract proportions
    corr_var_total <- sum(c(corr_var_A11, corr_var_C11, corr_var_D11, corr_var_E11), na.rm = TRUE)

    corr_var_prop1 <- corr_var_A11 / corr_var_total
    corr_var_prop2 <- corr_var_C11 / corr_var_total
    corr_var_prop3 <- corr_var_D11 / corr_var_total
    corr_var_prop4 <- corr_var_E11 / corr_var_total

    # Add to combined results
    combined_results$Corr_Var_Est <- c(corr_var_A11, corr_var_C11, corr_var_D11, corr_var_E11,
                                      corr_var_prop1, corr_var_prop2, corr_var_prop3, corr_var_prop4)
  } else {
    combined_results$Corr_Var_Est <- NA
  }
} else if (correlation_based_model == "ADE") {
  # Get ADE model results
  if ("ADE" %in% names(fits_path)) {
    corr_model_path <- fits_path[["ADE"]]

    # Extract path coefficients
    corr_path_a11 <- if (comp_is_free(corr_model_path, "a2")) sqrt(extract_comp(corr_model_path, "a2")) else 0
    corr_path_c11 <- if (comp_is_free(corr_model_path, "c2")) sqrt(extract_comp(corr_model_path, "c2")) else 0
    corr_path_d11 <- if (comp_is_free(corr_model_path, "d2")) sqrt(extract_comp(corr_model_path, "d2")) else 0
    corr_path_e11 <- if (comp_is_free(corr_model_path, "e2")) sqrt(extract_comp(corr_model_path, "e2")) else 0

    # Extract proportions
    corr_path_a2 <- extract_comp(corr_model_path, "a2")
    corr_path_c2 <- extract_comp(corr_model_path, "c2")
    corr_path_d2 <- extract_comp(corr_model_path, "d2")
    corr_path_e2 <- extract_comp(corr_model_path, "e2")

    corr_path_total <- sum(c(corr_path_a2, corr_path_c2, corr_path_d2, corr_path_e2), na.rm = TRUE)

    corr_path_prop1 <- corr_path_a2 / corr_path_total
    corr_path_prop2 <- corr_path_c2 / corr_path_total
    corr_path_prop3 <- corr_path_d2 / corr_path_total
    corr_path_prop4 <- corr_path_e2 / corr_path_total

    # Add to combined results
    combined_results$Corr_Path_Est <- c(corr_path_a11, corr_path_c11, corr_path_d11, corr_path_e11,
                                       corr_path_prop1, corr_path_prop2, corr_path_prop3, corr_path_prop4)
  } else {
    combined_results$Corr_Path_Est <- NA
  }

  if ("ADE" %in% names(fits_var)) {
    corr_model_var <- fits_var[["ADE"]]

    # Extract variance components
    corr_var_A11 <- extract_comp(corr_model_var, "a2")
    corr_var_C11 <- extract_comp(corr_model_var, "c2")
    corr_var_D11 <- extract_comp(corr_model_var, "d2")
    corr_var_E11 <- extract_comp(corr_model_var, "e2")

    # Extract proportions
    corr_var_total <- sum(c(corr_var_A11, corr_var_C11, corr_var_D11, corr_var_E11), na.rm = TRUE)

    corr_var_prop1 <- corr_var_A11 / corr_var_total
    corr_var_prop2 <- corr_var_C11 / corr_var_total
    corr_var_prop3 <- corr_var_D11 / corr_var_total
    corr_var_prop4 <- corr_var_E11 / corr_var_total

    # Add to combined results
    combined_results$Corr_Var_Est <- c(corr_var_A11, corr_var_C11, corr_var_D11, corr_var_E11,
                                      corr_var_prop1, corr_var_prop2, corr_var_prop3, corr_var_prop4)
  } else {
    combined_results$Corr_Var_Est <- NA
  }
} else {
  combined_results$Corr_Path_Est <- NA
  combined_results$Corr_Var_Est <- NA
}

# Add metadata
combined_results$MZ_cor <- mz_cor
combined_results$DZ_cor <- dz_cor
combined_results$Best_Path_Model <- if (!is.null(best_path_model)) best_path_model$name else NA
combined_results$Best_Var_Model <- if (!is.null(best_var_model)) best_var_model$name else NA
combined_results$Corr_Based_Model <- correlation_based_model

# Save combined results
combined_results_file <- paste0("estimates/combined_twin_results_", timestamp, ".csv")
write.csv(combined_results, combined_results_file, row.names = FALSE)
cat("\nCombined results saved to:", combined_results_file, "\n")

# Store results for potential use by wrapper scripts
mma <- list()
if (!is.null(mma_path)) mma$path <- mma_path
if (!is.null(mma_var)) mma$var <- mma_var
mma$combined_results <- combined_results
mma$best_path_model <- best_path_model
mma$best_var_model <- best_var_model
mma$correlation_based_model <- correlation_based_model
