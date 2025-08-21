# ==============================================================================
# Multimodel Averaging For VIPBG Simulated Twins Data Study
# ==============================================================================

# Load Libraries & Options
require(OpenMx)
require(psych)
source("models/miFunctions.R")
options(width=245)

# ---------------------------------------------------------------------------------------------------------------------|
# PREPARE DATA

# Load Data for Relatives - MODIFIED TO USE QUEUE FILE
allVars <- c('zyg','zyg2','sex1','sex2','p1_t1','p1_t2','p2_t1','p2_t2')

# Check if running in parallel mode (variables set by wrapper)
if (exists("current_file_path")) {
  # Running in parallel - use assigned file with retry logic
  max_retries <- 5
  retry_count <- 0
  dataTwo <- NULL

  while (retry_count < max_retries && is.null(dataTwo)) {
    tryCatch({
      dataTwo <- read.csv(current_file_path, header=T, na.strings=".", col.names=allVars)
    }, error = function(e) {
      if (grepl("cannot open the connection", e$message)) {
        retry_count <<- retry_count + 1
        if (retry_count < max_retries) {
          Sys.sleep(runif(1, 0.1, 1.0))  # Random delay before retry
        }
      } else {
        stop(e)  # Different error, don't retry
      }
    })
  }

  if (is.null(dataTwo)) {
    stop("Failed to read file after ", max_retries, " attempts: ", current_file_path)
  }

} else {
  # Running standalone - use default file (UPDATED PATH)
  dataTwo <- read.csv("datasets/data_a11_lg/simulated_data_1.csv", header=T, na.strings=".", col.names=allVars)
}

dataMZm <- subset(dataTwo,zyg==3,)
dataMZf <- subset(dataTwo,zyg==1,)
dataDZm <- subset(dataTwo,zyg==4,)
dataDZf <- subset(dataTwo,zyg==2,)
dataDZo <- subset(dataTwo,zyg==5,)

# CRITICAL FIX: Reorder data to match model expectations
dataMZm <- dataMZm[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]
dataMZf <- dataMZf[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]
dataDZm <- dataDZm[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]
dataDZf <- dataDZf[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]
dataDZo <- dataDZo[, c("p1_t1", "p2_t1", "p1_t2", "p2_t2")]

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
    pct <- round((current_line / total_lines) * 100, 1)
    # Suppress progress output when running in parallel
    if (!exists("current_file_path")) {
      cat(sprintf("Running lines up to %d/%d (%.1f%%)\n", current_line, total_lines, pct))
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

  if (!exists("current_file_path")) {
    cat("✅ Finished running file:", file_path, "\n")
  }
}

quiet_run('models/multimodel_fitFn.R')

# Test mxModelAverage(), with free parameters only (which is faster), for all permutations of arguments 'include', 'SE', 'refAsBlock', & 'type'
mma <- mxModelAverage(reference=c(colPra,"Prop"), models=fitACE5, include="onlyFree", SE=NULL, refAsBlock=FALSE, type="AIC")

# ---------------------------------------------------------------------------------------------------------------------|
# SAVE RESULTS (ADDED FOR PARALLEL PROCESSING)

# Only save results if running in parallel mode
if (exists("current_file_path") && exists("current_condition") && exists("current_dataset")) {

  # Create output directory structure
  output_base <- file.path("estimates", current_condition, paste0("dataset_", current_dataset))
  dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

  # Save complete MMA object
  saveRDS(mma, file.path(output_base, "complete_mma.rds"))

  # Save each component in separate folders
  for (comp_name in names(mma)) {

    # Clean component name for folder
    folder_name <- gsub("[^A-Za-z0-9_-]", "_", comp_name)
    comp_dir <- file.path(output_base, folder_name)
    dir.create(comp_dir, showWarnings = FALSE)

    # Save component data
    comp_data <- mma[[comp_name]]

    # Save as RDS
    saveRDS(comp_data, file.path(comp_dir, paste0(folder_name, ".rds")))

    # Save as CSV if possible
    if (is.data.frame(comp_data) || is.matrix(comp_data)) {
      write.csv(comp_data, file.path(comp_dir, paste0(folder_name, ".csv")), row.names = TRUE)
    }
  }

  # Create summary file
  summary_info <- data.frame(
    condition = current_condition,
    dataset = current_dataset,
    timestamp = Sys.time(),
    n_components = length(mma),
    component_names = paste(names(mma), collapse = "; "),
    stringsAsFactors = FALSE
  )

  write.csv(summary_info, file.path(output_base, "summary.csv"), row.names = FALSE)
}
