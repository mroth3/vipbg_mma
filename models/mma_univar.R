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

# Load Data for Relatives


dataMZm   <- subset(dataTwo,zyg==3,)
dataMZf   <- subset(dataTwo,zyg==1,)
dataDZm   <- subset(dataTwo,zyg==4,)
dataDZf   <- subset(dataTwo,zyg==2,)
dataDZo   <- subset(dataTwo,zyg==5,)

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
    cat(sprintf("Running lines up to %d/%d (%.1f%%)\n", current_line, total_lines, pct))

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

    # Break if we’ve reached the last line
    if (current_line >= total_lines) break
  }

  cat("✅ Finished running file:", file_path, "\n")
}

quiet_run('models/multimodel_fitUnivar.R')



