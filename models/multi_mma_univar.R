# evaluate_all_versions.R
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(ggplot2)
})

# ----------------------------
# 1) Paths & IO setup
# ----------------------------
base_dir <- "datasets/N_300"
out_dir  <- "estimates"
plot_dir <- file.path(out_dir, "plots")

dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Find variation folders
variation_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
variation_dirs <- variation_dirs[grepl("variation_[1-7]$", variation_dirs)]

if (length(variation_dirs) == 0) {
  stop("No variation_* folders found under 'datasets/'.")
}

# ----------------------------
# 2) Truth map for the 7 targets
# (ensure they sum to 1; your note uses 0.34 where needed)
# ----------------------------
truth_map <- list(
  variation_1 = c(A=.25, C=.25, D=.25, E=.25),
  variation_2 = c(A=.33, C=.33, D=0.00, E=.34),
  variation_3 = c(A=.33, C=0.00, D=.33, E=.34),
  variation_4 = c(A=.50, C=0.00, D=0.00, E=.50),
  variation_5 = c(A=.80, C=.10, D=0.00, E=.10),
  variation_6 = c(A=.80, C=0.00, D=.10, E=.10),
  variation_7 = c(A=.90, C=0.00, D=0.00, E=.10)
)

# helper: map prop param -> component
prop_to_comp <- function(p) {
  switch(p,
         prop1 = "A",
         prop2 = "C",
         prop3 = "D",
         prop4 = "E",
         NA_character_)
}

# ----------------------------
# 3) Per-file runner
# ----------------------------
run_one_file <- function(csv_path, truth_vec, out_dir) {
  # Load CSV to "dataTwo" expected format
  # NOTE: your models/mma_univar2.R reads a global dataTwo
  # and creates "combined_results"
  allVars <- c('zyg','zyg2','sex1','sex2','p1_t1','p1_t2','p2_t1','p2_t2')
  dataTwo <<- read.csv(csv_path, header = TRUE, col.names = allVars)

  # Extract phenotype data for later plotting
  phenotype_data <- dataTwo %>%
    select(p1_t1, p1_t2) %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = recode(variable, p1_t1 = "t1", p1_t2 = "t2"),
           dataset = basename(csv_path))

  # run your existing analysis (must assign combined_results)
  source('models/mma_univar2.R', local = TRUE)

  if (!exists("combined_results")) {
    warning("combined_results not found after running mma_univar.R for: ", csv_path)
    return(list(long_data = NULL, phenotype_data = phenotype_data))
  }

  # Save per-dataset results with the same filename in estimates/
  dataset_name <- basename(csv_path)
  out_file     <- file.path(out_dir, dataset_name)
  write.csv(combined_results, out_file, row.names = FALSE)

  # Tidy and attach truth for bias calculations
  # We calculate metrics for proportion scale only:
  #   - PROPORTION scale: prop1,prop2,prop3,prop4 vs truth component
  cr <- combined_results

  # MZ/DZ cor (grab once)
  get_cor <- function(col) if (col %in% names(cr)) cr[[col]][1] else NA_real_
  mz_cor <- get_cor("MZ_cor")
  dz_cor <- get_cor("DZ_cor")

  # Build long frame for proportion scale (prop1..prop4)
  prop_long <- cr %>%
    filter(Parameter %in% c("prop1","prop2","prop3","prop4")) %>%
    transmute(
      Parameter,
      component = vapply(Parameter, prop_to_comp, character(1)),
      # Use MM_Path_Est or MM_Var_Est based on which contains the proportion estimates
      # They should be equal in the example, but we'll check both
      MM_est = if ("MM_Path_Est" %in% names(.)) MM_Path_Est else MM_Var_Est,
      MM_se  = if ("MM_Path_SE" %in% names(.)) MM_Path_SE else MM_Var_SE,
      SM_est = if ("SM_Var_Est" %in% names(.)) SM_Var_Est else NA_real_,
      SM_se  = if ("SM_Var_SE" %in% names(.)) SM_Var_SE else NA_real_,
      CB_est = if ("Corr_Var_Est" %in% names(.)) Corr_Var_Est else NA_real_,
      scale  = "proportion"
    ) %>%
    mutate(
      truth      = truth_vec[component],
      file       = dataset_name,
      MZ_cor     = mz_cor,
      DZ_cor     = dz_cor
    )

  return(list(long_data = prop_long, phenotype_data = phenotype_data))
}

# ----------------------------
# 4) Loop all variations/files, compute metrics, plot
# ----------------------------
all_results <- list()
all_phenotype_data <- list()

for (vdir in sort(variation_dirs)) {
  vname <- basename(vdir)
  if (!vname %in% names(truth_map)) {
    message("Skipping unknown variation folder: ", vname)
    next
  }
  cat("\n== Processing ", vname, " ==\n", sep="")

  truth_vec <- truth_map[[vname]]
  files <- list.files(vdir, pattern = "^sim\\d+\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    message("No sim*.csv files in ", vdir)
    next
  }

  v_out_dir  <- file.path(out_dir, vname)
  v_plot_dir <- file.path(plot_dir, vname)
  dir.create(v_out_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(v_plot_dir, recursive = TRUE, showWarnings = FALSE)

  # Run all files for this variation and collect results
  v_results <- map(files, ~run_one_file(.x, truth_vec, v_out_dir))

  # Extract long data and phenotype data
  v_res <- bind_rows(lapply(v_results, function(x) x$long_data))
  v_phenotype <- bind_rows(lapply(v_results, function(x) x$phenotype_data))

  if (nrow(v_res) == 0) next

  # Long-to-long: melt methods for one tidy table
  v_long <- v_res %>%
    pivot_longer(cols = c(MM_est, SM_est, CB_est),
                 names_to = "method", values_to = "estimate") %>%
    mutate(method = recode(method,
                           MM_est = "Multimodel Average",   # Multi-model average
                           SM_est = "Single Model",   # Single best model
                           CB_est = "rMZ/rDZ"  # correlation-based choice
    ))

  # Add the standard errors directly in the pivot operation to avoid join issues
  v_long <- v_long %>%
    mutate(se = case_when(
      method == "MM" & !is.na(MM_se) ~ MM_se,
      method == "SM" & !is.na(SM_se) ~ SM_se,
      TRUE ~ NA_real_
    )) %>%
    select(file, Parameter, component, scale, method, estimate, se, truth, MZ_cor, DZ_cor)

  # Compute metrics
  v_metrics <- v_long %>%
    group_by(component, method) %>%
    summarize(
      n         = sum(!is.na(estimate)),
      bias_mean = mean(estimate - truth, na.rm = TRUE),
      mae       = mean(abs(estimate - truth), na.rm = TRUE),
      rmse      = sqrt(mean((estimate - truth)^2, na.rm = TRUE)),
      cor_est   = suppressWarnings(cor(estimate, truth, use = "complete.obs")),
      .groups = "drop"
    ) %>%
    mutate(variation = vname)

  # Save per-variation summaries
  write_csv(v_long,    file.path(v_out_dir,  paste0(vname, "_long_estimates.csv")))
  write_csv(v_metrics, file.path(v_out_dir,  paste0(vname, "_metrics.csv")))

  # Store for global aggregation
  all_results[[vname]] <- list(long = v_long, metrics = v_metrics)
  all_phenotype_data[[vname]] <- v_phenotype

  # ----------------------------
  # Plots (per variation)
  # ----------------------------

  # 1) Violin plots of estimates vs truth (faceted by component)
  p_violin <- v_long %>%
    ggplot(aes(x = method, y = estimate, fill = method)) +
    geom_violin(trim = FALSE, alpha = .25) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = .5) +
    geom_hline(aes(yintercept = truth), color = "black", linetype = "dashed") +
    facet_wrap(~ component, ncol = 4) +
    labs(title = paste0("Proportion estimates by method \u2014 ", vname),
         y = "Proportion", x = "Method") +
    theme_bw()

  ggsave(file.path(v_plot_dir, paste0(vname, "_violin_estimates.png")),
         p_violin, width = 12, height = 6, dpi = 200)

  # 2) Bias violin (estimate - truth)
  p_bias <- v_long %>%
    mutate(bias = estimate - truth) %>%
    ggplot(aes(x = method, y = bias, fill = method)) +
    geom_violin(trim = FALSE, alpha = .25) +
    geom_boxplot(width = .15, outlier.shape = NA, alpha = .5) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    facet_wrap(~ component, ncol = 4) +
    labs(title = paste0("Proportion bias (estimate - truth) \u2014 ", vname),
         y = "Bias", x = "Method") +
    theme_bw()

  ggsave(file.path(v_plot_dir, paste0(vname, "_violin_bias.png")),
         p_bias, width = 12, height = 6, dpi = 200)

  # 3) RMSE bars
  p_rmse <- v_metrics %>%
    ggplot(aes(x = method, y = rmse, fill = method)) +
    geom_col(alpha = .8) +
    facet_wrap(~ component, ncol = 4) +
    labs(title = paste0("Proportion RMSE by method \u2014 ", vname),
         y = "RMSE", x = "Method") +
    theme_bw()

  ggsave(file.path(v_plot_dir, paste0(vname, "_rmse.png")),
         p_rmse, width = 12, height = 6, dpi = 200)

  # 4) Phenotype histograms + density
  if (nrow(v_phenotype) > 0) {
    p_phenotype_hist <- v_phenotype %>%
      ggplot(aes(x = value)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.7) +
      geom_density(color = "darkred", linewidth = 1) +
      facet_wrap(~ variable, ncol = 1) +
      labs(title = paste0("Phenotype distributions \u2014 ", vname),
           x = "Value", y = "Density") +
      theme_bw()

    ggsave(file.path(v_plot_dir, paste0(vname, "_phenotype_hist_density.png")),
           p_phenotype_hist, width = 10, height = 8, dpi = 200)

    # 5) Phenotype boxplots
    p_phenotype_box <- v_phenotype %>%
      ggplot(aes(x = variable, y = value, fill = variable)) +
      geom_boxplot() +
      labs(title = paste0("Phenotype distributions \u2014 ", vname),
           x = "Variable", y = "Value") +
      theme_bw()

    ggsave(file.path(v_plot_dir, paste0(vname, "_phenotype_box.png")),
           p_phenotype_box, width = 8, height = 6, dpi = 200)
  }

  cat("  Saved: ",
      file.path(v_out_dir,  paste0(vname, "_metrics.csv")), " and plots in ",
      v_plot_dir, "\n", sep="")
}

# ----------------------------
# 5) Global summary across variations
# ----------------------------
if (length(all_results) > 0) {
  all_long    <- bind_rows(lapply(all_results, `[[`, "long"), .id = "variation")
  all_metrics <- bind_rows(lapply(all_results, `[[`, "metrics"))

  write_csv(all_long,    file.path(out_dir, "ALL_long_estimates.csv"))
  write_csv(all_metrics, file.path(out_dir, "ALL_metrics.csv"))

  # Global RMSE plot across all variations
  p_global_rmse <- all_metrics %>%
    ggplot(aes(x = method, y = rmse, fill = method)) +
    geom_col(alpha = .85, position = position_dodge(width = .7)) +
    facet_wrap(~ component, ncol = 4) +
    labs(title = "Global Proportion RMSE across all variations",
         y = "RMSE", x = "Method") +
    theme_bw()

  ggsave(file.path(plot_dir, "ALL_rmse.png"),
         p_global_rmse, width = 12, height = 6, dpi = 200)

  cat("\n== DONE ==\n",
      "Saved overall tables in 'estimates/' and plots in 'estimates/plots/'.\n")
} else {
  cat("\nNo results were produced. Check dataset folders and model script paths.\n")
}
