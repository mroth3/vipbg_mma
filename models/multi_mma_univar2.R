# ============================
# Proportion-based evaluation
# ============================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

# ---- folders ----
est_dir   <- "estimates"          # where combined_results CSVs live
sum_dir   <- file.path(est_dir, "summary")
plot_dir  <- file.path(est_dir, "plots")

dir.create(sum_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- truth lookup for variation_1..7 ----
# Map: prop1=A, prop2=C, prop3=D, prop4=E
truth_map <- list(
  variation_1 = c(prop1=.25, prop2=.25, prop3=.25, prop4=.25),
  variation_2 = c(prop1=.33, prop2=.33, prop3=0.00, prop4=.34),
  variation_3 = c(prop1=.33, prop2=0.00, prop3=.33, prop4=.34),
  variation_4 = c(prop1=.50, prop2=0.00, prop3=0.00, prop4=.50),
  variation_5 = c(prop1=.80, prop2=.10, prop3=0.00, prop4=.10),
  variation_6 = c(prop1=.80, prop2=0.00, prop3=.10, prop4=.10),
  variation_7 = c(prop1=.90, prop2=0.00, prop3=0.00, prop4=.10)
)

comp_labels <- c(prop1="A", prop2="C", prop3="D", prop4="E")

# Which estimator columns to evaluate (if present)
est_cols <- c("MM_Path_Est","MM_Var_Est","SM_Path_Est","SM_Var_Est",
              "Corr_Path_Est","Corr_Var_Est")

# -------------------------------
# Helper: read one result CSV and
# return long tibble of proportions
# -------------------------------
read_one_results <- function(fpath) {
  df <- suppressMessages(read_csv(fpath, show_col_types = FALSE))

  # Keep only the four proportion rows
  df_props <- df %>%
    filter(Parameter %in% c("prop1","prop2","prop3","prop4"))

  # Pivot to long for all estimator columns that exist
  avail <- intersect(est_cols, names(df_props))
  if (length(avail) == 0) {
    return(NULL)
  }

  df_long <- df_props %>%
    select(Parameter, all_of(avail)) %>%
    pivot_longer(cols = all_of(avail), names_to = "estimator", values_to = "estimate")

  df_long
}

# -------------------------------
# Helper: attach truth & metrics
# -------------------------------
attach_truth_and_metrics <- function(tbl, variation_name) {
  truth_vec <- truth_map[[variation_name]]
  if (is.null(truth_vec)) stop("Unknown variation: ", variation_name)

  tbl %>%
    mutate(
      comp   = recode(Parameter, !!!comp_labels),
      truth  = truth_vec[Parameter],
      error  = estimate - truth,
      abs_error = abs(error)
    )
}

# -------------------------------
# Aggregate RMSE etc. for one variation
# -------------------------------
evaluate_variation <- function(variation_dir) {
  variation_name <- basename(variation_dir)

  files <- list.files(variation_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(NULL)

  # Read all files -> bind
  dat <- files %>%
    set_names(basename) %>%
    map(read_one_results) %>%
    compact() %>%
    bind_rows(.id = "dataset")

  if (nrow(dat) == 0) return(NULL)

  dat2 <- dat %>%
    attach_truth_and_metrics(variation_name)

  # Per dataset,estimator,component
  # add RMSE later in a grouped summarize
  perf <- dat2 %>%
    group_by(estimator, comp) %>%
    summarize(
      n        = dplyr::n(),
      mean_est = mean(estimate, na.rm=TRUE),
      mean_truth = mean(truth, na.rm=TRUE),  # constant
      bias     = mean(estimate - truth, na.rm=TRUE),
      mae      = mean(abs(estimate - truth), na.rm=TRUE),
      rmse     = sqrt(mean((estimate - truth)^2, na.rm=TRUE)),
      sd_est   = sd(estimate, na.rm=TRUE),
      .groups = "drop"
    ) %>%
    mutate(variation = variation_name)

  list(raw = dat2, perf = perf)
}

# -------------------------------
# Run across estimates/<variation_*> folders
# (We assume your earlier pipeline wrote files to estimates/variation_X/)
# -------------------------------
variation_dirs <- list.dirs(est_dir, recursive = FALSE, full.names = TRUE)
variation_dirs <- variation_dirs[grepl("variation_", basename(variation_dirs))]

if (length(variation_dirs) == 0) {
  cat("No variation_* folders found under 'estimates/'. Nothing to do.\n")
  quit(save="no", status=0)
}

all_raw  <- list()
all_perf <- list()

for (vd in variation_dirs) {
  cat("Evaluating proportions in:", vd, "\n")
  res <- evaluate_variation(vd)
  if (is.null(res)) next
  all_raw[[basename(vd)]]  <- res$raw
  all_perf[[basename(vd)]] <- res$perf
}

raw_tbl  <- bind_rows(all_raw, .id = "variation_folder")
perf_tbl <- bind_rows(all_perf)

# Save tables
write_csv(perf_tbl, file.path(sum_dir, "proportion_performance_by_variation.csv"))
write_csv(raw_tbl,  file.path(sum_dir, "proportion_estimates_long.csv"))

# -------------------------------
# Plots (one per variation)
# -------------------------------
plot_violin_one <- function(raw_df, variation_name) {
  ggplot(raw_df, aes(x = comp, y = estimate, fill = estimator)) +
    geom_violin(scale = "width", alpha = 0.25, trim = FALSE) +
    geom_boxplot(width = 0.15, position = position_dodge(width = 0.9), outlier.shape = NA) +
    # truth lines by component
    geom_hline(data = distinct(raw_df, comp, truth),
               aes(yintercept = truth), linetype = "dashed") +
    labs(
      title = paste0("Proportion Estimates by Component — ", variation_name),
      x = "Component (A,C,D,E)",
      y = "Proportion",
      fill = "Estimator"
    ) +
    coord_cartesian(ylim = c(0,1)) +
    theme_minimal(base_size = 12)
}

plot_pointCI_one <- function(raw_df, variation_name) {
  summ <- raw_df %>%
    group_by(estimator, comp) %>%
    summarize(
      mean = mean(estimate, na.rm=TRUE),
      lo   = mean - 1.96 * sd(estimate, na.rm=TRUE) / sqrt(n()),
      hi   = mean + 1.96 * sd(estimate, na.rm=TRUE) / sqrt(n()),
      truth= first(truth),
      .groups = "drop"
    )

  ggplot(summ, aes(x = comp, y = mean, color = estimator, group = estimator)) +
    geom_point(position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15,
                  position = position_dodge(width = 0.4)) +
    geom_point(aes(y = truth), shape = 4, size = 3, color = "black",
               position = position_dodge(width = 0.4)) +
    labs(
      title = paste0("Mean ± 95% CI vs Truth — ", variation_name),
      x = "Component (A,C,D,E)",
      y = "Proportion",
      color = "Estimator"
    ) +
    coord_cartesian(ylim = c(0,1)) +
    theme_minimal(base_size = 12)
}

for (vn in unique(perf_tbl$variation)) {
  sub_raw <- raw_tbl %>% filter(variation == vn)

  p1 <- plot_violin_one(sub_raw, vn)
  p2 <- plot_pointCI_one(sub_raw, vn)

  ggsave(file.path(plot_dir, paste0(vn, "_violin_proportions.png")),
         p1, width = 9, height = 5, dpi = 300)
  ggsave(file.path(plot_dir, paste0(vn, "_meanCI_proportions.png")),
         p2, width = 9, height = 5, dpi = 300)
}

# -------------------------------
# A compact leaderboard (RMSE) across variations
# -------------------------------
leader <- perf_tbl %>%
  group_by(estimator, comp) %>%
  summarize(
    avg_rmse = mean(rmse, na.rm=TRUE),
    avg_bias = mean(bias, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  arrange(comp, avg_rmse)

write_csv(leader, file.path(sum_dir, "proportion_leaderboard_overall.csv"))

cat("\nSaved:\n",
    "- estimates/summary/proportion_performance_by_variation.csv\n",
    "- estimates/summary/proportion_estimates_long.csv\n",
    "- estimates/summary/proportion_leaderboard_overall.csv\n",
    "- estimates/plots/*_violin_proportions.png\n",
    "- estimates/plots/*_meanCI_proportions.png\n")
