# batch_acdeSim.R
# Generate 1,000 datasets for each of 7 ACDE targets and save to datasets/variation_k

# ---- user knobs with defaults ------------------------------------------------
# Default values
n_per <- 10   # datasets per variation
Nmz   <- 100    # MZ pairs per dataset
Ndz   <- 200    # DZ pairs per dataset

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
for (a in args) {
  if (grepl("--n_per=", a)) {
    n_per <- as.numeric(sub("--n_per=", "", a))
  }
  if (grepl("--Nmz=", a)) {
    Nmz <- as.numeric(sub("--Nmz=", "", a))
  }
  if (grepl("--Ndz=", a)) {
    Ndz <- as.numeric(sub("--Ndz=", "", a))
  }
}

cat("⚙️ Simulation settings:\n")
cat("   n_per =", n_per, "\n")
cat("   Nmz   =", Nmz, "\n")
cat("   Ndz   =", Ndz, "\n")
seed0    <- 20250902  # base seed; ensures reproducibility but uniqueness across runs
source('direct_acde/acde_direct.R')


# ---- targets (must sum to 1.00) ---------------------------------------------
targets <- list(
  # 1
  c(A=.25, C=.25, D=.25, E=.25),
  # 2
  c(A=.33, C=.33, D=0.00, E=.34),
  # 3
  c(A=.33, C=0.00, D=.33, E=.34),
  # 4
  c(A=.50, C=0.00, D=0.00, E=.50),
  # 5
  c(A=.80, C=.10, D=0.00, E=.10),
  # 6
  c(A=.80, C=0.00, D=.10, E=.10),
  # 7
  c(A=.90, C=0.00, D=0.00, E=.10)
)

# ---- helpers -----------------------------------------------------------------

ensure_dir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE) }

# Convert a univariate ACDE list (with mz, dz matrices having cols t1,t2) to dataTwo schema
to_dataTwo <- function(dat_list) {
  # dat_list$mz and dat_list$dz should be matrices or data.frames with columns named t1, t2
  mz <- as.data.frame(dat_list$mz)
  dz <- as.data.frame(dat_list$dz)

  # Defensive colnames
  colnames(mz) <- c("t1", "t2")
  colnames(dz) <- c("t1", "t2")

  # Build rows
  mz_df <- data.frame(
    zyg  = 1L, zyg2 = 1L,
    sex1 = 2L, sex2 = 2L,            # no sex differences programmed; keep fixed
    p1_t1 = mz$t1, p1_t2 = mz$t2,
    p2_t1 = NA_real_,                # univariate: second trait columns are NA
    p2_t2 = NA_real_
  )

  dz_df <- data.frame(
    zyg  = 2L, zyg2 = 2L,
    sex1 = 2L, sex2 = 2L,
    p1_t1 = dz$t1, p1_t2 = dz$t2,
    p2_t1 = NA_real_,
    p2_t2 = NA_real_
  )

  rbind(mz_df, dz_df)
}

# Single dataset generate-and-save
save_one_dataset <- function(k, i, A, C, D, E, Nmz, Ndz, out_dir, seed0) {
  # Unique seed per variation/dataset
  set.seed(seed0 + k*1e6 + i)

  # Call acdeSim (must return a list with mz, dz, mzCov, dzCov; mz/dz have cols t1, t2)
  dat <- acdeSim(A=A, C=C, D=D, E=E, Nmz=Nmz, Ndz=Ndz)

  # Convert to dataTwo format for CSV
  dataTwo <- to_dataTwo(dat)

  # Paths
  csv_path <- file.path(out_dir, sprintf("sim%d.csv", i))
  rds_path <- file.path(out_dir, sprintf("sim%d_dat.rds", i))

  # Write
  utils::write.csv(dataTwo, csv_path, row.names = FALSE)
  saveRDS(dat, rds_path)

  invisible(list(csv=csv_path, rds=rds_path))
}

# ---- main loop ----------------------------------------------------------------
# Base directory
out_dir  <- "datasets"
N <- Nmz+Ndz
# Construct subdirectory name: e.g. "datasets/10002500variation_1"
out_root <- file.path(out_dir, paste0("N","_",N))

# Create directory if it does not exist
if (!dir.exists(out_root)) {
  dir.create(out_root, recursive = TRUE)
  cat("Created output directory:", out_root, "\n")
} else {
  cat("Using existing output directory:", out_root, "\n")
}
ensure_dir(out_root)

manifest_all <- list()

for (k in seq_along(targets)) {
  tgt <- targets[[k]]
  stopifnot(abs(sum(tgt) - 1) < 1e-8)

  variation_dir <- file.path(out_root, sprintf("variation_%d", k))
  ensure_dir(variation_dir)

  message(sprintf("[variation_%d] Target: A=%.2f C=%.2f D=%.2f E=%.2f  -> %s",
                  k, tgt["A"], tgt["C"], tgt["D"], tgt["E"], variation_dir))

  # Manifest entries
  items <- vector("list", n_per)

  for (i in seq_len(n_per)) {
    out <- save_one_dataset(
      k=k, i=i,
      A=tgt["A"], C=tgt["C"], D=tgt["D"], E=tgt["E"],
      Nmz=Nmz, Ndz=Ndz,
      out_dir=variation_dir, seed0=seed0
    )
    items[[i]] <- c(sim=i, csv=out$csv, rds=out$rds)
    if (i %% 50 == 0) message(sprintf("  ... %d / %d saved", i, n_per))
  }

  manifest <- list(
    target = unclass(tgt),
    Nmz = Nmz, Ndz = Ndz,
    n_datasets = n_per,
    files = items
  )

  # Save a JSON and CSV manifest for convenience
  jsonlite::write_json(manifest, file.path(variation_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)

  # Also a compact CSV index
  idx <- data.frame(
    sim = seq_len(n_per),
    csv = vapply(items, function(x) x["csv"], character(1)),
    rds = vapply(items, function(x) x["rds"], character(1)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(idx, file.path(variation_dir, "manifest.csv"), row.names = FALSE)

  manifest_all[[paste0("variation_", k)]] <- manifest
}

# Top-level manifest
jsonlite::write_json(manifest_all, file.path(out_root, "manifest_all.json"), auto_unbox = TRUE, pretty = TRUE)

message("\n✅ Done. Datasets are in 'datasets/variation_1' ... 'variation_7'.")
message("   Each folder contains sim1.csv ... sim1000.csv and sim*_dat.rds plus a manifest.")
