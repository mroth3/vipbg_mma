library(MASS)  # for mvrnorm
library(parallel)

simulate_bivariate_twin_data <- function(n_mz_each = 250, n_dz_each = 250, n_dz_os = 250,
                                         a11, a22, a12,
                                         c11, c22, c12,
                                         e11, e22, e12) {
  zyg <- c(rep(1, n_mz_each), rep(2, n_dz_each), rep(3, n_mz_each), rep(4, n_dz_each), rep(5, n_dz_os))
  zyg2 <- ifelse(zyg %in% c(1, 3), 1, 2)
  sex1 <- c(rep(0, n_mz_each), rep(0, n_dz_each), rep(1, n_mz_each), rep(1, n_dz_each), rep(0, n_dz_os))
  sex2 <- c(rep(0, n_mz_each), rep(0, n_dz_each), rep(1, n_mz_each), rep(1, n_dz_each), rep(1, n_dz_os))

  A_cov <- matrix(c(a11, a12, a12, a22), 2, 2)
  C_cov <- matrix(c(c11, c12, c12, c22), 2, 2)
  E_cov <- matrix(c(e11, e12, e12, e22), 2, 2)

  simulate_pair <- function(zyg_code) {
    A1 <- mvrnorm(1, mu = c(0, 0), Sigma = A_cov)
    C1 <- mvrnorm(1, mu = c(0, 0), Sigma = C_cov)
    E1 <- mvrnorm(1, mu = c(0, 0), Sigma = E_cov)
    A2 <- if (zyg_code %in% c(1, 3)) A1 else mvrnorm(1, mu = c(0, 0), Sigma = A_cov * 0.5)
    C2 <- C1
    E2 <- mvrnorm(1, mu = c(0, 0), Sigma = E_cov)
    c(A1[1] + C1[1] + E1[1], A2[1] + C2[1] + E2[1], A1[2] + C1[2] + E1[2], A2[2] + C2[2] + E2[2])
  }

  trait_data <- t(sapply(zyg, simulate_pair))
  colnames(trait_data) <- c("p1_t1", "p2_t1", "p1_t2", "p2_t2")
  data.frame(zyg, zyg2, sex1, sex2, trait_data)
}

# Create folders
if (!dir.exists("dataset")) dir.create("dataset")
if (!dir.exists("truth")) dir.create("truth")

# Multiprocessing simulation
n_sims <- 100
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, varlist = ls())
parLapply(cl, 1:n_sims, function(i) {
  library(MASS)
  set.seed(as.integer(Sys.time()) + i)
  a11 <- runif(1, 0.2, 0.6)
  a22 <- runif(1, 0.2, 0.6)
  a12 <- runif(1, 0.0, min(a11, a22))
  c11 <- runif(1, 0.05, 0.3)
  c22 <- runif(1, 0.05, 0.3)
  c12 <- runif(1, 0.0, min(c11, c22))
  e11 <- runif(1, 0.2, 0.6)
  e22 <- runif(1, 0.2, 0.6)
  e12 <- runif(1, 0.0, min(e11, e22))

  dat <- simulate_bivariate_twin_data(
    a11 = a11, a22 = a22, a12 = a12,
    c11 = c11, c22 = c22, c12 = c12,
    e11 = e11, e22 = e22, e12 = e12
  )
  write.csv(dat, file = sprintf("dataset/data_%d.csv", i), row.names = FALSE)
  write.csv(data.frame(a11, a22, a12, c11, c22, c12, e11, e22, e12),
            file = sprintf("truth/truth_%d.csv", i), row.names = FALSE)
})
stopCluster(cl)
