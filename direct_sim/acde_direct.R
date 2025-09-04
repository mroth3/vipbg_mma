# ==============================================================================
# DIRECT ACDE TWIN DATA SIMULATION
# ==============================================================================
# Provides direct simulation of ACDE twin data using multivariate normal
# distributions, without rejection sampling. This is a more efficient approach
# for generating datasets with fixed ACDE proportions.

library(MASS)  # For mvrnorm

# ==============================================================================
# 1) Univariate ACDE data generator (direct, no accept/reject)
# ==============================================================================

acdeSim <- function(A, C, D, E, Nmz, Ndz, mean=0) {
  stopifnot(A >= 0, C >= 0, D >= 0, E >= 0)
  tot <- A + C + D + E
  if (tot <= 0) stop("All ACDE components are zero.")

  # Variance scale doesn't have to be 1, but this keeps things standard.
  # If you give standardized inputs, tot==1.
  mzMat <- matrix(c(tot, A + C + D, A + C + D, tot), 2, 2, byrow = TRUE)
  dzMat <- matrix(c(tot, 0.5*A + C + 0.25*D, 0.5*A + C + 0.25*D, tot), 2, 2, byrow = TRUE)

  mzData <- MASS::mvrnorm(n = Nmz, mu = c(mean, mean), Sigma = mzMat, empirical = FALSE)
  dzData <- MASS::mvrnorm(n = Ndz, mu = c(mean, mean), Sigma = dzMat, empirical = FALSE)

  colnames(mzData) <- colnames(dzData) <- c("t1","t2")

  list(mz = mzData, dz = dzData, mzCov = mzMat, dzCov = dzMat)
}


# ==============================================================================
# 2) Build a direct-variance ACDE model and run it (OpenMx)
# ==============================================================================

fitACDE_var <- function(dat, start = c(A=.3, C=.1, D=.05, E=.55), lb = c(A=0, C=0, D=0, E=1e-6)) {
  selVars <- c("t1","t2")
  MZdata <- mxData(dat$mz, type="raw")
  DZdata <- mxData(dat$dz, type="raw")

  A <- mxMatrix("Symm",1,1,free=TRUE, values=start["A"], lbound=lb["A"], labels="A11", name="A")
  C <- mxMatrix("Symm",1,1,free=TRUE, values=start["C"], lbound=lb["C"], labels="C11", name="C")
  D <- mxMatrix("Symm",1,1,free=TRUE, values=start["D"], lbound=lb["D"], labels="D11", name="D")
  E <- mxMatrix("Symm",1,1,free=TRUE, values=start["E"], lbound=lb["E"], labels="E11", name="E")

  Mean <- mxMatrix("Full",1,1,free=TRUE, values=0, labels="mu", name="Mean")
  expMean <- mxAlgebra(cbind(Mean,Mean), name="expMean")

  expMZ <- mxAlgebra(rbind(cbind(A+C+D+E, A+C+D), cbind(A+C+D, A+C+D+E)), name="expCovMZ")
  expDZ <- mxAlgebra(rbind(cbind(A+C+D+E, 0.5%x%A + C + 0.25%x%D), cbind(0.5%x%A + C + 0.25%x%D, A+C+D+E)), name="expCovDZ")

  mxMZ <- mxModel("MZ", A,C,D,E, Mean, expMean, expMZ, MZdata,
                  mxExpectationNormal(covariance="expCovMZ", means="expMean", dimnames=selVars),
                  mxFitFunctionML())

  mxDZ <- mxModel("DZ", A,C,D,E, Mean, expMean, expDZ, DZdata,
                  mxExpectationNormal(covariance="expCovDZ", means="expMean", dimnames=selVars),
                  mxFitFunctionML())

  top <- mxModel("ACDE_var", mxMZ, mxDZ, mxFitFunctionMultigroup(c("MZ","DZ")))

  suppressWarnings(mxRun(top, silent=TRUE))
}

# ==============================================================================
# 3) Convenience: fit ACE_var and ADE_var (also direct-variance)
# ==============================================================================

fitACE_var <- function(dat, start = c(A=.3, C=.2, E=.5)) {
  selVars <- c("t1","t2")
  MZdata <- mxData(dat$mz, type="raw")
  DZdata <- mxData(dat$dz, type="raw")

  A <- mxMatrix("Symm",1,1,free=TRUE, values=start["A"], lbound=0, labels="A11", name="A")
  C <- mxMatrix("Symm",1,1,free=TRUE, values=start["C"], lbound=0, labels="C11", name="C")
  E <- mxMatrix("Symm",1,1,free=TRUE, values=start["E"], lbound=1e-6, labels="E11", name="E")

  Mean <- mxMatrix("Full",1,1,free=TRUE, values=0, labels="mu", name="Mean")
  expMean <- mxAlgebra(cbind(Mean,Mean), name="expMean")

  expMZ <- mxAlgebra(rbind(cbind(A+C+E, A+C), cbind(A+C, A+C+E)), name="expCovMZ")
  expDZ <- mxAlgebra(rbind(cbind(A+C+E, 0.5%x%A + C), cbind(0.5%x%A + C, A+C+E)), name="expCovDZ")

  mxMZ <- mxModel("MZ", A,C,E, Mean, expMean, expMZ, MZdata,
                  mxExpectationNormal(covariance="expCovMZ", means="expMean", dimnames=selVars),
                  mxFitFunctionML())

  mxDZ <- mxModel("DZ", A,C,E, Mean, expMean, expDZ, DZdata,
                  mxExpectationNormal(covariance="expCovDZ", means="expMean", dimnames=selVars),
                  mxFitFunctionML())

  top <- mxModel("ACE_var", mxMZ, mxDZ, mxFitFunctionMultigroup(c("MZ","DZ")))

  suppressWarnings(mxRun(top, silent=TRUE))
}

fitADE_var <- function(dat, start = c(A=.3, D=.2, E=.5)) {
  selVars <- c("t1","t2")
  MZdata <- mxData(dat$mz, type="raw")
  DZdata <- mxData(dat$dz, type="raw")

  A <- mxMatrix("Symm",1,1,free=TRUE, values=start["A"], lbound=0, labels="A11", name="A")
  D <- mxMatrix("Symm",1,1,free=TRUE, values=start["D"], lbound=0, labels="D11", name="D")
  E <- mxMatrix("Symm",1,1,free=TRUE, values=start["E"], lbound=1e-6, labels="E11", name="E")

  Mean <- mxMatrix("Full",1,1,free=TRUE, values=0, labels="mu", name="Mean")
  expMean <- mxAlgebra(cbind(Mean,Mean), name="expMean")

  expMZ <- mxAlgebra(rbind(cbind(A+D+E, A+D), cbind(A+D, A+D+E)), name="expCovMZ")
  expDZ <- mxAlgebra(rbind(cbind(A+D+E, 0.5%x%A + 0.25%x%D), cbind(0.5%x%A + 0.25%x%D, A+D+E)), name="expCovDZ")

  mxMZ <- mxModel("MZ", A,D,E, Mean, expMean, expMZ, MZdata,
                  mxExpectationNormal(covariance="expCovMZ", means="expMean", dimnames=selVars),
                  mxFitFunctionML())

  mxDZ <- mxModel("DZ", A,D,E, Mean, expMean, expDZ, DZdata,
                  mxExpectationNormal(covariance="expCovDZ", means="expMean", dimnames=selVars),
                  mxFitFunctionML())

  top <- mxModel("ADE_var", mxMZ, mxDZ, mxFitFunctionMultigroup(c("MZ","DZ")))

  suppressWarnings(mxRun(top, silent=TRUE))
}

# ==============================================================================
# Example usage
# ==============================================================================
#
# # Simulate standardized ACDE = (.25,.25,.25,.25) with 500/750 pairs
# dat <- acdeSim(A=.25, C=.25, D=.25, E=.25, Nmz=500, Ndz=750)
#
# # Fit models
# fit_acde <- fitACDE_var(dat)
# fit_ace <- fitACE_var(dat)
# fit_ade <- fitADE_var(dat)
#
# # View results
# summary(fit_acde)
# summary(fit_ace)
# summary(fit_ade)
