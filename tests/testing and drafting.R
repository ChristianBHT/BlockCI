# tests and drafts script

library(CausalBlocks)
sim_1 <- function(n = 500, p = 6) {
  Z1 <- matrix(stats::rnorm(n * p), n, p)
  Z2 <- matrix(stats::rnorm(n * p), n, p)
  Z  <- cbind(Z1, Z2)
  
  X  <- Z %*% matrix(stats::rnorm(2 * p * p), 2 * p, p) + matrix(stats::rnorm(n * p), n, p)
  Y  <- Z %*% matrix(stats::rnorm(2 * p * p), 2 * p, p) + matrix(stats::rnorm(n * p), n, p)
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
}
set.seed(1)
dat <- sim_1()
tic <- Sys.time()
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "pls")
toc <- Sys.time()
print(as.numeric(toc) - as.numeric(tic))

tic <- Sys.time()
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "rosa")
toc <- Sys.time()
print(as.numeric(toc) - as.numeric(tic))

Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "RF")
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "XGB")


library(grf)
library(grf)

# Cross-fitted prediction of M from Z: linear PLS first, then a multi-output
# random forest on the linear residuals. Returns the predictions for all rows,
# each made by models trained on the other folds.
pls_rf_crossfit <- function(M, Z, folds, num.trees = 2000, min.node.size = 3) {
  pred <- matrix(NA_real_, nrow(M), ncol(M))
  for (k in unique(folds)) {
    train <- which(folds != k)
    test  <- which(folds == k)
    
    # 1. linear part: PLS with all components (= OLS) on the training rows
    dtr <- data.frame(M = I(M[train, , drop = FALSE]), Z = I(Z[train, , drop = FALSE]))
    dte <- data.frame(Z = I(Z[test, , drop = FALSE]))
    lin <- pls::plsr(M ~ Z, ncomp = ncol(Z), data = dtr)
    lin_train <- predict(lin, ncomp = ncol(Z))[, , 1]                 # fitted, training rows
    lin_test  <- predict(lin, newdata = dte, ncomp = ncol(Z))[, , 1]  # held-out rows
    
    # 2. nonlinear remainder: forest on the training residuals of the linear fit
    rf <- multi_regression_forest(
      X = Z[train, , drop = FALSE],
      Y = M[train, , drop = FALSE] - lin_train,
      num.trees = num.trees,
      min.node.size = min.node.size
    )
    pred[test, ] <- lin_test + predict(rf, newdata = Z[test, , drop = FALSE])$predictions
  }
  pred
}

rv <- function(A, B) {
  Sab <- crossprod(A, B)
  sum(Sab^2) / sqrt(sum(crossprod(A)^2) * sum(crossprod(B)^2))
}

n_sim  <- 40
B      <- 1000
pvalues <- numeric(n_sim)

for (i in seq_len(n_sim)) {
  data <- sim_1(1000, p = 6)
  data$X <- scale(data$X, center = TRUE, scale = FALSE)
  data$Y <- scale(data$Y, center = TRUE, scale = FALSE)
  Z <- scale(cbind(data$Z1), center = TRUE, scale = FALSE)
  
  n <- nrow(data$X)
  folds <- sample(rep(1:5, length.out = n))      # same folds for X and Y
  
  resid_X <- data$X - pls_rf_crossfit(data$X, Z, folds)
  resid_Y <- data$Y - pls_rf_crossfit(data$Y, Z, folds)
  
  rv_obs <- rv(resid_X, resid_Y)
  null_dist <- replicate(B, rv(resid_X, resid_Y[sample(n), , drop = FALSE]))
  pvalues[i] <- (1 + sum(null_dist >= rv_obs)) / (B + 1)
  
  cat(sprintf("run %2d: RV = %.4f, null mean = %.4f, p = %.3f\n",
              i, rv_obs, mean(null_dist), pvalues[i]))
}

mean(pvalues < 0.05)   # rejection rate under H0: should be about 0.05



sim_non_linear <- function(n = 500, p = 6, linear = 1, nonlinear = 1,
                           noise_sd = 1, effect = 0) {
  Z1 <- matrix(stats::rnorm(n * p), n, p)
  Z2 <- matrix(stats::rnorm(n * p), n, p)
  Z  <- cbind(Z1, Z2)
  k  <- ncol(Z)                                   # 2p conditioning variables
  
  # nonlinear features of Z (shared by X and Y)
  a <- sample(k); b <- c(a[-1], a[1])             # k random pairs of different columns
  Fz <- cbind(Z^2,                                # quadratic
              sin(2 * Z),                         # oscillating
              (Z > 0) * 1,                        # threshold / step
              Z[, a] * Z[, b])                    # pairwise interactions
  Fz <- scale(Fz)                                 # unit variance per feature
  m  <- ncol(Fz)
  
  # random coefficients; both parts scaled to variance about k per column
  lin_part <- function() Z %*% matrix(stats::rnorm(k * p), k, p)
  nl_part  <- function() Fz %*% matrix(stats::rnorm(m * p), m, p) * sqrt(k / m)
  
  X <- linear * lin_part() + nonlinear * nl_part() +
    matrix(stats::rnorm(n * p, sd = noise_sd), n, p)
  Y <- linear * lin_part() + nonlinear * nl_part() +
    matrix(stats::rnorm(n * p, sd = noise_sd), n, p)
  
  if (effect != 0) Y <- Y + effect * X            # direct X -> Y effect: H0 false
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
}


dat <- sim_non_linear()
Block.CCI.test(Y ~ X | Z1 + Z2, data = dat, method = "RF",
               method_args = list(num.threads = NULL, num.trees = 500))
Block.CCI.test(Y ~ X | Z1, data = dat, method = "RF",
               method_args = list(num.threads = NULL, num.trees = 500))

dat <- sim_non_linear(effect = 1)
Block.CCI.test(Y ~ X | Z1 + Z2, data = dat, 
               method = "RF", 
               method_args = list(num.threads = NULL, num.trees = 500))

n_sim  <- 40
B      <- 1000
pvalues <- numeric(n_sim)
for (i in seq_len(n_sim)) {
  data <- sim_non_linear(1000, p = 6)
  data$X <- scale(data$X, center = TRUE, scale = FALSE)
  data$Y <- scale(data$Y, center = TRUE, scale = FALSE)
  Z <- scale(cbind(data$Z1, data$Z2), center = TRUE, scale = FALSE)
  
  n <- nrow(data$X)
  folds <- sample(rep(1:5, length.out = n))      # same folds for X and Y
  
  resid_X <- data$X - pls_rf_crossfit(data$X, Z, folds)
  resid_Y <- data$Y - pls_rf_crossfit(data$Y, Z, folds)
  
  rv_obs <- rv(resid_X, resid_Y)
  null_dist <- replicate(B, rv(resid_X, resid_Y[sample(n), , drop = FALSE]))
  pvalues[i] <- (1 + sum(null_dist >= rv_obs)) / (B + 1)
  
  cat(sprintf("run %2d: RV = %.4f, null mean = %.4f, p = %.3f\n",
              i, rv_obs, mean(null_dist), pvalues[i]))
}

mean(pvalues < 0.05)   # rejection rate under H0: should be about 0.05
