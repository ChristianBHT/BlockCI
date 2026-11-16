# Simulation study

devtools::install()
library(weightedGCM)
library(CausalBlocks)
simple_data <- function(n = 500, effect = 0){
  Z <- rnorm(n)
  X <- rnorm(n,Z)
  Y <- rnorm(n,Z + effect*X)
  data.frame(Y,X,Z)
}
dat <- simple_data()
Block.resid.test(Y ~ X | Z, data = dat, method = "PLS")
Block.resid.test(Y ~ X | Z, data = dat, method = "ROSA")
Block.resid.test(Y ~ X | Z, data = dat, method = "RF")
Block.resid.test(Y ~ X | Z, data = dat, method = "XGB")
Block.resid.test(Y ~ X | Z, 
                 data = dat, 
                 method = "RF", 
                 stat = "wgcm", 
                 k0 = 3, 
                 method_args = list(num.threads = NULL, num.trees = 500))

Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "PLS")
Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "SOPLS")
Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "ROSA")
Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "RF",
               method_args = list(num.threads = NULL, num.trees = 500))
Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "XGB",
               method_args = list(num.threads = NULL, num.trees = 500))

dat <- simple_data(effect = 1)
Block.resid.test(Y ~ X | Z, data = dat, method = "PLS")
Block.resid.test(Y ~ X | Z, data = dat, method = "ROSA")
Block.resid.test(Y ~ X | Z, data = dat, method = "RF")
Block.resid.test(Y ~ X | Z, data = dat, method = "XGB")
Block.resid.test(Y ~ X | Z, 
                 data = dat, 
                 method = "RF", 
                 stat = "wgcm", 
                 k0 = 3, 
                 method_args = list(num.threads = NULL, num.trees = 500))

Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "PLS")
Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "SOPLS")
Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "ROSA")
Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "RF",
               method_args = list(num.threads = NULL, num.trees = 500))
Block.CCI.test(Y ~ X | Z, 
               data = dat, 
               method = "XGB",
               method_args = list(num.threads = NULL, num.trees = 500))


linear_blocks <- function(n = 500, p = 6, effect = 0) {
  Z1 <- matrix(stats::rnorm(n * p), n, p)
  Z2 <- matrix(stats::rnorm(n * p), n, p)
  Z  <- cbind(Z1, Z2)
  
  X  <- Z %*% matrix(stats::rnorm(2 * p * p), 2 * p, p) + matrix(stats::rnorm(n * p), n, p)
  Y  <- Z %*% matrix(stats::rnorm(2 * p * p), 2 * p, p) + matrix(stats::rnorm(n * p), n, p)
  
  if (effect != 0) Y <- Y + effect * X  
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
}

dat <- linear_blocks()
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "PLS")
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "ROSA")
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "RF", 
                 method_args = list(num.threads = NULL, num.trees = 500))
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "XGB", 
                 method_args = list(num.threads = NULL, num.trees = 500))
Block.resid.test(Y ~ X | Z1 + Z2, 
                 data = dat, 
                 method = "RF", 
                 stat = "wgcm", 
                 k0 = 3, 
                 method_args = list(num.threads = NULL, num.trees = 500))

Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat, 
               method = "PLS")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat, 
               method = "SOPLS")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat, 
               method = "ROSA")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat, 
               method = "RF",
               method_args = list(num.threads = NULL, num.trees = 500))
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat, 
               method = "XGB",
               method_args = list(num.threads = NULL, num.trees = 500))

dat <- linear_blocks(effect = 1)
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "PLS")
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "ROSA")
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "RF", method_args = list(num.threads = NULL, num.trees = 500))
Block.resid.test(Y ~ X | Z1 + Z2, 
                 data = dat, 
                 method = "RF", 
                 stat = "wgcm", 
                 k0 = 3, 
                 method_args = list(num.threads = NULL, num.trees = 500))

Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat,
               nperm = 100,
               pvalue = "parametric",
               method = "PLS")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat,
               nperm = 100,
               pvalue = "parametric",
               method = "SOPLS")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat,
               nperm = 100,
               pvalue = "parametric",
               method = "ROSA")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat, 
               method = "RF",
               nperm = 100,
               pvalue = "parametric",
               method_args = list(num.threads = NULL, num.trees = 500))


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
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "PLS")
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "ROSA")
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "RF", method_args = list(num.threads = NULL, num.trees = 500))
Block.resid.test(Y ~ X | Z1 + Z2, 
                 data = dat, 
                 method = "RF", 
                 stat = "wgcm", 
                 k0 = 3, 
                 method_args = list(num.threads = NULL, num.trees = 500))


Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat,
               nperm = 100,
               pvalue = "parametric",
               method = "PLS")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat,
               nperm = 100,
               pvalue = "parametric",
               method = "SOPLS")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat,
               nperm = 100,
               pvalue = "parametric",
               method = "ROSA")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat, 
               method = "RF",
               nperm = 100,
               pvalue = "parametric",
               method_args = list(num.threads = NULL, num.trees = 500))

dat <- sim_non_linear(effect = 1)
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "PLS")
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "ROSA")
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, method = "RF", method_args = list(num.threads = NULL, num.trees = 500))
Block.resid.test(Y ~ X | Z1 + Z2, 
                 data = dat, 
                 method = "RF", 
                 stat = "wgcm", 
                 k0 = 3, 
                 method_args = list(num.threads = NULL, num.trees = 500))


Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat,
               nperm = 100,
               pvalue = "parametric",
               method = "PLS")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat,
               nperm = 100,
               pvalue = "parametric",
               method = "SOPLS")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat,
               nperm = 100,
               pvalue = "parametric",
               method = "ROSA")
Block.CCI.test(Y ~ X | Z1 + Z2, 
               data = dat, 
               method = "RF",
               nperm = 100,
               pvalue = "parametric",
               method_args = list(num.threads = NULL, num.trees = 500))


# weightedGCM is extremely slow!!!

X <- dat$X
Y <- dat$Y
Z_full <- cbind(dat$Z1, dat$Z2)

wgcm.fix(X, Y, Z_full, regr.meth = "xgboost") 
