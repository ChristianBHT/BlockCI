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

dat <- sim_1()
multixg_ci_test(Y ~ X | Z1 + Z2, data = dat)



Block.CCI.test(Y ~ X | Z1 + Z2, data = dat, method = "XGB")
