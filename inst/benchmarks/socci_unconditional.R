# Unconditional case (no V): X and Y independent vs dependent
suppressMessages(library(CausalBlocks))
set.seed(3); n <- 200
X <- matrix(rnorm(n * 4), n); Y0 <- matrix(rnorm(n * 2), n); Y1 <- Y0; Y1[, 1] <- Y1[, 1] + X[, 1]
print(SOPLS.CCI.test(Y ~ X, list(Y = Y0, X = X), nperm = 49, seed = 1, progress = FALSE))
print(SOPLS.CCI.test(Y ~ X, list(Y = Y1, X = X), nperm = 49, seed = 1, progress = FALSE))
