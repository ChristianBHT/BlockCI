# Sanity checks: type I error under H0 and power under H1, auto ncomp, rotation vs permutation
suppressMessages(library(CausalBlocks))
sim <- function(n, effect) {
  V <- matrix(rnorm(n * 5), n)
  X <- V %*% matrix(rnorm(5 * 4), 5) + matrix(rnorm(n * 4), n)    # X depends on V
  Y <- cbind(V %*% rnorm(5) + rnorm(n), V %*% rnorm(5) + rnorm(n)) # Y depends on V
  Y[, 1] <- Y[, 1] + effect * (X[, 1] - V %*% qr.solve(V, X[, 1]))  # direct X -> Y (effect=0: H0)
  list(Y = Y, X = X, V = V)
}
set.seed(42); sims <- 30; R <- 39
t0 <- Sys.time()
res <- sapply(c(H0 = 0, H1_0.6 = 0.6, H1_1.5 = 1.5), function(e)
  sapply(c("rotation", "permutation"), function(m)
    mean(sapply(seq_len(sims), function(i) {
      d <- sim(200, e); SOPLS.CCI.test(Y ~ X | V, d, nperm = R, ncomp = "auto", method = m, progress = FALSE)$p.value
    }) < 0.05)))
cat("ncomp = auto, rejection rate at 0.05 (", sims, "sims, R =", R, ")\n"); print(round(res, 2))
cat("elapsed:", format(Sys.time() - t0), "\n")
d <- sim(200, 1.5); print(SOPLS.CCI.test(Y ~ X | V, d, nperm = 49, seed = 1, progress = FALSE))
cat("Unconditional (no V):\n"); print(SOPLS.CCI.test(Y ~ X, d, nperm = 49, seed = 1, progress = FALSE))
