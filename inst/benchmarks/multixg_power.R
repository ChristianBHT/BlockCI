# Type I error and power of multixg_ci_test() on sim_1, sim_3 and sim_5.
# For every dataset the cross-fitted XGBoost residuals are computed once and the
# permutation test is run with both permutation types ("conditional", "simple").
#
# Truth (n rows, blocks Y, X, Z1, Z2):
#   sim_1: X, Y both generated from (Z1, Z2)        -> H0 true given Z1,Z2; false given Z1
#   sim_3: chain X <- Z1 -> Z2 -> Y                  -> H0 true given Z1,Z2 and given Z1; false given Z2
#   sim_5: nonlinear, heteroskedastic, from (Z1, Z2) -> H0 true given Z1,Z2; false given Z1
# Usage: Rscript multixg_power.R   (set reps, n, B below; uses all but one core)

suppressMessages(library(parallel))
library(CausalBlocks)

reps <- 50; n <- 500; B <- 199; alpha <- 0.05

cells <- list(
  list(sim = "sim_1", f = Y ~ X | Z1 + Z2, H0 = TRUE),
  list(sim = "sim_1", f = Y ~ X | Z1,      H0 = FALSE),
  list(sim = "sim_3", f = Y ~ X | Z1 + Z2, H0 = TRUE),
  list(sim = "sim_3", f = Y ~ X | Z1,      H0 = TRUE),
  list(sim = "sim_3", f = Y ~ X | Z2,      H0 = FALSE),
  list(sim = "sim_5", f = Y ~ X | Z1 + Z2, H0 = TRUE),
  list(sim = "sim_5", f = Y ~ X | Z1,      H0 = FALSE))

one_rep <- function(i, cells, n, B) {
  ns <- asNamespace("CausalBlocks")
  dat <- list()
  do.call(rbind, lapply(seq_along(cells), function(j) {
    ce <- cells[[j]]
    if (is.null(dat[[ce$sim]])) dat[[ce$sim]] <<- get(ce$sim, ns)(n)
    d <- dat[[ce$sim]]
    cc <- ns$.canonical(ce$f, d)
    Z <- do.call(cbind, cc$V)
    rs <- ns$.multixg_residuals(cc$X, cc$Y, Z, nfolds = 5, nrounds = 300,
                                early_stopping_rounds = 10, val_fraction = 0.2, params = list())
    data.frame(cell = j, rep = i,
               p_conditional = ns$.resid_perm_test(rs$res_x, rs$res_y, Z, "rv", "conditional", B, 10)$p,
               p_simple      = ns$.resid_perm_test(rs$res_x, rs$res_y, Z, "rv", "simple", B, 10)$p)
  }))
}

cl <- makeCluster(max(1, detectCores() - 1))
clusterEvalQ(cl, suppressMessages(library(CausalBlocks)))
clusterSetRNGStream(cl, 2026)
res <- do.call(rbind, parLapply(cl, seq_len(reps), one_rep, cells = cells, n = n, B = B))
stopCluster(cl)

tab <- do.call(rbind, lapply(seq_along(cells), function(j) {
  r <- res[res$cell == j, ]
  data.frame(scenario = cells[[j]]$sim, conditioning = deparse(cells[[j]]$f[[3]]),
             H0_true = cells[[j]]$H0, reps = nrow(r),
             reject_conditional = mean(r$p_conditional < alpha),
             reject_simple = mean(r$p_simple < alpha))
}))
print(tab, row.names = FALSE)
write.csv(res, file.path(tempdir(), "multixg_power_raw.csv"), row.names = FALSE)
