# Type I error and power of the three Monte Carlo CCI tests that differ only in the
# regression: Block.CCI.test with method = "SOPLS", "PLS" (plain PLS) and
# "XGB" (multi-target XGBoost). Rotation of X, empirical p-value.
#
# Truth (n rows, blocks Y, X, Z1, Z2) -- same cells as multixg_power.R:
#   sim_1: X, Y both generated from (Z1, Z2)        -> H0 true given Z1,Z2; false given Z1
#   sim_3: chain X <- Z1 -> Z2 -> Y                  -> H0 true given Z1,Z2 and given Z1; false given Z2
#   sim_5: nonlinear, heteroskedastic, from (Z1, Z2) -> H0 true given Z1,Z2; false given Z1
# Usage: Rscript cci_learners_power.R   (uses up to 7 cores; writes cci_learners_power_results.csv)

suppressMessages(library(parallel))
suppressMessages(library(CausalBlocks))

n <- 500; alpha <- 0.05
settings <- list(sopls = list(reps = 40, B = 99),
                 pls   = list(reps = 40, B = 99),
                 xgb   = list(reps = 24, B = 49))

cells <- list(
  list(sim = "sim_1", f = Y ~ X | Z1 + Z2, H0 = TRUE),
  list(sim = "sim_1", f = Y ~ X | Z1,      H0 = FALSE),
  list(sim = "sim_3", f = Y ~ X | Z1 + Z2, H0 = TRUE),
  list(sim = "sim_3", f = Y ~ X | Z1,      H0 = TRUE),
  list(sim = "sim_3", f = Y ~ X | Z2,      H0 = FALSE),
  list(sim = "sim_5", f = Y ~ X | Z1 + Z2, H0 = TRUE),
  list(sim = "sim_5", f = Y ~ X | Z1,      H0 = FALSE))

one_task <- function(task, cells, n) {
  suppressMessages(library(CausalBlocks))
  ns <- asNamespace("CausalBlocks")
  B <- task$B
  dat <- list()
  do.call(rbind, lapply(seq_along(cells), function(j) {
    ce <- cells[[j]]
    if (is.null(dat[[ce$sim]])) dat[[ce$sim]] <<- get(ce$sim, ns)(n)
    d <- dat[[ce$sim]]
    t0 <- Sys.time()
    r <- suppressMessages(switch(task$method,
      sopls = Block.CCI.test(ce$f, d, method = "SOPLS", nperm = B, progress = FALSE),
      pls   = Block.CCI.test(ce$f, d, method = "PLS", nperm = B, progress = FALSE),
      xgb   = Block.CCI.test(ce$f, d, method = "XGB", nperm = B, progress = FALSE)))
    data.frame(method = task$method, cell = j, rep = task$rep, p = r$p.value, T = unname(r$statistic),
               secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }))
}

tasks <- do.call(c, lapply(names(settings), function(m)
  lapply(seq_len(settings[[m]]$reps), function(i) list(method = m, rep = i, B = settings[[m]]$B))))
tasks <- tasks[order(vapply(tasks, function(t) t$method != "xgb", TRUE))]   # slow tasks first

cl <- makeCluster(min(7, max(1, detectCores() - 1)))
clusterSetRNGStream(cl, 2027)
res <- do.call(rbind, parLapplyLB(cl, tasks, one_task, cells = cells, n = n))
stopCluster(cl)

tab <- do.call(rbind, lapply(names(settings), function(m) do.call(rbind, lapply(seq_along(cells), function(j) {
  r <- res[res$method == m & res$cell == j, ]
  data.frame(method = m, scenario = cells[[j]]$sim, conditioning = deparse(cells[[j]]$f[[3]]),
             H0_true = cells[[j]]$H0, reps = nrow(r), B = settings[[m]]$B,
             reject = mean(r$p < alpha), mean_secs = mean(r$secs))
}))))
print(tab, row.names = FALSE)
out <- if (dir.exists("inst/benchmarks")) "inst/benchmarks" else "."
write.csv(tab, file.path(out, "cci_learners_power_results.csv"), row.names = FALSE)
