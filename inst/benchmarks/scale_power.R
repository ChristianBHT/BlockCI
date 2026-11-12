# Type I error and power of pls_ci_test() and rosa_ci_test() on large N and high-dimensional
# data with different noise (sim_scale(): gaussian, t3, laplace, skewed, hetero).
#   H0: effect = 0 (Y _||_ X | Z1, Z2 true)       H1: local alternative effect = 0.06 * sqrt(1000 / n)
# Cells:
#   "high N"   : n = 1e3, 1e4, 1e5, p_z = 20, setting "fast" (ncomp 12, in-sample residuals, simple perm)
#   "high dim" : n = 500, p_z = 200 and 1000, setting "cv" (CV-chosen ncomp, CV residuals, 5 segments)
#   "dense"    : n = 500, p_z = 200, signal in random directions of Z (hard case), setting "cv"
# B = 199 permutations. Results: scale_power_results.csv (one row per test) in the working directory.
# Usage: Rscript scale_power.R   (package CausalBlocks installed; uses `workers` cores)

suppressMessages(library(parallel))
library(CausalBlocks)
workers <- 6
alpha <- 0.05
effect_for <- function(n) 0.06 * sqrt(1000 / n)

fast <- list(ncomp = 12, resid = "fitted", perm = "simple", B = 199)
cv   <- list(ncomp = "auto", max_ncomp = 12, resid = "cv", perm = "simple", B = 199, segments = 5)
all_noise <- c("gaussian", "t3", "laplace", "skewed", "hetero")
sub_noise <- c("gaussian", "t3", "hetero")

grid <- rbind(
  expand.grid(cell = "high N",   n = 1e3, p_z = 20,   signal = "factor", noise = all_noise, reps = 80, setting = "fast", stringsAsFactors = FALSE),
  expand.grid(cell = "high N",   n = 1e4, p_z = 20,   signal = "factor", noise = all_noise, reps = 40, setting = "fast", stringsAsFactors = FALSE),
  expand.grid(cell = "high N",   n = 1e5, p_z = 20,   signal = "factor", noise = sub_noise, reps = 16, setting = "fast", stringsAsFactors = FALSE),
  expand.grid(cell = "high dim", n = 500, p_z = 200,  signal = "factor", noise = sub_noise, reps = 24, setting = "cv",   stringsAsFactors = FALSE),
  expand.grid(cell = "high dim", n = 500, p_z = 1000, signal = "factor", noise = sub_noise, reps = 24, setting = "cv",   stringsAsFactors = FALSE),
  expand.grid(cell = "dense",    n = 500, p_z = 200,  signal = "dense",  noise = "gaussian", reps = 24, setting = "cv",   stringsAsFactors = FALSE))

# one job = one dataset, both tests
jobs <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[rep(i, grid$reps[i] * 2), ]
  g$effect_on <- rep(c(FALSE, TRUE), each = grid$reps[i]); g$rep <- rep(seq_len(grid$reps[i]), 2); g
}))
jobs <- jobs[order(-jobs$n * jobs$p_z), ]                 # heavy jobs first (load balancing)
jobs$id <- seq_len(nrow(jobs))

run_job <- function(j, jobs, fast, cv, effect_for) {
  g <- jobs[j, ]
  set.seed(g$id * 7919L)
  d <- sim_scale(g$n, g$p_z, noise = g$noise, signal = g$signal,
                 effect = if (g$effect_on) effect_for(g$n) else 0)
  cfg <- if (g$setting == "fast") fast else cv
  f <- Y ~ X | Z1 + Z2
  do.call(rbind, lapply(c("pls", "rosa"), function(te) {
    fn <- if (te == "pls") pls_ci_test else rosa_ci_test
    t <- system.time(p <- do.call(fn, c(list(f, d), cfg))$p.value)[["elapsed"]]
    data.frame(g[, c("cell", "n", "p_z", "signal", "noise", "setting", "effect_on", "rep")], test = te, p.value = p, seconds = t)
  }))
}

cl <- makeCluster(workers)
clusterEvalQ(cl, suppressMessages(library(CausalBlocks)))
t0 <- Sys.time()
res <- do.call(rbind, parLapplyLB(cl, seq_len(nrow(jobs)), run_job, jobs = jobs, fast = fast, cv = cv, effect_for = effect_for))
stopCluster(cl)
write.csv(res, "scale_power_results.csv", row.names = FALSE)

tab <- aggregate(p.value ~ cell + n + p_z + signal + noise + test + effect_on, data = res,
                 FUN = function(p) c(rate = mean(p < alpha), reps = length(p)))
tab <- do.call(data.frame, tab)
names(tab)[names(tab) == "p.value.rate"] <- "reject"; names(tab)[names(tab) == "p.value.reps"] <- "reps"
tab$truth <- ifelse(tab$effect_on, "power (H1)", "type I (H0)")
print(tab[order(tab$cell, tab$n, tab$p_z, tab$noise, tab$test, tab$effect_on),
          c("cell", "n", "p_z", "signal", "noise", "test", "truth", "reps", "reject")], row.names = FALSE)
cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
