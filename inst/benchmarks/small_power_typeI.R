# Small simulation: type I error and power of the implemented CI tests.
# Usage (from package root): Rscript inst/benchmarks/small_power_typeI.R [reps] [nperm] [n] [sopls|resid]
# H0: Y _||_ X | Z1,Z2 (sim_2 linear, sim_5 nonlinear). H1: add effect * X-part to Y.
args  <- commandArgs(TRUE)
reps  <- if (length(args) >= 1) as.integer(args[1]) else 100
nperm <- if (length(args) >= 2) as.integer(args[2]) else 39
n     <- if (length(args) >= 3) as.integer(args[3]) else 300
set    <- if (length(args) >= 4) args[4] else "sopls"   # "sopls" or "resid" (PLS/ROSA residual tests)

lib <- file.path(tempdir(), "lib"); dir.create(lib)
system2(file.path(R.home("bin"), "Rcmd.exe"),
        c("INSTALL", "--no-multiarch", paste0("--library=", shQuote(lib)), shQuote(getwd())),
        stdout = FALSE, stderr = FALSE)
.libPaths(c(lib, .libPaths())); Sys.setenv(R_LIBS = lib)
suppressMessages(library(CausalBlocks))

scen    <- c(linear = "sim_2", nonlinear = "sim_5")
effects <- c(0, 0.3, 0.6)
tests   <- if (set == "resid") c("PLS-cond", "PLS-simple", "ROSA-cond", "ROSA-simple") else c("MC-rotation", "MC-permutation", "fast-NB", "fast-none")
out     <- if (set == "resid") "inst/benchmarks/small_power_typeI_resid_results.csv" else "inst/benchmarks/small_power_typeI_results.csv"
grid <- expand.grid(scenario = names(scen), effect = effects, rep = seq_len(reps),
                    test = tests, stringsAsFactors = FALSE)
grid$id <- as.integer(factor(paste(grid$scenario, grid$effect, grid$rep)))   # same data for all tests

one <- function(i) {
  g <- grid[i, ]
  set.seed(g$id)
  d <- get(scen[[g$scenario]], envir = asNamespace("CausalBlocks"))(n)
  if (g$effect > 0) {                      # direct X -> Y effect given Z
    W <- matrix(rnorm(ncol(d$X) * ncol(d$Y)), ncol(d$X), ncol(d$Y))
    S <- scale(d$X %*% W)
    d$Y <- d$Y + g$effect * sd(d$Y) * S * 1.5
  }
  set.seed(g$id + 1e6)
  f <- Y ~ X | Z1 + Z2
  t0 <- Sys.time()
  p <- switch(g$test,
    "MC-rotation"    = SOPLS.CCI.test(f, d, nperm = nperm, method = "rotation", progress = FALSE)$p.value,
    "MC-permutation" = SOPLS.CCI.test(f, d, nperm = nperm, method = "permutation", progress = FALSE)$p.value,
    "fast-NB"        = SOPLS.CCI.fast(f, d, nsplit = 10)$p.value,
    "fast-none"      = SOPLS.CCI.fast(f, d, nsplit = 10, correction = "none")$p.value,
    "PLS-cond"       = pls_ci_test(f, d, B = nperm, perm = "conditional")$p.value,
    "PLS-simple"     = pls_ci_test(f, d, B = nperm, perm = "simple")$p.value,
    "ROSA-cond"      = rosa_ci_test(f, d, B = nperm, perm = "conditional")$p.value,
    "ROSA-simple"    = rosa_ci_test(f, d, B = nperm, perm = "simple")$p.value)
  cbind(g, p.value = p, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
nc <- max(1, parallel::detectCores() - 2)
cl <- parallel::makeCluster(nc); on.exit(parallel::stopCluster(cl))
parallel::clusterExport(cl, c("grid", "scen", "n", "nperm", "one"))
parallel::clusterEvalQ(cl, suppressMessages(library(CausalBlocks)))
res <- do.call(rbind, parallel::parLapplyLB(cl, sample(nrow(grid)), one))
write.csv(res, out, row.names = FALSE)

summ <- do.call(rbind, lapply(split(res, list(res$scenario, res$effect, res$test), drop = TRUE), function(d) {
  r <- mean(d$p.value <= 0.05, na.rm = TRUE); m <- sum(!is.na(d$p.value))
  data.frame(scenario = d$scenario[1], effect = d$effect[1], test = d$test[1], reps = m,
             reject05 = r, mc_se = sqrt(r * (1 - r) / m), secs = mean(d$secs))
}))
rownames(summ) <- NULL
print(summ[order(summ$scenario, summ$effect, summ$test), ], digits = 3)
