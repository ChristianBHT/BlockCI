# Type I error of SOPLS.CCI.test in four H0 scenarios: Y _||_ X | Z1, Z2
# Usage (from the package root):  Rscript inst/benchmarks/h0_scenarios.R [reps] [nperm]
# Installs the package into a temporary library so the working copy is what is tested.
args  <- commandArgs(TRUE)
reps  <- if (length(args) >= 1) as.integer(args[1]) else 40
nperm <- if (length(args) >= 2) as.integer(args[2]) else 49
ns <- c(200, 500, 1000, 2000); scen <- 1:4; methods <- c("rotation", "permutation")

lib <- file.path(tempdir(), "lib"); dir.create(lib)
system2(file.path(R.home("bin"), "Rcmd.exe"),
        c("INSTALL", "--no-multiarch", paste0("--library=", shQuote(lib)), shQuote(getwd())),
        stdout = FALSE, stderr = FALSE)
.libPaths(c(lib, .libPaths())); Sys.setenv(R_LIBS = lib)
suppressMessages(library(CausalBlocks))

grid <- expand.grid(scenario = scen, n = ns, method = methods, rep = seq_len(reps), stringsAsFactors = FALSE)
grid <- grid[order(-grid$n), ]            # slowest jobs first for load balancing
one <- function(i) {
  g <- grid[i, ]; set.seed(1000 * g$scenario + 10 * g$rep + match(g$n, ns))
  d <- CausalBlocks::sim_h0_scenario(g$scenario, g$n)
  t0 <- Sys.time()
  r <- CausalBlocks::SOPLS.CCI.test(Y ~ X | Z1 + Z2, d, nperm = nperm, method = g$method, progress = FALSE)
  cbind(g, p.value = r$p.value, a = r$ncomp_selected[["V"]], b = r$ncomp_selected[["X"]],
        secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
nc <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(nc); on.exit(parallel::stopCluster(cl))
parallel::clusterExport(cl, c("grid", "ns", "nperm", "one"))
parallel::clusterEvalQ(cl, suppressMessages(library(CausalBlocks)))
res <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(grid)), one))
write.csv(res, "inst/benchmarks/h0_scenarios_results.csv", row.names = FALSE)

summ <- do.call(rbind, lapply(split(res, list(res$scenario, res$n, res$method), drop = TRUE), function(d)
  data.frame(scenario = d$scenario[1], n = d$n[1], method = d$method[1], reps = nrow(d),
             reject05 = mean(d$p.value <= 0.05), mean_p = mean(d$p.value),
             mean_b = mean(d$b), secs = mean(d$secs))))
rownames(summ) <- NULL; print(summ[order(summ$scenario, summ$method, summ$n), ], digits = 3)
