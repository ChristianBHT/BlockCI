# Two conditioning sets on the same data: Y _||_ X | Z1,Z2   and   Y _||_ X | Z1
# Truth: first holds in all 4 scenarios. Second holds ONLY in scenario 2 (X <- Z1 -> Z2 -> Y);
# in scenarios 1, 3, 4 Z2 is a common cause, so it is false (power case).
# Usage (package root): Rscript inst/benchmarks/h0_conditioning_sets.R [reps] [n] [nperm]
args  <- commandArgs(TRUE)
reps  <- if (length(args) >= 1) as.integer(args[1]) else 100
n     <- if (length(args) >= 2) as.integer(args[2]) else 500
nperm <- if (length(args) >= 3) as.integer(args[3]) else 49

lib <- file.path(tempdir(), "lib"); dir.create(lib)
system2(file.path(R.home("bin"), "Rcmd.exe"),
        c("INSTALL", "--no-multiarch", paste0("--library=", shQuote(lib)), shQuote(getwd())),
        stdout = FALSE, stderr = FALSE)
.libPaths(c(lib, .libPaths())); Sys.setenv(R_LIBS = lib)
suppressMessages(library(CausalBlocks))

grid <- expand.grid(scenario = 1:4, cset = c("Z1,Z2", "Z1"), rep = seq_len(reps), stringsAsFactors = FALSE)
one <- function(i) {
  g <- grid[i, ]
  set.seed(100000 + 100 * g$scenario + g$rep)        # same data for both conditioning sets
  d <- CausalBlocks::sim_h0_scenario(g$scenario, n)
  f <- if (g$cset == "Z1,Z2") Y ~ X | Z1 + Z2 else Y ~ X | Z1
  r <- CausalBlocks::SOPLS.CCI.test(f, d, nperm = nperm, method = "rotation", progress = FALSE)
  cbind(g, n = n, H0_true = g$cset == "Z1,Z2" || g$scenario == 2, p.value = r$p.value)
}
cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1)); on.exit(parallel::stopCluster(cl))
parallel::clusterExport(cl, c("grid", "n", "nperm", "one"))
parallel::clusterEvalQ(cl, suppressMessages(library(CausalBlocks)))
res <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(nrow(grid)), one))
write.csv(res, "inst/benchmarks/h0_conditioning_sets_results.csv", row.names = FALSE)

summ <- do.call(rbind, lapply(split(res, list(res$scenario, res$cset)), function(d)
  data.frame(scenario = d$scenario[1], cset = d$cset[1], H0_true = d$H0_true[1], reps = nrow(d),
             reject05 = mean(d$p.value <= 0.05), mean_p = mean(d$p.value))))
rownames(summ) <- NULL; print(summ[order(summ$scenario, summ$cset), ], digits = 3)
