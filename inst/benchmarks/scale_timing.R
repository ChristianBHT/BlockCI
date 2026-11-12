# Computation time of ONE test for pls_ci_test() and rosa_ci_test() as the data grow:
#   A) high N        : n = 1e3 ... 3e6 rows,       Z has 20 variables, X and Y 4 each
#   B) many Z vars   : p_z = 20 ... 5000 variables, n = 1000
#   C) many X/Y vars : q = 4 ... 1000 variables in X and in Y, n = 1000, p_z = 20
# Two settings per test (B = 199 permutations in both):
#   fast: ncomp = 12, in-sample residuals, simple permutation (no cross-validation)
#   cv  : ncomp chosen by CV (max 12), cross-validated residuals (5 segments), simple permutation
# Each size is run once, sequentially. A chain (study x setting x test) stops at the first
# run that takes longer than `cutoff` seconds, fails (e.g. out of memory), or would need more
# than `max_cells` numbers of data. Rows are appended to scale_timing_results.csv.
# Usage: Rscript scale_timing.R   (working directory = this folder; package CausalBlocks installed)

suppressMessages(library(CausalBlocks))
cutoff <- 120          # seconds
max_cells <- 6e7       # n * (p_z + q_x + q_y) numbers of data (6e7 doubles = 480 MB)
out <- "scale_timing_results.csv"
if (file.exists(out)) file.remove(out)

settings <- list(
  fast = list(ncomp = 12, resid = "fitted", perm = "simple", B = 199),
  cv   = list(ncomp = "auto", max_ncomp = 12, resid = "cv", perm = "simple", B = 199, segments = 5))
tests <- list(pls = pls_ci_test, rosa = rosa_ci_test)
studies <- list(
  A_high_N      = list(var = "n",   sizes = c(1e3, 1e4, 1e5, 3e5, 1e6, 3e6),
                       args = function(s) list(n = s, p_z = 20, q_x = 4, q_y = 4)),
  B_many_Z_vars = list(var = "p_z", sizes = c(20, 100, 500, 1000, 3000, 5000),
                       args = function(s) list(n = 1000, p_z = s, q_x = 4, q_y = 4)),
  C_many_XY_vars = list(var = "q",  sizes = c(4, 25, 100, 300, 1000),
                       args = function(s) list(n = 1000, p_z = 20, q_x = s, q_y = s)))

for (st in names(studies)) for (se in names(settings)) for (te in names(tests)) {
  S <- studies[[st]]
  for (s in S$sizes) {
    a <- S$args(s)
    cells <- a$n * (a$p_z + a$q_x + a$q_y)
    if (cells > max_cells) { cat(sprintf("%s %s %s %s=%g skipped (%.1e numbers of data)\n", st, se, te, S$var, s, cells)); break }
    set.seed(1)
    secs <- tryCatch({
      d <- do.call(sim_scale, a)
      gc(FALSE)
      p <- NA
      t <- system.time(p <- do.call(tests[[te]], c(list(Y ~ X | Z1 + Z2, d), settings[[se]]))$p.value)[["elapsed"]]
      rm(d); gc(FALSE)
      t
    }, error = function(e) { cat("  error:", conditionMessage(e), "\n"); NA_real_ })
    cat(sprintf("%s %-5s %-4s %s = %-8g : %s s\n", st, se, te, S$var, s, if (is.na(secs)) "NA" else format(round(secs, 2))))
    write.table(data.frame(study = st, setting = se, test = te, variable = S$var, size = s, n = a$n,
                           p_z = a$p_z, q = a$q_x, seconds = secs),
                out, sep = ",", row.names = FALSE, col.names = !file.exists(out), append = file.exists(out))
    if (is.na(secs) || secs > cutoff) break
  }
}
cat("done\n")
