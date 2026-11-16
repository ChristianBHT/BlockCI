# Simulation study: type I error and power of the block CI tests under different noise
# ----------------------------------------------------------------------------------------
# Data designs. X and Y depend only on the OBSERVED Z blocks plus their own noise, so with
# effect = 0 the null  Y _||_ X | Z1, Z2  holds exactly; effect > 0 adds Y <- Y + effect * X.
#   "linear"    : Z1, Z2 with 4 columns each; X and Y (4 columns) linear in Z.
#   "nonlinear" : as "linear" plus nonlinear functions of Z (squares, sin, step, pairwise
#                 interactions), shared by X and Y with their own random coefficients.
#   "highdim"   : Z1, Z2 with 100 columns each (200 Z variables, more than n at n = 200)
#                 from a 3-factor model; X and Y (10 columns) depend linearly and
#                 nonlinearly on the 3 factor directions of the observed Z.
# Noise on X and Y (all scaled to unit variance):
#   gaussian, t3 (heavy tails), laplace, skewed (centred exponential), outliers (5 % of
#   values from a normal with sd 5), hetero (sd depends on Z; the null still holds),
#   corr (the columns of the noise are AR(1)-correlated, rho = 0.6).
# Sample sizes 200, 500, 1000; effects 0 (type I error) and `effect_H1` (power).
#
# Methods (see `run_method()`):
#   CCI_RF, CCI_PLS, CCI_SOPLS : Block.CCI.test with method "RF", "PLS", "SOPLS"
#   RES_PLS, RES_RF            : Block.resid.test (RV statistic) with PLS, and PLS + RF
#   RES_RF_wgcm                : Block.resid.test, PLS + RF residuals, weighted GCM (k0 = 3)
#
# How to run: install the package (devtools::install()), set the working directory to
# where the results should go, then source this file. Section 1 is a quick check that
# everything runs (the simplest cases once); section 2 is the full study. Results are
# saved per dataset in `cfg$out_dir`, so an interrupted run continues where it stopped.
# Start with a pilot (cfg$reps = 5) to see the time per dataset before the full run.

library(CausalBlocks)
library(parallel)

cfg <- list(
  n_grid    = c(200, 500, 1000),
  designs   = c("linear", "nonlinear", "highdim"),
  noises    = c("gaussian", "t3", "laplace", "skewed", "outliers", "hetero", "corr"),
  effects   = c(0, 0.25),          # 0 = H0 (type I error), 0.25 = H1 (power)
  reps      = 100,                 # datasets per cell; MC s.e. of a 5 % rate: 0.022 (100), 0.015 (200)
  methods   = c("CCI_RF", "CCI_PLS", "CCI_SOPLS", "RES_PLS", "RES_RF", "RES_RF_wgcm"),
  nperm     = 99,                  # null draws in Block.CCI.test (p resolution 0.01)
  B         = 199,                 # permutations / Gaussian draws in Block.resid.test
  num_trees = 500,                 # trees per forest (RF methods)
  snr       = 2,                   # variance of the Z-part of each X / Y column (noise variance = 1)
  alpha     = 0.05,
  workers   = max(1, detectCores() - 1),
  out_dir   = "sim_noise_results",
  seed      = 20261116
)

# ---- data generators ----------------------------------------------------------------

# Unit-variance noise matrix (n x q). `s` is a standardised function of Z (mean of s^2
# about 1) that sets the noise level per row for type "hetero".
make_noise <- function(n, q, type, s) {
  m <- n * q
  e <- switch(type,
    gaussian = rnorm(m),
    t3       = rt(m, df = 3) / sqrt(3),
    laplace  = { u <- runif(m) - 0.5; -sign(u) * log(1 - 2 * abs(u)) / sqrt(2) },
    skewed   = rexp(m) - 1,
    outliers = ifelse(runif(m) < 0.05, rnorm(m, sd = 5), rnorm(m)) / sqrt(0.95 + 0.05 * 25),
    hetero   = rnorm(m) * sqrt(0.2 + 0.8 * s^2),   # s is recycled down every column: row i gets s_i
    corr     = rnorm(m),
    stop("unknown noise type ", type))
  E <- matrix(e, n, q)
  if (type == "corr") E <- E %*% chol(0.6^abs(outer(seq_len(q), seq_len(q), "-")))
  E
}

# Z-part of a block with q columns from features Fm (n x m, roughly unit variance):
# random coefficients scaled so that each column has variance about v.
mean_part <- function(Fm, q, v) {
  A <- matrix(rnorm(ncol(Fm) * q, sd = 1 / sqrt(ncol(Fm))), ncol(Fm), q)
  (Fm %*% A) * sqrt(v)
}

# Nonlinear features of Z: squares, sin(2z), steps at 0 and products of random pairs.
nl_features <- function(Z) {
  k <- ncol(Z)
  a <- sample(k); b <- c(a[-1], a[1])
  scale(cbind(Z^2, sin(2 * Z), (Z > 0) * 1, Z[, a, drop = FALSE] * Z[, b, drop = FALSE]))
}

gen_data <- function(design, n, noise, effect, snr = 2) {
  if (design == "highdim") {
    pz <- 100; q <- 10; r <- 3
    L  <- matrix(rnorm(r * 2 * pz), r, 2 * pz)                   # factor loadings
    Z  <- matrix(rnorm(n * r), n, r) %*% L + matrix(rnorm(n * 2 * pz), n, 2 * pz)
    U  <- t(L) / sqrt(rowSums(L^2))[col(t(L))]                   # unit factor directions
    S  <- scale(Z %*% U)                                         # 3 directions of the OBSERVED Z
    lin <- S; nl <- nl_features(S)
  } else {
    pz <- 4; q <- 4
    Z  <- matrix(rnorm(n * 2 * pz), n, 2 * pz)
    lin <- Z; nl <- if (design == "nonlinear") nl_features(Z) else NULL
  }
  # Z-part: the same features for X and Y, with their own random coefficients
  mu <- function() {
    if (is.null(nl)) mean_part(lin, q, snr) else mean_part(lin, q, snr / 2) + mean_part(nl, q, snr / 2)
  }
  s <- as.numeric(scale(Z %*% rnorm(ncol(Z))))                   # for "hetero" noise
  X <- mu() + make_noise(n, q, noise, s)
  Y <- mu() + make_noise(n, q, noise, s)
  if (effect != 0) Y <- Y + effect * X                           # direct X -> Y effect: H0 false
  list(Y = Y, X = X, Z1 = Z[, seq_len(pz), drop = FALSE], Z2 = Z[, pz + seq_len(pz), drop = FALSE])
}

# ---- one method on one dataset --------------------------------------------------------

run_method <- function(m, f, d, cfg) {
  rf_args <- list(num.trees = cfg$num_trees, num.threads = 1)
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(switch(m,
    CCI_RF      = Block.CCI.test(f, d, method = "RF", nperm = cfg$nperm,
                                 method_args = rf_args, progress = FALSE),
    CCI_PLS     = Block.CCI.test(f, d, method = "PLS", nperm = cfg$nperm, progress = FALSE),
    CCI_SOPLS   = Block.CCI.test(f, d, method = "SOPLS", nperm = cfg$nperm, progress = FALSE),
    RES_PLS     = Block.resid.test(f, d, method = "PLS", B = cfg$B),
    RES_RF      = Block.resid.test(f, d, method = "RF", B = cfg$B, nfolds = 5,
                                   method_args = rf_args),
    RES_RF_wgcm = Block.resid.test(f, d, method = "RF", stat = "wgcm", k0 = 3, B = cfg$B,
                                   nfolds = 5, method_args = rf_args),
    stop("unknown method ", m)),
    error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t0
  if (inherits(r, "error"))
    return(data.frame(method = m, p = NA_real_, p_par = NA_real_, seconds = secs,
                      error = conditionMessage(r)))
  data.frame(method = m, p = r$p.value,
             p_par = if (is.null(r$p.parametric)) NA_real_ else r$p.parametric,
             seconds = secs, error = NA_character_)
}

# ---- one dataset (job) with all methods; saved to its own file ---------------------------

run_job <- function(id, grid, cfg) {
  file <- file.path(cfg$out_dir, sprintf("job_%06d.rds", id))
  if (file.exists(file)) return(readRDS(file))
  g <- grid[grid$id == id, ]
  set.seed(cfg$seed + id)
  d <- gen_data(g$design, g$n, g$noise, g$effect, cfg$snr)
  res <- do.call(rbind, lapply(cfg$methods, run_method, f = Y ~ X | Z1 + Z2, d = d, cfg = cfg))
  res <- cbind(g[rep(1, nrow(res)), c("id", "design", "n", "noise", "effect", "rep")], res)
  rownames(res) <- NULL
  saveRDS(res, file)
  res
}

# =========================================================================================
# 1. Quick check: the simplest cases once, then one small dataset per design
# =========================================================================================
run_check <- TRUE
if (run_check) {
  cfg_small <- modifyList(cfg, list(nperm = 19, B = 49, num_trees = 100))

  # simplest case (univariate Y, X, Z), as in tests/Simulation study.R
  simple_data <- function(n = 500, effect = 0) {
    Z <- rnorm(n); X <- rnorm(n, Z); Y <- rnorm(n, Z + effect * X)
    data.frame(Y, X, Z)
  }
  set.seed(1)
  for (e in c(0, 1)) {
    cat("\nSimple univariate case, effect =", e, "\n")
    print(do.call(rbind, lapply(cfg$methods, run_method, f = Y ~ X | Z,
                                d = simple_data(300, e), cfg = cfg_small)))
  }

  # generators: dimensions and noise variance for every design and noise type
  cat("\nGenerators (noise variance of X given Z should be about 1):\n")
  for (des in cfg$designs) for (nz in cfg$noises) {
    set.seed(2)
    d0 <- gen_data(des, 1000, nz, effect = 0, snr = cfg$snr)
    Zall <- cbind(d0$Z1, d0$Z2)
    resid_var <- if (des == "linear") mean(apply(d0$X, 2, function(x) var(lm.fit(cbind(1, Zall), x)$residuals))) else NA
    cat(sprintf("  %-9s %-8s  Y %s  X %s  Z1 %s  Z2 %s  var(X) %.2f  noise var (linear only) %s\n",
                des, nz, paste(dim(d0$Y), collapse = "x"), paste(dim(d0$X), collapse = "x"),
                paste(dim(d0$Z1), collapse = "x"), paste(dim(d0$Z2), collapse = "x"),
                mean(apply(d0$X, 2, var)), if (is.na(resid_var)) "-" else sprintf("%.2f", resid_var)))
  }

  # all methods once per design (gaussian noise, n = 200, H0 true)
  for (des in cfg$designs) {
    set.seed(3)
    cat("\nDesign", des, "(n = 200, gaussian, H0 true):\n")
    print(do.call(rbind, lapply(cfg$methods, run_method, f = Y ~ X | Z1 + Z2,
                                d = gen_data(des, 200, "gaussian", 0, cfg$snr), cfg = cfg_small)))
  }
}

# =========================================================================================
# 2. Full study
# =========================================================================================
run_study <- FALSE      # set to TRUE when the check above looks right
if (run_study) {
  grid <- expand.grid(rep = seq_len(cfg$reps), effect = cfg$effects, noise = cfg$noises,
                      n = cfg$n_grid, design = cfg$designs, stringsAsFactors = FALSE)
  grid$id <- seq_len(nrow(grid))                               # fixed id: seed and file name
  order_ids <- grid$id[order(-(grid$n * ifelse(grid$design == "highdim", 10, 1)))]  # heavy jobs first
  dir.create(cfg$out_dir, showWarnings = FALSE)
  cat(nrow(grid), "datasets x", length(cfg$methods), "methods on", cfg$workers, "workers\n")

  cl <- makeCluster(cfg$workers)
  clusterEvalQ(cl, suppressMessages(library(CausalBlocks)))
  clusterExport(cl, c("gen_data", "make_noise", "mean_part", "nl_features", "run_method", "run_job"))
  t0 <- Sys.time()
  res <- do.call(rbind, parLapplyLB(cl, order_ids, run_job, grid = grid, cfg = cfg))
  stopCluster(cl)
  cat("Elapsed:", format(Sys.time() - t0), "\n")
  saveRDS(res, "sim_noise_results.rds")

  # ---- summary: rejection rates with Monte Carlo standard errors ----
  res$reject     <- res$p < cfg$alpha
  res$reject_par <- res$p_par < cfg$alpha
  key <- c("design", "noise", "n", "effect", "method")
  summ <- do.call(rbind, lapply(split(res, res[key], drop = TRUE), function(s) {
    k <- sum(!is.na(s$p)); r <- mean(s$reject, na.rm = TRUE)
    data.frame(s[1, key], reps = k, rate = r, mc_se = sqrt(r * (1 - r) / max(k, 1)),
               rate_par = mean(s$reject_par, na.rm = TRUE),
               mean_sec = mean(s$seconds), errors = sum(!is.na(s$error)))
  }))
  rownames(summ) <- NULL
  write.csv(summ, "sim_noise_summary.csv", row.names = FALSE)

  wide <- function(eff) {
    w <- reshape(summ[summ$effect == eff, c("design", "noise", "n", "method", "rate")],
                 idvar = c("design", "noise", "n"), timevar = "method", direction = "wide")
    names(w) <- sub("^rate\\.", "", names(w))
    w[order(w$design, w$noise, w$n), ]
  }
  cat("\nType I error (effect = 0), nominal", cfg$alpha, ":\n"); print(wide(0), digits = 2, row.names = FALSE)
  cat("\nPower (effect =", max(cfg$effects), "):\n"); print(wide(max(cfg$effects)), digits = 2, row.names = FALSE)

  # cells where the type I error is clearly above the nominal level
  bad <- summ[summ$effect == 0 & summ$rate > cfg$alpha + 2 * sqrt(cfg$alpha * (1 - cfg$alpha) / summ$reps), ]
  cat("\nCells with type I error above alpha + 2 s.e.:\n"); print(bad[, c(key, "rate", "reps")], row.names = FALSE)

  # ---- plot: rejection rate against n, one panel per design x noise ----
  pdf("sim_noise_plots.pdf", width = 11, height = 8.5)
  cols <- setNames(seq_along(cfg$methods), cfg$methods)
  for (eff in cfg$effects) for (des in cfg$designs) {
    op <- par(mfrow = c(2, 4), mar = c(4, 4, 2.5, 1))
    for (nz in cfg$noises) {
      s <- summ[summ$effect == eff & summ$design == des & summ$noise == nz, ]
      plot(NA, xlim = range(cfg$n_grid), ylim = if (eff == 0) c(0, 0.3) else c(0, 1),
           xlab = "n", ylab = "rejection rate", main = paste(des, nz, sep = ", "))
      abline(h = if (eff == 0) cfg$alpha else 1, lty = 2, col = "grey60")
      for (m in cfg$methods) {
        sm <- s[s$method == m, ]; sm <- sm[order(sm$n), ]
        lines(sm$n, sm$rate, type = "b", pch = 19, col = cols[m])
      }
    }
    plot.new(); legend("center", legend = cfg$methods, col = cols, lwd = 2, pch = 19, bty = "n",
                       title = if (eff == 0) "Type I error" else paste("Power, effect", eff))
    par(op)
  }
  dev.off()
}
