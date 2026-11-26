# Simulation studies: type I error and power of the block CI tests under different noise
# ----------------------------------------------------------------------------------------
# One function per data-generating process (a "study"):
#
#   study_<name>(testing_method, effect = 0, n = 500, n_tests = 100, ...)
#
# generates `n_tests` datasets with `n` rows, runs the test(s) on each, shows a progress
# bar with the speed and the estimated time left, and reports the rejection rate (with
# its Monte Carlo standard error), the median time per test and any errors.
#
#   testing_method  one or more of the names in `test_methods` below (e.g. "CCI_RF"), or
#                   a function(formula, data, settings) returning an object with $p.value.
#                   Several methods are run on the SAME datasets (paired comparison).
#   effect          0: the null Y _||_ X | Z1, Z2 holds (type I error);
#                   > 0: Y <- Y + effect * X, the null is false (power).
#   n, n_tests      rows per dataset, number of datasets.
#   ...             passed to run_study(): design ("linear" or "nonlinear"), alpha,
#                   settings (list overriding `default_settings`), seed, save_dir, verbose.
#
# Studies: study_gaussian_noise, study_t3_noise, study_laplace_noise, study_skewed_noise,
# study_poisson_noise, study_outlier_noise, study_hetero_noise, study_correlated_noise,
# study_highdim (nonlinear = TRUE: the limit case; FALSE: the linear high-dimensional
# check, also as study_highdim_linear). In every study X and Y depend only on the OBSERVED Z blocks plus their
# own noise, so with effect = 0 the null holds exactly. All noise has variance 1.
#
# Each call returns (invisibly) list(results = one row per test, summary = one row per
# method, settings). With save_dir it is also saved as an .rds file; collect_studies()
# reads all saved files into one table. A relative save_dir is a folder inside
# results_root ("simulation_results/", ignored by Git and by the package build), so run
# from the project folder. To run several studies at the same time, start them as
# separate RStudio Background Jobs (one R process each).
#
# Usage: install the package (devtools::install()), source this file, then call the study
# functions (examples at the bottom).

library(CausalBlocks)

# ---- where results are saved -----------------------------------------------------------

results_root <- "simulation_results"

# a relative folder name lies inside results_root; a full path is used as it is
results_path <- function(dir) {
  if (grepl("^([A-Za-z]:|/|\\\\|~)", dir)) dir else file.path(results_root, dir)
}

# ---- settings and test methods ---------------------------------------------------------

default_settings <- list(
  nperm     = 99,      # null draws in Block.CCI.test (p resolution 0.01)
  p         = 0.5,     # training fraction in Block.CCI.test (older results: 0.75)
  npval     = 1,       # observed splits whose p-values are averaged (> 1: experimental)
  B         = 199,     # permutations / Gaussian draws in Block.resid.test
  num_trees = 500,     # trees per forest (RF methods)
  threads   = NULL,    # forest threads: NULL = all cores (one test at a time); 1 in parallel jobs
  nfolds    = 5,       # cross-fitting folds in Block.resid.test
  k0        = 3,       # weight functions per Z column for the weighted GCM
  k_perm    = 5        # neighbours in the local permutation (resampling = "conditional")
)

test_methods <- list(
  CCI_RF      = function(f, d, s) Block.CCI.test(f, d, method = "RF", nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE,
                                                 method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  CCI_PLS     = function(f, d, s) Block.CCI.test(f, d, method = "PLS", nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE),
  # two-fold cross-fitted statistic and null draws (experimental; always p = 0.5, npval = 1)
  CCI_PLS_cf  = function(f, d, s) Block.CCI.test(f, d, method = "PLS", nperm = s$nperm, crossfit = TRUE, progress = FALSE),
  CCI_SOPLS   = function(f, d, s) Block.CCI.test(f, d, method = "SOPLS", nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE),
  CCI_ROSA    = function(f, d, s) Block.CCI.test(f, d, method = "ROSA", nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE),
  FAST_RF     = function(f, d, s) Block.CCI.fast(f, d, method = "RF", progress = FALSE,
                                                 method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  RES_PLS     = function(f, d, s) Block.resid.test(f, d, method = "PLS", B = s$B, nfolds = s$nfolds),
  RES_RF      = function(f, d, s) Block.resid.test(f, d, method = "RF", B = s$B, nfolds = s$nfolds,
                                                   method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  RES_RF_wgcm = function(f, d, s) Block.resid.test(f, d, method = "RF", stat = "wgcm", k0 = s$k0, B = s$B,
                                                   nfolds = s$nfolds,
                                                   method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  CCI_KNN     = function(f, d, s) Block.CCI.test(f, d, method = "KNN", nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE),
  RES_KNN     = function(f, d, s) Block.resid.test(f, d, method = "KNN", B = s$B, nfolds = s$nfolds),
  RES_KNN_wgcm = function(f, d, s) Block.resid.test(f, d, method = "KNN", stat = "wgcm", k0 = s$k0, B = s$B,
                                                   nfolds = s$nfolds),
  # with dimension reduction of Z (for the high-dimensional designs)
  RES_PLS_pca = function(f, d, s) Block.resid.test(f, d, method = "PLS", reduce = "pca", B = s$B, nfolds = s$nfolds),
  RES_RF_pls  = function(f, d, s) Block.resid.test(f, d, method = "RF", reduce = "pls", B = s$B, nfolds = s$nfolds,
                                                   method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  RES_RF_pca  = function(f, d, s) Block.resid.test(f, d, method = "RF", reduce = "pca", B = s$B, nfolds = s$nfolds,
                                                   method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  RES_RF_pca_wgcm = function(f, d, s) Block.resid.test(f, d, method = "RF", reduce = "pca", stat = "wgcm", k0 = s$k0,
                                                   B = s$B, nfolds = s$nfolds,
                                                   method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  CCI_PLS_pca = function(f, d, s) Block.CCI.test(f, d, method = "PLS", reduce = "pca", nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE),
  CCI_RF_pca  = function(f, d, s) Block.CCI.test(f, d, method = "RF", reduce = "pca", nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE,
                                                 method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  CCI_RF_pls  = function(f, d, s) Block.CCI.test(f, d, method = "RF", reduce = "pls", nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE,
                                                 method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  # resampling of X for the null distribution: permutation (unconditional: breaks X-Z) and
  # conditional (greedy local kNN) permutation; the entries above use rotation
  CCI_PLS_perm = function(f, d, s) Block.CCI.test(f, d, method = "PLS", resampling = "permutation",
                                                  nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE),
  CCI_PLS_cond = function(f, d, s) Block.CCI.test(f, d, method = "PLS", resampling = "conditional",
                                                  k_perm = s$k_perm, nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE),
  CCI_RF_perm  = function(f, d, s) Block.CCI.test(f, d, method = "RF", resampling = "permutation",
                                                  nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE,
                                                  method_args = list(num.trees = s$num_trees, num.threads = s$threads)),
  CCI_RF_cond  = function(f, d, s) Block.CCI.test(f, d, method = "RF", resampling = "conditional",
                                                  k_perm = s$k_perm, nperm = s$nperm, p = s$p, npval = s$npval, progress = FALSE,
                                                  method_args = list(num.trees = s$num_trees, num.threads = s$threads))
)

# ---- noise and data generators ---------------------------------------------------------

# Unit-variance noise matrix (n x q). `s` is a standardised function of Z (mean of s^2
# about 1) that sets the noise level per row for type "hetero".
make_noise <- function(n, q, type, s = NULL) {
  m <- n * q
  e <- switch(type,
    gaussian = rnorm(m),
    t3       = rt(m, df = 3) / sqrt(3),                                      # heavy tails
    laplace  = { u <- runif(m) - 0.5; -sign(u) * log(1 - 2 * abs(u)) / sqrt(2) },
    skewed   = rexp(m) - 1,                                                  # right-skewed
    poisson  = rpois(m, 1) - 1,                                              # discrete counts
    outliers = ifelse(runif(m) < 0.05, rnorm(m, sd = 5), rnorm(m)) / sqrt(0.95 + 0.05 * 25),
    hetero   = rnorm(m) * sqrt(0.2 + 0.8 * s^2),   # s recycled down every column: row i gets s_i
    corr     = rnorm(m),
    stop("unknown noise type ", type))
  E <- matrix(e, n, q)
  if (type == "corr") E <- E %*% chol(0.6^abs(outer(seq_len(q), seq_len(q), "-")))   # AR(1) columns
  E
}

# Z-part of a block with q columns from features Fm (roughly unit variance): random
# coefficients scaled so that each column has variance about v.
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

# Blocks Y, X, Z1, Z2.
#   design "linear"    : Z1, Z2 with pz columns each; X and Y (q columns) linear in Z
#   design "nonlinear" : plus nonlinear functions of Z (half of the Z-part variance)
#   design "highdim"   : Z1, Z2 with 100 columns each from a 3-factor model; X and Y
#                        (10 columns) depend linearly and nonlinearly on the 3 factor
#                        directions of the observed Z. A deliberate LIMIT CASE: X and Y
#                        share nonlinear functions of a few directions hidden in 200
#                        noisy variables, which a regression must find from n rows
#                        (Shah and Peters 2020: no test can hold the level uniformly
#                        when the regressions cannot be learned).
#   design "highdim_linear": as "highdim", but X and Y depend only linearly on the
#                        factor directions (a linear method should hold the level)
#   snr                : variance of the Z-part of each X / Y column (noise variance 1)
gen_blocks <- function(n, noise = "gaussian", effect = 0, design = "linear",
                       snr = 2, pz = 4, q = 4) {
  design <- match.arg(design, c("linear", "nonlinear", "highdim", "highdim_linear"))
  if (design %in% c("highdim", "highdim_linear")) {
    pz <- 100; q <- 10; r <- 3
    L  <- matrix(rnorm(r * 2 * pz), r, 2 * pz)                      # factor loadings
    Z  <- matrix(rnorm(n * r), n, r) %*% L + matrix(rnorm(n * 2 * pz), n, 2 * pz)
    U  <- t(L) / sqrt(rowSums(L^2))[col(t(L))]                      # unit factor directions
    S  <- scale(Z %*% U)                                            # directions of the OBSERVED Z
    lin <- S; nl <- if (design == "highdim") nl_features(S) else NULL
  } else {
    Z   <- matrix(rnorm(n * 2 * pz), n, 2 * pz)
    lin <- Z
    nl  <- if (design == "nonlinear") nl_features(Z) else NULL
  }
  mu <- function()                                                  # same features, own coefficients
    if (is.null(nl)) mean_part(lin, q, snr) else mean_part(lin, q, snr / 2) + mean_part(nl, q, snr / 2)
  s <- as.numeric(scale(Z %*% rnorm(ncol(Z))))                      # noise level for "hetero"
  X <- mu() + make_noise(n, q, noise, s)
  Y <- mu() + make_noise(n, q, noise, s)
  if (effect != 0) Y <- Y + effect * X                              # direct X -> Y effect: H0 false
  list(Y = Y, X = X, Z1 = Z[, seq_len(pz), drop = FALSE], Z2 = Z[, pz + seq_len(pz), drop = FALSE])
}

# ---- progress bar with speed and time left ---------------------------------------------

.fmt_time <- function(s) {
  if (!is.finite(s)) return("--:--:--")
  s <- round(s)
  sprintf("%02d:%02d:%02d", s %/% 3600, (s %% 3600) %/% 60, s %% 60)
}

progress_start <- function(total, label) {
  pb <- new.env()
  pb$total <- total; pb$label <- label; pb$start <- proc.time()[["elapsed"]]
  cat(sprintf("%s: %d tests\n", label, total))
  pb
}

progress_tick <- function(pb, k, last = "") {
  el <- proc.time()[["elapsed"]] - pb$start
  rate <- k / el                                                    # tests per second
  width <- 25
  done <- round(width * k / pb$total)
  cat(sprintf("\r[%s%s] %d/%d | %.1f tests/min | elapsed %s | left ~%s | %s      ",
              strrep("=", done), strrep(" ", width - done), k, pb$total, 60 * rate,
              .fmt_time(el), .fmt_time((pb$total - k) / rate), last))
  if (k == pb$total) cat("\n")
  utils::flush.console()
}

# ---- the study runner (shared by all study functions) -----------------------------------

# Rejection rates per method (Monte Carlo s.e. = sqrt(rate * (1 - rate) / tests)).
summarise_tests <- function(res, alpha = 0.05) {
  out <- do.call(rbind, lapply(split(res, res$method), function(s) {
    ok <- !is.na(s$p); k <- sum(ok); r <- mean(s$p[ok] < alpha)
    data.frame(study = s$study[1], design = s$design[1], method = s$method[1], n = s$n[1],
               effect = s$effect[1], tests = k, rejection_rate = r,
               mc_se = sqrt(r * (1 - r) / max(k, 1)),
               rejection_rate_par = if (all(is.na(s$p_par))) NA_real_ else mean(s$p_par[ok] < alpha, na.rm = TRUE),
               median_sec = stats::median(s$seconds), errors = sum(!ok))
  }))
  rownames(out) <- NULL
  out
}

run_study <- function(study, generator, testing_method, effect, n, n_tests,
                      design = "linear", alpha = 0.05, settings = list(),
                      seed = NULL, save_dir = NULL, verbose = TRUE) {
  s <- utils::modifyList(default_settings, settings)
  if (is.function(testing_method)) {
    methods <- list(custom = testing_method)
  } else {
    bad <- setdiff(testing_method, names(test_methods))
    if (length(bad)) stop("unknown testing_method: ", paste(bad, collapse = ", "),
                          "; available: ", paste(names(test_methods), collapse = ", "))
    methods <- test_methods[testing_method]
  }
  if (!is.null(seed)) set.seed(seed)
  # One seed per dataset and one per test, all drawn up front. The generator and the tests
  # then do not share a random number stream, so dataset i depends only on `seed` (not on
  # which methods are run) and a method's result does not depend on the other methods or
  # their order. Filled column by column: a larger n_tests keeps the earlier seeds.
  seeds <- matrix(sample.int(.Machine$integer.max, 2 * n_tests), nrow = 2,
                  dimnames = list(c("data", "test"), NULL))
  total <- n_tests * length(methods)
  pb <- progress_start(total, sprintf("%s (%s), n = %d, effect = %g, %s", study, design, n, effect,
                                      paste(names(methods), collapse = " + ")))
  rows <- vector("list", total); k <- 0
  for (i in seq_len(n_tests)) {
    set.seed(seeds["data", i])
    d <- generator(n, effect)
    f <- if (is.null(attr(d, "formula"))) Y ~ X | Z1 + Z2 else attr(d, "formula")
    for (m in names(methods)) {
      set.seed(seeds["test", i])             # not the data seed: draws independent of the data
      t0 <- proc.time()[["elapsed"]]
      r <- tryCatch(methods[[m]](f, d, s), error = function(e) e)
      secs <- proc.time()[["elapsed"]] - t0
      failed <- inherits(r, "error")
      k <- k + 1
      rows[[k]] <- data.frame(
        study = study, design = design, method = m, n = n, effect = effect, test = i,
        data_seed = seeds[["data", i]],
        p = if (failed) NA_real_ else r$p.value,
        p_par = if (failed || is.null(r$p.parametric)) NA_real_ else r$p.parametric,
        seconds = secs, error = if (failed) conditionMessage(r) else NA_character_)
      progress_tick(pb, k, sprintf("last: %s %.1f s", m, secs))
    }
  }
  res <- do.call(rbind, rows)
  summ <- summarise_tests(res, alpha)
  out <- list(results = res, summary = summ,
              settings = c(s, list(study = study, design = design, n = n, effect = effect,
                                   n_tests = n_tests, alpha = alpha, seed = seed,
                                   seed_scheme = "per_dataset")))   # absent in older files
  if (!is.null(save_dir)) {
    save_dir <- results_path(save_dir)
    dir.create(save_dir,showWarnings = FALSE, recursive = TRUE)
    file <- file.path(save_dir, sprintf("%s_%s_%s_n%d_effect%s.rds", study, design,
                                        paste(names(methods), collapse = "-"), n, format(effect)))
    saveRDS(out, file)
  }
  if (verbose) {
    print(summ[, c("method", "n", "effect", "tests", "rejection_rate", "mc_se",
                   "rejection_rate_par", "median_sec", "errors")], digits = 3, row.names = FALSE)
    if (any(summ$errors > 0)) cat("First error:", res$error[which(!is.na(res$error))[1]], "\n")
  }
  invisible(out)
}

# ---- the studies (one data-generating process each) -------------------------------------

# Gaussian noise: the reference case.
study_gaussian_noise <- function(testing_method, effect = 0, n = 500, n_tests = 100, design = "linear", ...)
  run_study("gaussian_noise", function(n, effect) gen_blocks(n, "gaussian", effect, design),
            testing_method, effect, n, n_tests, design = design, ...)

# Heavy-tailed noise: Student t with 3 degrees of freedom.
study_t3_noise <- function(testing_method, effect = 0, n = 500, n_tests = 100, design = "linear", ...)
  run_study("t3_noise", function(n, effect) gen_blocks(n, "t3", effect, design),
            testing_method, effect, n, n_tests, design = design, ...)

# Laplace noise: moderately heavy tails.
study_laplace_noise <- function(testing_method, effect = 0, n = 500, n_tests = 100, design = "linear", ...)
  run_study("laplace_noise", function(n, effect) gen_blocks(n, "laplace", effect, design),
            testing_method, effect, n, n_tests, design = design, ...)

# Skewed noise: centred exponential.
study_skewed_noise <- function(testing_method, effect = 0, n = 500, n_tests = 100, design = "linear", ...)
  run_study("skewed_noise", function(n, effect) gen_blocks(n, "skewed", effect, design),
            testing_method, effect, n, n_tests, design = design, ...)

# Discrete noise: centred Poisson counts (lambda = 1, so skewed with many ties).
study_poisson_noise <- function(testing_method, effect = 0, n = 500, n_tests = 100, design = "linear", ...)
  run_study("poisson_noise", function(n, effect) gen_blocks(n, "poisson", effect, design),
            testing_method, effect, n, n_tests, design = design, ...)

# Outliers: 5 % of the noise values come from a normal with sd 5.
study_outlier_noise <- function(testing_method, effect = 0, n = 500, n_tests = 100, design = "linear", ...)
  run_study("outlier_noise", function(n, effect) gen_blocks(n, "outliers", effect, design),
            testing_method, effect, n, n_tests, design = design, ...)

# Heteroskedastic noise: the noise level depends on Z (the null still holds).
study_hetero_noise <- function(testing_method, effect = 0, n = 500, n_tests = 100, design = "linear", ...)
  run_study("hetero_noise", function(n, effect) gen_blocks(n, "hetero", effect, design),
            testing_method, effect, n, n_tests, design = design, ...)

# Correlated noise: the columns of the X and the Y noise are AR(1)-correlated (rho = 0.6).
study_correlated_noise <- function(testing_method, effect = 0, n = 500, n_tests = 100, design = "linear", ...)
  run_study("correlated_noise", function(n, effect) gen_blocks(n, "corr", effect, design),
            testing_method, effect, n, n_tests, design = design, ...)

# High-dimensional Z (200 variables in Z1, Z2) with a factor structure; X and Y have 10
# columns. `noise` chooses the noise type (default gaussian). nonlinear = TRUE (default)
# is the limit case "highdim"; nonlinear = FALSE the linear check "highdim_linear".
study_highdim <- function(testing_method, effect = 0, n = 500, n_tests = 100, noise = "gaussian",
                          nonlinear = TRUE, ...) {
  design <- if (nonlinear) "highdim" else "highdim_linear"
  run_study(paste0("highdim_", noise), function(n, effect) gen_blocks(n, noise, effect, design),
            testing_method, effect, n, n_tests, design = design, ...)
}

study_highdim_linear <- function(testing_method, effect = 0, n = 500, n_tests = 100, ...)
  study_highdim(testing_method, effect, n, n_tests, nonlinear = FALSE, ...)

# ---- reading saved studies ---------------------------------------------------------------

# Summary rows of all studies saved in `dir` (e.g. "t3_cci_results"), as one table.
collect_studies <- function(dir) {
  dir <- results_path(dir)
  files <-list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  out <- do.call(rbind, lapply(files, function(f) readRDS(f)$summary))
  out[order(out$study, out$design, out$method, out$effect, out$n), ]
}
# Single test
for (n in c(100, 200)) for (e in c(0, 0.5))
  study_gaussian_noise(c("RES_PLS"), effect = e, n = n, n_tests = 50, seed = n + 10 * e, save_dir = "test", design = "nonlinear")

# ---- Studies Residual based tests ------------------------------------------------------------------------------
# Linear cases
# for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
#   study_gaussian_noise(c("RES_PLS", "RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "gaussian_residual_results")
# 
# for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
#   study_t3_noise(c("RES_PLS", "RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "t3_residual_results")
# 
# for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
#   study_laplace_noise(c("RES_PLS", "RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "laplace_residual_results")
# 
# for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
#   study_skewed_noise(c("RES_PLS", "RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "skewed_residual_results")
# 
# for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
#   study_poisson_noise(c("RES_PLS", "RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "poisson_residual_results")
# 
# for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
#   study_outlier_noise(c("RES_PLS", "RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "outlier_residual_results")
# 
# for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
#   study_hetero_noise(c("RES_PLS", "RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "hetero_residual_results")
# 
# for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
#   study_correlated_noise(c("RES_PLS", "RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "correlated_residual_results")

# Non linear cases (Remove PLS based method since it is only designed for linear cases)
# for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
#   study_gaussian_noise(c("RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, design = "nonlinear", save_dir = "gaussian_residual_nonlin_results")

for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
  study_t3_noise(c("RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, design = "nonlinear", save_dir = "t3_residual_nonlin_results")

for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
  study_laplace_noise(c("RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, design = "nonlinear", save_dir = "laplace_residual_nonlin_results")

for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
  study_skewed_noise(c("RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, design = "nonlinear", save_dir = "skewed_residual_nonlin_results")

for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
  study_poisson_noise(c("RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, design = "nonlinear", save_dir = "poisson_residual_nonlin_results")

for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
  study_outlier_noise(c("RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, design = "nonlinear", save_dir = "outlier_residual_nonlin_results")

for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
  study_hetero_noise(c("RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, design = "nonlinear", save_dir = "hetero_residual_nonlin_results")

for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
  study_correlated_noise(c("RES_RF", "RES_RF_wgcm"), effect = e, n = n, n_tests = 100, seed = n + 10 * e, design = "nonlinear", save_dir = "correlated_residual_nonlin_results")

# ---- Studies CCI based tests ------------------------------------------------------------------------------
cci_methods <- c("CCI_PLS", "CCI_ROSA", "CCI_RF")


  for (n in c(200, 500, 1000)) 
    for (e in c(0, 0.5))
      study_gaussian_noise(cci_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "gaussian_cci_results")
  
  for (n in c(200, 500, 1000)) 
    for (e in c(0, 0.5))
      study_t3_noise(cci_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "t3_cci_results")
  
  for (n in c(200, 500, 1000)) 
    for (e in c(0, 0.5))
      study_laplace_noise(cci_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "laplace_cci_results")
  
  for (n in c(1000)) 
    for (e in c(0.5))
      study_skewed_noise(cci_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "skewed_cci_results")
  
  for (n in c(200, 500, 1000)) 
    for (e in c(0, 0.5))
      study_poisson_noise(cci_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "poisson_cci_results")
  
  for (n in c(200, 500, 1000)) 
    for (e in c(0, 0.5))
      study_outlier_noise(cci_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "outlier_cci_results")
  
  for (n in c(200, 500, 1000)) 
    for (e in c(0, 0.5))
      study_hetero_noise(cci_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "hetero_cci_results")
  
  for (n in c(200, 500, 1000)) 
    for (e in c(0, 0.5))
      study_correlated_noise(cci_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e, save_dir = "correlated_cci_results")
  
cci_methods <- c("CCI_RF")

studies <- list(
  gaussian   = study_gaussian_noise,
  t3         = study_t3_noise,
  laplace    = study_laplace_noise,
  skewed     = study_skewed_noise,
  poisson    = study_poisson_noise,
  outlier    = study_outlier_noise,
  hetero     = study_hetero_noise,
  correlated = study_correlated_noise
)

for (s in names(studies))
  for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
    studies[[s]](cci_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e, design = "nonlinear",
                 save_dir = paste0(s, "_nonlin_cci_results"))



# ---- Studies High Dimensions both residual and CCI  ------------------------------------------------------------------------------
study_highdim_linear("RES_PLS", effect = 0, n = 400, n_tests = 100)
study_highdim_linear("RES_PLS", effect = 0.5, n = 400, n_tests = 100)

study_highdim_linear("RES_RF_wgcm", effect = 0, n = 400, n_tests = 100)
study_highdim_linear("RES_PLS", effect = 0.5, n = 400, n_tests = 100)

study_highdim("RES_PLS", effect = 0, n = 400, n_tests = 10)
study_highdim("CCI_PLS", effect = 0, n = 400, n_tests = 5)



# ---- High-dimensional studies: linear check and the nonlinear limit case --------------------
# Linear high-dimensional design: linear methods should hold the level.
for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
   study_highdim_linear(c("RES_PLS", "RES_PLS_pca", "CCI_PLS", "CCI_PLS_pca"), effect = e, n = n,
                        n_tests = 100, seed = n + 10 * e, save_dir = "highdim_results_linear")
# Nonlinear limit case: tree learners on reduced Z (expected to stay liberal at large n).
for (n in c(200, 500, 1000)) for (e in c(0, 0.5))
  study_highdim(c("RES_RF_pls", "RES_RF_pca", "RES_RF_pca_wgcm", "CCI_RF_pca", "CCI_RF_pls"), effect = e, n = n,
                       n_tests = 100, seed = n + 10 * e, save_dir = "highdim_results_nonlinear")

# ---- Resampling of X: rotation vs permutation vs conditional permutation ---------------------
# Same seeds as the rotation studies above (seed = n + 10 * e). Studies saved with
# seed_scheme = "per_dataset" and the same seed, noise, design, n and effect use the same
# datasets and can be compared pair by pair. Results saved before that scheme (no
# seed_scheme in `settings`) share only the first dataset: compare their rates as
# independent samples. Hypothesis to test: conditional permutation controls the
# type I error less well than unconditional permutation in some settings.
# Priority: effect = 0 (type I error) first; then effect = 0.5 (power).
resampling_methods <- c("CCI_PLS_perm", "CCI_PLS_cond", "CCI_RF_perm", "CCI_RF_cond")
resampling_noise <- list(skewed = study_skewed_noise, outlier = study_outlier_noise, hetero = study_hetero_noise)

for (des in c("linear", "nonlinear"))
  for (s in names(resampling_noise))
    for (n in c(500, 1000)) for (e in c(0))
      resampling_noise[[s]](resampling_methods, effect = e, n = n, n_tests = 100, seed = n + 10 * e,
                            design = des,
                            save_dir = paste0(s, if (des == "nonlinear") "_nonlin", "_cci_resampling_results"))

# conditional permutation with a larger neighbourhood (default k_perm = 5)
 for (s in c("gaussian", "outlier")) for (n in c(500, 1000))
   resampling_noise[[s]](c("CCI_PLS_cond", "CCI_RF_cond"), effect = 0, n = n, n_tests = 100,
                         seed = n, design = "nonlinear", settings = list(k_perm = 10),
                         save_dir = paste0(s, "_nonlin_cci_cond_k10_results"))

