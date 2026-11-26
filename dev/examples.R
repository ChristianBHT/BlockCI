# Playground: try the block conditional independence tests on different scenarios
# =============================================================================
# All tests test H0: Y is independent of X given the blocks after the bar.
#   p-value small -> reject H0 (evidence that X and Y are related beyond Z)
#   effect = 0    -> H0 is true: small p-values are false alarms
#   effect > 0    -> H0 is false (Y <- Y + effect * X): small p-values are what we want
#
# How to use: open BlockCI.Rproj, run section 0, then run any other section line by
# line (Ctrl+Enter). In RStudio the sections are listed in the outline (Ctrl+Shift+O).
# Sections 1-13 need only section 0. Slow parts are skipped unless run_slow = TRUE.
#
# One p-value says little about whether a test works. To judge a test, repeat it on
# many datasets (section 13).

# ---- 0. Setup (run this first) ----------------------------------------------

devtools::load_all()          # the current source code; no reinstall needed

run_slow <- FALSE             # TRUE: also run random forest, XGBoost and other slow parts
nperm    <- 99                # null draws in the CCI tests (smallest p-value 1/(nperm + 1))
B        <- 199               # permutations / Gaussian draws in the residual tests

has <- function(pkg) requireNamespace(pkg, quietly = TRUE)    # optional packages

# The data generators of the simulation study. Only these function definitions are
# read from the script; no study is started.
local({
  wanted <- c("make_noise", "mean_part", "nl_features", "gen_blocks")
  for (e in parse("inst/benchmarks/simulation_study_noise.R"))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && as.character(e[[2]])[1] %in% wanted)
      eval(e, globalenv())
})

# One dataset: list(Y, X, Z1, Z2).
#   design: "linear", "nonlinear", "highdim_linear" or "highdim" (200 Z columns, 10
#           columns in X and Y; "highdim" is nonlinear and a deliberately hard case)
#   noise : "gaussian", "t3", "laplace", "skewed", "poisson", "outliers", "hetero"
#           (noise level depends on Z) or "corr" (correlated columns)
# The same arguments and seed always give the same data.
scenario <- function(design = "linear", noise = "gaussian", effect = 0, n = 500, seed = 1) {
  set.seed(10000 + seed)      # kept apart from the seeds used by the tests
  gen_blocks(n, noise = noise, effect = effect, design = design)
}

# Run several tests on the same data. `tests` is a named list of functions of the
# data. Returns one row per test with the p-value and the time in seconds.
compare <- function(d, tests, seed = 1) {
  rows <- lapply(names(tests), function(name) {
    set.seed(seed)
    time <- system.time(r <- tryCatch(tests[[name]](d), error = function(e) e))[["elapsed"]]
    failed <- inherits(r, "error")
    if (failed) message(name, ": ", conditionMessage(r))
    data.frame(test = name, p = if (failed) NA else round(r$p.value, 4),
               seconds = round(time, 1))
  })
  do.call(rbind, rows)
}

# Short forms of the two main tests, used in the sections below
cci <- function(d, ..., formula = Y ~ X | Z1 + Z2)
  Block.CCI.test(formula, d, nperm = nperm, progress = FALSE, ...)
res <- function(d, ..., formula = Y ~ X | Z1 + Z2)
  Block.resid.test(formula, d, B = B, ...)


# ---- 1. Quick start ---------------------------------------------------------

dat <- sim_1()                # linear Gaussian, 6 columns per block, H0 true given Z1, Z2

Block.CCI.test(Y ~ X | Z1 + Z2, data = dat, nperm = nperm)       # H0 true
Block.resid.test(Y ~ X | Z1 + Z2, data = dat, B = B)

Block.CCI.test(Y ~ X | Z2, data = dat, nperm = nperm)            # Z1 left out: H0 false
Block.resid.test(Y ~ X | Z2, data = dat, B = B)


# ---- 2. Look inside one CCI test --------------------------------------------

d1 <- scenario("linear", effect = 0.15)
r  <- Block.CCI.test(Y ~ X | Z1 + Z2, d1, nperm = nperm, seed = 1)
r
plot(r)                       # null statistics and the observed one (reject if far left)

r$statistic                   # T = RMSEP(Y ~ Z + X) - RMSEP(Y ~ Z); negative: X helps
r$p.empirical
r$p.parametric                # normal approximation of the null distribution
r$complexity                  # here: PLS components without and with X
summary(r$null.distribution)

# The two building blocks: one statistic, and the null statistics
cci_test_gen(Y ~ X | Z1 + Z2, d1)
T0 <- cci_null_dist(Y ~ X | Z1 + Z2, d1, nperm = 20)
T0


# ---- 3. Pick a scenario -----------------------------------------------------
# Change the arguments and rerun the tests below.

d <- scenario(design = "linear", noise = "gaussian", effect = 0, n = 500, seed = 1)
str(d)

compare(d, list(
  "CCI PLS"             = function(d) cci(d),
  "residual PLS, RV"    = function(d) res(d),
  "residual PLS, wGCM"  = function(d) res(d, stat = "wgcm")
))


# ---- 4. Regression methods in the CCI test ----------------------------------

d <- scenario("linear", effect = 0.2)

methods_cci <- list(
  PLS   = function(d) cci(d, method = "PLS"),
  SOPLS = function(d) cci(d, method = "SOPLS"),     # one component per Z block by default
  ROSA  = function(d) cci(d, method = "ROSA")
)
if (has("FNN"))
  methods_cci$KNN <- function(d) cci(d, method = "KNN")
if (run_slow && has("grf"))
  methods_cci$RF <- function(d) cci(d, method = "RF",
                                    method_args = list(num.trees = 500, num.threads = NULL))
if (run_slow && has("xgboost"))
  methods_cci$XGB <- function(d) cci(d, method = "XGB", method_args = list(nrounds = 100))

compare(d, methods_cci)

# Method settings
cci(d, method = "PLS", ncomp = c(8, 10))            # fixed components without / with X
cci(d, method = "SOPLS", ncomp = c(4, 4, 2))        # Z1, Z2, then X
if (has("FNN")) cci(d, method = "KNN", method_args = list(k = 15))


# ---- 5. How X is resampled for the null distribution ------------------------
# "rotation" (default): keeps the linear relation between X and Z
# "permutation"       : shuffles all rows, which breaks the relation between X and Z
# "conditional"       : local kNN permutation (Runge 2018); k_perm = neighbours

d <- scenario("linear", effect = 0)                 # H0 true

compare(d, list(
  rotation               = function(d) cci(d, resampling = "rotation"),
  permutation            = function(d) cci(d, resampling = "permutation"),
  "conditional, k = 5"   = function(d) cci(d, resampling = "conditional"),
  "conditional, k = 10"  = function(d) cci(d, resampling = "conditional", k_perm = 10)
))

# See what the three schemes do: Z has one column, X depends on it nonlinearly
set.seed(2)
z <- matrix(rnorm(300))
x <- matrix(sin(2 * z) + rnorm(300, sd = 0.2))
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(z, x, main = "observed X")
plot(z, rotate_X(x, z), main = "rotation")          # only the linear part survives
plot(z, permute_X(x), main = "permutation")
plot(z, permute_X_cond(x, z, k_perm = 5), main = "conditional")
par(op)


# ---- 6. Training fraction, observed splits and cross-fitting ----------------

d <- scenario("linear", effect = 0.12)

compare(d, list(
  "single split, p = 0.5 (default)" = function(d) cci(d),
  "single split, p = 0.75"          = function(d) cci(d, p = 0.75),
  "cross-fitted (experimental)"     = function(d) cci(d, crossfit = TRUE),
  "npval = 5 (experimental)"        = function(d) cci(d, npval = 5)
))
# crossfit: promising with PLS on linear data; too many false alarms with RF.
# npval > 1: averaged p-values were far too cautious in simulations.

# How much does the p-value depend on the random numbers? Same data, 10 test seeds
p_over_seeds <- function(d, test, seeds = 1:10)
  sapply(seeds, function(s) { set.seed(s); test(d)$p.value })

summary(p_over_seeds(d, function(d) cci(d)))
summary(p_over_seeds(d, function(d) cci(d, crossfit = TRUE)))


# ---- 7. Residual tests ------------------------------------------------------
# Z is removed from X and from Y by regression; then: are the residuals related?

d <- scenario("linear", effect = 0.1)

# regression methods
methods_res <- list(
  PLS  = function(d) res(d, method = "PLS"),
  ROSA = function(d) res(d, method = "ROSA")
)
if (has("FNN"))
  methods_res$"PLS + KNN" <- function(d) res(d, method = "KNN")
if (run_slow && has("grf"))
  methods_res$"PLS + RF" <- function(d) res(d, method = "RF", nfolds = 5,
                                            method_args = list(num.trees = 500, num.threads = NULL))
if (run_slow && has("xgboost"))
  methods_res$"PLS + XGB" <- function(d) res(d, method = "XGB", nfolds = 5,
                                             method_args = list(nrounds = 100))
compare(d, methods_res)

# statistics and permutation schemes
stats_res <- list(
  "RV, conditional (default)" = function(d) res(d),
  "RV, simple permutation"    = function(d) res(d, perm = "simple"),
  "first canonical corr."     = function(d) res(d, stat = "cc1"),
  "weighted GCM"              = function(d) res(d, stat = "wgcm")
)
if (has("energy"))
  stats_res$"distance corr." <- function(d) res(d, stat = "dcor")
compare(d, stats_res)

out <- res(d)
out
hist(out$perm_stats, breaks = 30, main = "residual test: permutation null",
     xlab = "RV", xlim = range(out$perm_stats, out$statistic))
abline(v = out$statistic, col = "red", lwd = 2)

# KNN without the linear PLS step (gave many false alarms in simulations)
if (has("FNN")) res(scenario("linear"), method = "KNN", linear_first = FALSE)


# ---- 8. Noise that depends on Z (heteroskedastic) ---------------------------
# H0 is true. In simulations the RV residual test gave far too many false alarms here;
# the weighted GCM held the level much better.

dh <- scenario("linear", noise = "hetero", effect = 0)

compare(dh, list(
  "residual PLS, RV"   = function(d) res(d),
  "residual PLS, wGCM" = function(d) res(d, stat = "wgcm"),
  "CCI PLS"            = function(d) cci(d)
))


# ---- 9. Other noise distributions -------------------------------------------

noises <- c("gaussian", "t3", "laplace", "skewed", "poisson", "outliers", "hetero", "corr")

do.call(rbind, lapply(noises, function(noise) {
  d <- scenario("linear", noise = noise, effect = 0)
  cbind(noise = noise, compare(d, list(
    "CCI PLS"            = function(d) cci(d),
    "residual PLS, wGCM" = function(d) res(d, stat = "wgcm"))))
}))


# ---- 10. Nonlinear relations ------------------------------------------------
# H0 is true, but linear regressions cannot remove Z, so they raise false alarms.
# No method here has held the level reliably at n = 500 in all our simulations.

dn <- scenario("nonlinear", effect = 0)

nonlinear_tests <- list(
  "CCI PLS"          = function(d) cci(d, method = "PLS"),
  "residual PLS, RV" = function(d) res(d)
)
if (has("FNN"))
  nonlinear_tests$"CCI KNN" <- function(d) cci(d, method = "KNN")
if (run_slow && has("grf")) {
  nonlinear_tests$"CCI RF" <- function(d) cci(d, method = "RF",
                              method_args = list(num.trees = 500, num.threads = NULL))
  nonlinear_tests$"residual PLS + RF, wGCM" <- function(d) res(d, method = "RF", stat = "wgcm",
                              nfolds = 5, method_args = list(num.trees = 500, num.threads = NULL))
}
compare(dn, nonlinear_tests)


# ---- 11. Many conditioning variables ----------------------------------------
# 200 Z columns built from 3 factors. reduce = "pca" / "pls" shrinks Z first.

dz <- scenario("highdim_linear", effect = 0.2, n = 500)
sapply(dz, dim)

compare(dz, list(
  "residual PLS"        = function(d) res(d),
  "residual PLS on PCA" = function(d) res(d, reduce = "pca"),
  "CCI PLS"             = function(d) cci(d),
  "CCI PLS on PCA"      = function(d) cci(d, reduce = "pca"),
  "CCI PLS on PLS"      = function(d) cci(d, reduce = "pls")
))


# ---- 12. Other data generators and your own method --------------------------

# sim_2: common causes, sim_3: chain, sim_4: nearly collinear X, sim_5: nonlinear
# (in all of them H0 is true given Z1 and Z2)
sims <- list(sim_2 = sim_2, sim_3 = sim_3, sim_4 = sim_4, sim_5 = sim_5)
do.call(rbind, lapply(names(sims), function(name) {
  set.seed(10001)
  cbind(data = name, compare(sims[[name]](n = 500), list("CCI PLS" = function(d) cci(d),
                                                          "residual PLS" = function(d) res(d))))
}))

# sim_scale: choose rows, number of Z columns, noise and effect
ds <- sim_scale(n = 2000, p_z = 40, q_x = 6, q_y = 3, noise = "t3", effect = 0.05)
compare(ds, list("residual PLS" = function(d) res(d), "CCI PLS" = function(d) cci(d)))

# Your own data: a named list of matrices with the same rows, used like this:
#   mydata <- list(Y = Ymat, X = Xmat, A = Amat, B = Bmat)
#   Block.CCI.test(Y ~ X | A + B, mydata)

# Your own regression in the CCI test: predict the test rows without and with X
ols_cci <- function(Ys, X, Vl, train, test, opts) {
  V <- do.call(cbind, Vl)
  fit <- function(F) {
    F <- cbind(1, F)
    F[test, , drop = FALSE] %*% qr.solve(F[train, , drop = FALSE], Ys[train, , drop = FALSE])
  }
  list(yh0 = fit(V), yh1 = fit(cbind(V, X)), a = ncol(V), b = ncol(V) + ncol(X))
}

# ... and in the residual test: return cross-fitted residuals of M given Z
ols_res <- function(M, Zl, folds, opts) {
  Z <- cbind(1, do.call(cbind, Zl))
  out <- M
  for (k in unique(folds)) {
    te <- folds == k
    out[te, ] <- M[te, , drop = FALSE] -
      Z[te, , drop = FALSE] %*% qr.solve(Z[!te, , drop = FALSE], M[!te, , drop = FALSE])
  }
  list(res = out, complexity = ncol(Z))
}

d <- scenario("linear", effect = 0.15)
compare(d, list("CCI, own OLS"      = function(d) cci(d, method = ols_cci),
                "residual, own OLS" = function(d) res(d, method = ols_res)))


# ---- 13. Repeat over many datasets: false alarm rate and power --------------
# For each repetition a new dataset is drawn and the test is run on it. The share of
# p-values below alpha is the type I error (effect = 0) or the power (effect > 0).
# `se` is its Monte Carlo standard error: with 20 repetitions it is large (about 0.05
# at a rate of 0.05). Use several hundred repetitions before drawing conclusions.

rejection_rate <- function(tests, reps = 20, alpha = 0.05, ...) {
  p <- sapply(seq_len(reps), function(i) {
    d <- scenario(..., seed = i)             # ... : design, noise, effect, n
    out <- compare(d, tests, seed = i)$p
    cat("\r", i, "/", reps)
    out
  })
  cat("\r")
  p <- matrix(p, nrow = length(tests))
  rate <- rowMeans(p < alpha, na.rm = TRUE)
  ok <- rowSums(!is.na(p))
  data.frame(test = names(tests), rejection_rate = round(rate, 3),
             se = round(sqrt(rate * (1 - rate) / ok), 3), datasets = ok)
}

my_tests <- list(
  "CCI PLS"            = function(d) cci(d),
  "residual PLS, RV"   = function(d) res(d),
  "residual PLS, wGCM" = function(d) res(d, stat = "wgcm")
)

rejection_rate(my_tests, reps = 20, design = "linear", noise = "hetero", effect = 0)
rejection_rate(my_tests, reps = 20, design = "linear", noise = "gaussian", effect = 0.1)

if (run_slow) {
  # power curve over the effect size
  effects <- c(0, 0.05, 0.1, 0.15, 0.2)
  curve_res <- do.call(rbind, lapply(effects, function(e)
    cbind(effect = e, rejection_rate(my_tests, reps = 100, design = "linear", effect = e))))
  print(curve_res)
  tab <- with(curve_res, tapply(rejection_rate, list(effect, test), identity))
  matplot(effects, tab, type = "b", pch = 1:ncol(tab), lty = 1, col = 1:ncol(tab),
          ylim = c(0, 1), xlab = "effect", ylab = "share with p < 0.05")
  abline(h = 0.05, lty = 3)
  legend("topleft", colnames(tab), pch = 1:ncol(tab), col = 1:ncol(tab), lty = 1, bty = "n")
}

# The full simulation study, with saved results, is in
# inst/benchmarks/simulation_study_noise.R. Running that whole file starts the studies.
