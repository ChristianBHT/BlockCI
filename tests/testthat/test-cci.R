# Block.CCI.test, cci_test_gen, cci_null_dist and Block.CCI.fast
library(testthat)
# OLS learner with the CCI learner interface (for the custom-method test)
ols_learner <- function(Ys, X, Vl, train, test, opts) {
  fitp <- function(F) {
    F1 <- cbind(1, F)
    F1[test, , drop = FALSE] %*% qr.solve(F1[train, , drop = FALSE], Ys[train, , drop = FALSE])
  }
  V <- if (is.null(Vl)) NULL else do.call(cbind, Vl)
  yh0 <- if (is.null(V)) matrix(0, length(test), ncol(Ys)) else fitp(V)
  list(yh0 = yh0, yh1 = fitp(cbind(V, X)), a = if (is.null(V)) 0 else ncol(V), b = ncol(cbind(V, X)))
}

test_that("Block.CCI.test (SOPLS) detects dependence and not independence", {
  set.seed(3); n <- 200
  V <- matrix(rnorm(n * 4), n); X <- V[, 1:2] + matrix(rnorm(n * 2), n)
  Y0 <- cbind(V[, 1] + rnorm(n), rnorm(n))
  Y1 <- Y0; Y1[, 1] <- Y1[, 1] + 1.5 * (X[, 1] - V[, 1])
  r0 <- Block.CCI.test(Y ~ X | V, list(Y = Y0, X = X, V = V), nperm = 39, ncomp = c(4, 2), progress = FALSE, seed = 1)
  r1 <- Block.CCI.test(Y ~ X | V, list(Y = Y1, X = X, V = V), nperm = 39, ncomp = c(4, 2), progress = FALSE, seed = 1)
  expect_s3_class(r1, "BlockCCI")
  expect_s3_class(r1, "htest")
  expect_gt(r0$p.value, 0.05)
  expect_lt(r1$p.value, 0.05)
})

test_that("Block.CCI.test returns the documented fields", {
  set.seed(4)
  d <- sim_1(n = 150, p = 3)
  r <- Block.CCI.test(Y ~ X | Z1 + Z2, d, method = "PLS", nperm = 19, progress = FALSE, seed = 1)
  expect_named(r$complexity, c("M0", "M1"))
  expect_length(r$null.distribution, 19)
  expect_equal(unname(r$statistic), r$M1 - r$M0)
  expect_equal(r$p.value, r$p.empirical)
  expect_equal(r$p.empirical, (sum(r$null.distribution <= r$statistic) + 1) / 20)
  expect_equal(r$method, "PLS")
  expect_equal(r$reduce, "none")
})

test_that("Block.CCI.test with PLS and ROSA detects a false null (Z2 missing)", {
  set.seed(8)
  d <- sim_1(n = 300, p = 3)
  for (m in c("PLS", "ROSA")) {
    r <- Block.CCI.test(Y ~ X | Z1, d, method = m, nperm = 39, progress = FALSE, seed = 1)
    expect_lt(r$p.value, 0.05)
  }
})

test_that("resampling options, parametric p-value and the unconditional test run", {
  set.seed(5); n <- 100
  d <- sim_1(n = n, p = 2)
  for (rs in c("permutation", "conditional")) {
    r <- Block.CCI.test(Y ~ X | Z1 + Z2, d, method = "PLS", resampling = rs, nperm = 9,
                        progress = FALSE, seed = 1)
    expect_equal(r$resampling, rs)
    expect_true(r$p.value > 0 && r$p.value <= 1)
  }
  rp <- Block.CCI.test(Y ~ X | Z1 + Z2, d, method = "PLS", pvalue = "parametric", nperm = 19,
                       progress = FALSE, seed = 1)
  expect_equal(rp$p.value, rp$p.parametric)
  X <- matrix(rnorm(n * 3), n); Y <- matrix(rnorm(n * 2), n)
  ru <- Block.CCI.test(Y ~ X, list(Y = Y, X = X), nperm = 9, resampling = "permutation", progress = FALSE)
  expect_true(ru$p.value > 0 && ru$p.value <= 1)
})

test_that("Block.CCI.test accepts a learner function", {
  set.seed(6)
  d <- sim_1(n = 200, p = 2)
  r <- Block.CCI.test(Y ~ X | Z1, d, method = ols_learner, nperm = 39, progress = FALSE, seed = 1)
  expect_equal(r$method, "custom")
  expect_lt(r$p.value, 0.05)
  expect_error(Block.CCI.test(Y ~ X | Z1, d, method = "nope", progress = FALSE), "unknown method")
})

test_that("Block.CCI.test runs with RF and XGB when installed", {
  set.seed(7)
  d <- sim_1(n = 150, p = 2)
  skip_if_not_installed("grf")
  r <- Block.CCI.test(Y ~ X | Z1 + Z2, d, method = "RF", nperm = 9, progress = FALSE, seed = 1,
                      method_args = list(num.trees = 100))
  expect_equal(r$method, "RF")
  expect_true(r$p.value > 0 && r$p.value <= 1)
  skip_if_not_installed("xgboost")
  r <- Block.CCI.test(Y ~ X | Z1 + Z2, d, method = "XGB", nperm = 9, progress = FALSE, seed = 1,
                      method_args = list(nrounds = 20))
  expect_equal(r$method, "XGB")
  expect_true(r$p.value > 0 && r$p.value <= 1)
})

test_that("print and plot methods for BlockCCI work", {
  set.seed(8)
  d <- sim_1(n = 120, p = 2)
  r <- Block.CCI.test(Y ~ X | Z1 + Z2, d, method = "PLS", nperm = 9, progress = FALSE, seed = 1)
  expect_output(print(r), "Block CCI test")
  rr <- Block.CCI.test(Y ~ X | Z1 + Z2, d, method = "PLS", nperm = 9, progress = FALSE, seed = 1,
                       reduce = "pca")
  expect_output(print(rr), "reduced by PCA")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_error(plot(r), NA)
})

test_that("cci_test_gen and cci_null_dist return the statistic and the null draws", {
  set.seed(9)
  d <- sim_1(n = 150, p = 3)
  g <- cci_test_gen(Y ~ X | Z1 + Z2, d, method = "PLS")
  expect_named(g, c("M0", "M1", "T", "a", "b"))
  expect_equal(g[["T"]], g[["M1"]] - g[["M0"]])
  T0 <- cci_null_dist(Y ~ X | Z1 + Z2, d, method = "PLS", nperm = 5)
  expect_length(T0, 5)
  expect_true(all(is.finite(T0)))
})

test_that("cci_test_gen reuses M0 predictions passed in yh0", {
  set.seed(9)
  d <- sim_1(n = 150, p = 3)
  s <- sample.int(150, 100)
  for (m in c("PLS", "ROSA")) {
    r1 <- cci_test_gen(Y ~ X | Z1 + Z2, d, method = m, split = s, losses = TRUE)
    d2 <- d; d2$X <- permute_X(d$X)
    r2 <- cci_test_gen(Y ~ X | Z1 + Z2, d2, method = m, split = s, losses = TRUE, yh0 = r1$yh0)
    expect_equal(r2$stats[["M0"]], r1$stats[["M0"]])
    expect_equal(r2$l0, r1$l0)
    expect_true(is.na(r2$stats[["a"]]))
  }
  expect_error(cci_test_gen(Y ~ X | Z1 + Z2, d, method = "PLS", split = s[-1], yh0 = r1$yh0),
               "does not match")
})

test_that("Block.CCI.fast runs with SOPLS, PLS, ROSA and RF", {
  set.seed(7)
  d <- sim_1(n = 200, p = 3)
  for (m in c("SOPLS", "PLS", "ROSA")) {
    r <- Block.CCI.fast(Y ~ X | Z1 + Z2, d, method = m, nsplit = 3, seed = 1, progress = FALSE)
    expect_s3_class(r, "BlockCCIfast")
    expect_equal(r$method, m)
    expect_length(r$delta, 3)
    expect_true(is.na(r$p.value) || (r$p.value >= 0 && r$p.value <= 1))
  }
  skip_if_not_installed("grf")
  r <- Block.CCI.fast(Y ~ X | Z1 + Z2, d, method = "RF", nsplit = 3,
                      method_args = list(num.trees = 100), seed = 1, progress = FALSE)
  expect_equal(r$method, "RF")
})

test_that("Block.CCI.fast: Nadeau-Bengio correction widens the se; false null is detected", {
  set.seed(10)
  d <- sim_1(n = 300, p = 3)
  r  <- Block.CCI.fast(Y ~ X | Z1, d, method = "PLS", nsplit = 10, nrot = 2, seed = 1, progress = FALSE)
  rc <- Block.CCI.fast(Y ~ X | Z1, d, method = "PLS", nsplit = 10, nrot = 2, seed = 1, progress = FALSE,
                       correction = "nadeau-bengio")
  expect_equal(r$delta, rc$delta)                      # same seed, same splits
  expect_gt(rc$se, r$se)
  expect_lt(r$p.value, 0.05)
  expect_output(print(r), "p-value")
  expect_error(Block.CCI.fast(Y ~ X | Z1, d, method = "PLS", nsplit = 1, progress = FALSE), "nsplit")
})

