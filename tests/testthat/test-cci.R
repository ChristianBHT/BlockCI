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
  expect_equal(r$npval, 1)                                          # the default
  expect_length(r$T.splits, 1)
  expect_equal(r$p.empirical, (sum(r$null.distribution <= r$statistic) + 1) / 20)
  expect_equal(r$method, "PLS")
  expect_equal(r$reduce, "none")
})

test_that("Block.CCI.test: npval sets the number of observed splits", {
  set.seed(4)
  d <- sim_1(n = 150, p = 3)
  run <- function(k) Block.CCI.test(Y ~ X | Z1 + Z2, d, method = "PLS", nperm = 9, npval = k,
                                    progress = FALSE, seed = 1)
  r1 <- run(1); r3 <- run(3); r9 <- run(9)                          # 9 = nperm is allowed
  expect_length(r3$T.splits, 3)
  expect_length(r9$p.empirical.splits, 9)
  # the reported values are averages over the splits
  expect_equal(unname(r3$statistic), mean(r3$T.splits))
  expect_equal(unname(r3$statistic), r3$M1 - r3$M0)
  expect_equal(r3$p.empirical.splits,
               sapply(r3$T.splits, function(t) (sum(r3$null.distribution <= t) + 1) / 10))
  expect_equal(r3$p.value, mean(r3$p.empirical.splits))
  expect_equal(r3$p.parametric, mean(r3$p.parametric.splits))
  # same seed: the null draws and the first observed split do not depend on npval
  expect_identical(r3$null.distribution, r1$null.distribution)
  expect_identical(r3$T.splits[1], r1$T.splits)
  # one split: the single-split p-value
  expect_equal(r1$p.value, (sum(r1$null.distribution <= r1$statistic) + 1) / 10)
  expect_output(print(r3), "average of 3 p-values")
  expect_false(any(grepl("average of", capture.output(print(r1)))))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_error(plot(r3), NA)
  expect_error(plot(r1), NA)
  for (bad in list(0, 10, 1.5, -1, NA, c(1, 2), "2")) expect_error(run(bad), "npval")
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


test_that("KNN regression: predictions are neighbour means; k is chosen by CV", {
  skip_if_not_installed("FNN")
  set.seed(20)
  n <- 120; Fm <- matrix(rnorm(n * 3), n); Y <- matrix(rnorm(n * 2), n)
  train <- 1:80; test <- 81:120
  f <- .knn_fit_predict(Fm, Y, train, test, k = 4, scale = FALSE)
  D <- as.matrix(stats::dist(Fm))[test, train]
  brute <- t(apply(D, 1, function(d) colMeans(Y[train[order(d)[1:4]], , drop = FALSE])))
  expect_equal(f$pred, unname(brute))
  expect_equal(f$k, 4L)
  # features are standardised with the training rows: the units of a column do not matter
  a <- .knn_fit_predict(Fm, Y, train, test, k = 4)
  b <- .knn_fit_predict(sweep(Fm, 2, c(1000, 1, 0.01), "*"), Y, train, test, k = 4)
  expect_equal(a$pred, b$pred)
  means <- .knn_means(Y, matrix(1:6, 2), c(1, 3))
  expect_equal(means[[1]], Y[1:2, ])
  expect_equal(means[[2]][1, ], colMeans(Y[c(1, 3, 5), ]))
  # pure noise favours many neighbours, a sharp signal few
  x <- matrix(runif(400), 400)
  expect_gte(.knn_choose_k(x, matrix(rnorm(400)), .knn_default_grid), 20)
  expect_lte(.knn_choose_k(x, matrix(sin(20 * x) + rnorm(400, sd = 0.05)), .knn_default_grid), 10)
  expect_lte(.knn_choose_k(x[1:6, , drop = FALSE], matrix(rnorm(6)), .knn_default_grid), 4)
})

test_that("Block.CCI.test, Block.CCI.fast and the building blocks run with KNN", {
  skip_if_not_installed("FNN")
  set.seed(21)
  n <- 300; Z <- matrix(rnorm(n * 2), n)
  X <- Z + matrix(rnorm(n * 2), n)
  d0 <- list(Y = Z + matrix(rnorm(n * 2), n), X = X, Z = Z)
  d1 <- d0; d1$Y <- d0$Y + 1.5 * X
  for (rs in c("rotation", "conditional")) {
    r1 <- Block.CCI.test(Y ~ X | Z, d1, method = "knn", resampling = rs, nperm = 39,
                         progress = FALSE, seed = 1)
    expect_equal(r1$method, "KNN")
    expect_lt(r1$p.value, 0.05)
  }
  r <- Block.CCI.test(Y ~ X | Z, d0, method = "KNN", nperm = 19, progress = FALSE, seed = 1,
                      method_args = list(k = 12))
  expect_equal(unname(r$complexity), c(12, 12))
  expect_output(print(r), "neighbours")
  g <- cci_test_gen(Y ~ X | Z, d0, method = "KNN", method_args = list(k_grid = c(4, 8)))
  expect_true(all(g[c("a", "b")] %in% c(4, 8)))
  expect_length(cci_null_dist(Y ~ X | Z, d0, method = "KNN", nperm = 3), 3)
  ru <- Block.CCI.test(Y ~ X, d1, method = "KNN", resampling = "permutation", nperm = 19,
                       progress = FALSE, seed = 1)
  expect_lt(ru$p.value, 0.1)
  f <- Block.CCI.fast(Y ~ X | Z, d1, method = "KNN", nsplit = 4, progress = FALSE, seed = 1)
  expect_equal(f$method, "KNN")
  expect_lt(f$p.value, 0.05)
  for (bad in list(list(k = 0), list(k = 2.5), list(neighbours = 5), list(k_grid = c(0, 5))))
    expect_error(cci_test_gen(Y ~ X | Z, d0, method = "KNN", method_args = bad), "KNN")
})

test_that("the cross-fitted statistic pools the errors of both halves", {
  set.seed(30)
  d <- sim_1(n = 101, p = 2)                                        # odd number of rows
  set.seed(1); cf <- .cci_crossfit(Y ~ X | Z1 + Z2, d, method = "PLS", ncomp = c(2, 3))
  expect_named(cf, c("M0", "M1", "T", "a", "b"))
  set.seed(1); half <- sample.int(101, 50)
  folds <- lapply(list(half, setdiff(1:101, half)), function(s)
    cci_test_gen(Y ~ X | Z1 + Z2, d, method = "PLS", ncomp = c(2, 3), split = s, losses = TRUE))
  expect_setequal(c(folds[[1]]$test, folds[[2]]$test), 1:101)       # every row predicted once
  l0 <- c(folds[[1]]$l0, folds[[2]]$l0); l1 <- c(folds[[1]]$l1, folds[[2]]$l1)
  expect_equal(cf[["M0"]], sqrt(sum(l0) / (101 * 2)))
  expect_equal(cf[["M1"]], sqrt(sum(l1) / (101 * 2)))
  expect_equal(cf[["T"]], cf[["M1"]] - cf[["M0"]])
  expect_equal(unname(cf[c("a", "b")]), c(2, 3))
  d$Y[1, 1] <- NA                                                   # incomplete row dropped first
  expect_true(is.finite(.cci_crossfit(Y ~ X | Z1 + Z2, d, method = "PLS")[["T"]]))
})

test_that("Block.CCI.test(crossfit = TRUE) cross-fits the observed statistic and the null", {
  set.seed(31)
  d <- sim_1(n = 200, p = 3)
  r <- Block.CCI.test(Y ~ X | Z1, d, nperm = 39, crossfit = TRUE, progress = FALSE, seed = 1)
  expect_true(r$crossfit)
  expect_lt(r$p.value, 0.05)                                        # Z2 missing: H0 false
  expect_equal(r$p.value, (sum(r$null.distribution <= r$statistic) + 1) / 40)
  expect_equal(unname(r$statistic), r$M1 - r$M0)
  # observed and null are the cross-fitted statistic, in this order of random numbers
  set.seed(1)
  obs <- .cci_crossfit(Y ~ X | Z1, d, method = "PLS")
  T0 <- cci_null_dist(Y ~ X | Z1, d, nperm = 39, crossfit = TRUE)
  expect_equal(unname(r$statistic), obs[["T"]])
  expect_equal(r$null.distribution, T0)
  expect_output(print(r), "cross-fitted")
  expect_output(print(r), "mean of the 2 folds")
  for (rs in c("permutation", "conditional")) {
    r2 <- Block.CCI.test(Y ~ X | Z1 + Z2, d, nperm = 9, crossfit = TRUE, resampling = rs,
                         progress = FALSE, seed = 1)
    expect_true(r2$p.value > 0 && r2$p.value <= 1)
  }
  # off by default, and the default test is unchanged by the new argument
  a <- Block.CCI.test(Y ~ X | Z1, d, nperm = 9, progress = FALSE, seed = 1)
  b <- Block.CCI.test(Y ~ X | Z1, d, nperm = 9, crossfit = FALSE, progress = FALSE, seed = 1)
  expect_false(a$crossfit)
  expect_identical(a$null.distribution, b$null.distribution)
  expect_false(any(grepl("cross-fitted", capture.output(print(a)))))
  expect_error(Block.CCI.test(Y ~ X | Z1, d, crossfit = TRUE, p = 0.75, progress = FALSE), "must be 0.5")
  expect_error(Block.CCI.test(Y ~ X | Z1, d, crossfit = TRUE, npval = 2, progress = FALSE), "npval")
  expect_error(Block.CCI.test(Y ~ X | Z1, d, crossfit = "yes", progress = FALSE), "TRUE or FALSE")
  expect_error(cci_null_dist(Y ~ X | Z1, d, nperm = 3, crossfit = TRUE, p = 0.75), "must be 0.5")
})
