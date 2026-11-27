test_that("Block.resid.test (PLS, ROSA) rejects a false null and returns an htest", {
  set.seed(5)
  d <- sim_1(n = 300, p = 3)
  for (m in c("PLS", "ROSA")) {
    r <- Block.resid.test(Y ~ X | Z1, d, method = m, B = 99, seed = 1)   # Z2 missing: H0 false
    expect_s3_class(r, "BlockResid")
    expect_s3_class(r, "htest")
    expect_lt(r$p.value, 0.05)
    expect_named(r$complexity, c("X", "Y"))
  }
})

test_that("Block.resid.test runs with fitted residuals, simple permutation and fixed ncomp", {
  set.seed(6)
  d <- sim_1(n = 150, p = 3)
  r <- Block.resid.test(Y ~ X | Z1 + Z2, d, resid = "fitted", perm = "simple",
                        ncomp = 3, B = 19, seed = 1)
  expect_true(r$p.value > 0 && r$p.value <= 1)
  expect_equal(unname(r$complexity), c(3L, 3L))
})

test_that("Block.resid.test works with XGB and RF (PLS step first) when installed", {
  d <- sim_1(n = 150, p = 3)
  skip_if_not_installed("xgboost")
  r <- Block.resid.test(Y ~ X | Z1 + Z2, d, method = "XGB", B = 19, nfolds = 3,
                        method_args = list(nrounds = 20), seed = 1)
  expect_true(r$p.value > 0 && r$p.value <= 1)
  expect_true(r$linear_first)
  expect_true(all(r$linear_ncomp >= 1))
  r0 <- Block.resid.test(Y ~ X | Z1 + Z2, d, method = "XGB", B = 19, nfolds = 3,
                         linear_first = FALSE, method_args = list(nrounds = 20), seed = 1)
  expect_false(r0$linear_first)
  expect_true(all(is.na(r0$linear_ncomp)))
  skip_if_not_installed("grf")
  r <- Block.resid.test(Y ~ X | Z1 + Z2, d, method = "RF", B = 19, nfolds = 3,
                        method_args = list(num.trees = 100), seed = 1)
  expect_true(r$p.value > 0 && r$p.value <= 1)
  expect_true(all(r$linear_ncomp >= 1))
  r0 <- Block.resid.test(Y ~ X | Z1 + Z2, d, method = "RF", B = 19, linear_first = FALSE,
                         method_args = list(num.trees = 100), seed = 1)   # out-of-bag
  expect_true(r0$p.value > 0 && r0$p.value <= 1)
})

test_that("Block.resid.test accepts a learner function and rejects unknown methods", {
  d <- sim_1(n = 100, p = 2)
  ols <- function(M, Zl, folds, opts) {
    Z <- cbind(1, do.call(cbind, Zl))
    list(res = M - Z %*% qr.solve(Z, M), complexity = ncol(Z))
  }
  r <- Block.resid.test(Y ~ X | Z1 + Z2, d, method = ols, B = 19, seed = 1)
  expect_equal(r$regression, "custom")
  expect_error(Block.resid.test(Y ~ X | Z1, d, method = "nope"), "unknown method")
})

test_that(".wgcm_weights builds k0 * ncol(Z) + 1 weights", {
  set.seed(10)
  Z <- matrix(rnorm(200 * 3), 200, 3)
  W <- .wgcm_weights(Z, k0 = 7)
  expect_equal(ncol(W), 7 * 3 + 1)
  expect_true(all(W[, 1] == 1))
  expect_true(all(W %in% c(-1, 0, 1)))
})

test_that("Block.resid.test with stat = 'wgcm' rejects a false null and runs under H0", {
  set.seed(11)
  d <- sim_1(n = 300, p = 2)
  r <- Block.resid.test(Y ~ X | Z1, d, stat = "wgcm", B = 199, seed = 1)   # Z2 missing: H0 false
  expect_s3_class(r, "BlockResid")
  expect_lt(r$p.value, 0.05)
  expect_equal(r$n_stats, 2 * 2 * (7 * 2 + 1))
  r0 <- Block.resid.test(Y ~ X | Z1 + Z2, d, stat = "wgcm", B = 99, k0 = 3, seed = 1)
  expect_true(r0$p.value > 0 && r0$p.value <= 1)
  expect_equal(r0$n_stats, 2 * 2 * (3 * 4 + 1))
})

test_that("Block.resid.test runs with the cc1 and dcor statistics", {
  set.seed(12)
  d <- sim_1(n = 200, p = 2)
  r <- Block.resid.test(Y ~ X | Z1, d, stat = "cc1", B = 99, seed = 1)   # Z2 missing: H0 false
  expect_equal(r$stat, "cc1")
  expect_lt(r$p.value, 0.05)
  skip_if_not_installed("energy")
  r <- Block.resid.test(Y ~ X | Z1 + Z2, d, stat = "dcor", B = 19, seed = 1)
  expect_true(r$p.value > 0 && r$p.value <= 1)
})

test_that("Block.resid.test: permutation statistics, print method and argument checks", {
  set.seed(13)
  d <- sim_1(n = 150, p = 2)
  r <- Block.resid.test(Y ~ X | Z1 + Z2, d, B = 49, seed = 1)
  expect_length(r$perm_stats, 49)
  expect_true(r$p.value > 0 && r$p.value <= 1)
  expect_output(print(r), "Block residual test")
  expect_error(Block.resid.test(Y ~ X, d), "conditioning set")
  expect_error(Block.resid.test(Y ~ X | Z1, d, nfolds = 1), "nfolds")
})

test_that(".right_tail_p is accurate for a skewed null and handles edge cases", {
  set.seed(14)
  T0 <- rgamma(2e5, shape = 3, scale = 2) + 1          # a shifted gamma is Pearson III
  for (q in c(0.9, 0.95, 0.99))
    expect_equal(.right_tail_p(qgamma(q, shape = 3, scale = 2) + 1, T0), 1 - q,
                 tolerance = 0.05)
  S <- c(-1, 1) %o% seq_len(50)                        # symmetric: normal fallback
  expect_equal(.right_tail_p(10, as.vector(S)),
               pnorm(10, 0, sd(as.vector(S)), lower.tail = FALSE))
  expect_true(is.na(.right_tail_p(2, c(1, 1, 1))))     # zero spread
})

test_that("Block.resid.test returns the empirical and the parametric p-value", {
  set.seed(15)
  d <- sim_1(n = 150, p = 2)
  r  <- Block.resid.test(Y ~ X | Z1 + Z2, d, B = 99, seed = 1)
  rp <- Block.resid.test(Y ~ X | Z1 + Z2, d, B = 99, seed = 1, pvalue = "parametric")
  expect_equal(r$pvalue, "empirical")
  expect_equal(r$p.value, r$p.empirical)
  expect_equal(r$p.empirical, (1 + sum(r$perm_stats >= r$statistic)) / 100)
  expect_equal(rp$p.value, rp$p.parametric)
  # the choice changes only which p-value is reported
  expect_equal(c(rp$p.empirical, rp$p.parametric), c(r$p.empirical, r$p.parametric))
  expect_equal(r$p.parametric, .right_tail_p(r$statistic[[1]], r$perm_stats))
  expect_output(print(rp), "parametric p-value")
  for (a in list(list(stat = "cc1"), list(perm = "simple"))) {
    x <- do.call(Block.resid.test, c(list(Y ~ X | Z1 + Z2, d, B = 49, seed = 1), a))
    expect_true(x$p.parametric >= 0 && x$p.parametric <= 1)
  }
  w <- Block.resid.test(Y ~ X | Z1 + Z2, d, stat = "wgcm", B = 49, seed = 1)
  expect_true(is.na(w$p.parametric))
  expect_equal(w$p.value, w$p.empirical)
  expect_error(Block.resid.test(Y ~ X | Z1, d, stat = "wgcm", pvalue = "parametric"),
               "not available")
})

test_that("Block.resid.test runs with KNN, with and without the linear step", {
  skip_if_not_installed("FNN")
  d <- sim_1(n = 200, p = 3)
  r <- Block.resid.test(Y ~ X | Z1 + Z2, d, method = "KNN", B = 19, nfolds = 5, seed = 1)
  expect_equal(r$regression, "KNN")
  expect_true(r$linear_first)
  expect_true(all(r$linear_ncomp >= 1))
  expect_named(r$complexity, c("X", "Y"))
  expect_true(all(r$complexity >= 3 & r$complexity <= 100))
  expect_match(r$method, "PLS + KNN", fixed = TRUE)
  r0 <- Block.resid.test(Y ~ X | Z1 + Z2, d, method = "KNN", linear_first = FALSE, B = 19,
                         nfolds = 5, method_args = list(k = 8), seed = 1)
  expect_false(r0$linear_first)
  expect_equal(unname(r0$complexity), c(8, 8))
  for (args in list(list(reduce = "pls"), list(reduce = "pca"), list(stat = "wgcm"),
                    list(resid = "fitted"))) {
    r <- do.call(Block.resid.test, c(list(Y ~ X | Z1 + Z2, d, method = "KNN", B = 19,
                                          nfolds = 5, seed = 1), args))
    expect_true(r$p.value > 0 && r$p.value <= 1)
  }
  r1 <- Block.resid.test(Y ~ X | Z1, d, method = "KNN", B = 99, seed = 1)   # Z2 missing: H0 false
  expect_lt(r1$p.value, 0.05)
})
