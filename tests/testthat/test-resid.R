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
