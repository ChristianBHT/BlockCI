# dimension reduction of the conditioning blocks (reduce = "pca" / "pls")

sim_factor <- function(n, p = 60, r = 2, effect = 0) {
  L <- matrix(rnorm(r * p), r, p)
  Fz <- matrix(rnorm(n * r), n, r)
  V <- Fz %*% L + matrix(rnorm(n * p), n, p)
  X <- Fz %*% matrix(rnorm(r * 2), r, 2) + matrix(rnorm(n * 2), n, 2)
  Y <- Fz %*% matrix(rnorm(r * 2), r, 2) + effect * X + matrix(rnorm(n * 2), n, 2)
  list(Y = Y, X = X, V = V)
}

test_that(".reduce_pca finds the factors and .reduce_pls returns n x k scores", {
  set.seed(11)
  d <- sim_factor(300, r = 2)
  tr <- 1:200
  S <- .reduce_pca(d$V, tr)
  expect_equal(dim(S), c(300L, 2L))                   # two factors above the MP edge
  expect_equal(ncol(.reduce_pca(d$V, tr, k = 5)), 5L)
  P <- .reduce_pls(d$V, scale(d$Y), tr, k = 3)
  expect_equal(dim(P), c(300L, 3L))
  expect_true(ncol(.reduce_pls(d$V, scale(d$Y), tr)) >= 1)
})

test_that("Block.CCI.test and Block.CCI.fast run with reduce", {
  set.seed(12)
  d0 <- sim_factor(200)
  d1 <- sim_factor(200, effect = 1)
  for (rd in c("pca", "pls")) {
    r0 <- Block.CCI.test(Y ~ X | V, d0, method = "PLS", nperm = 39, reduce = rd,
                         progress = FALSE, seed = 1)
    r1 <- Block.CCI.test(Y ~ X | V, d1, method = "PLS", nperm = 39, reduce = rd,
                         progress = FALSE, seed = 1)
    expect_equal(r1$reduce, rd)
    expect_gt(r0$p.value, 0.01)
    expect_lt(r1$p.value, 0.05)
  }
  f <- Block.CCI.fast(Y ~ X | V, d1, method = "PLS", nsplit = 3, reduce = "pca",
                      progress = FALSE, seed = 1)
  expect_true(f$p.value >= 0 && f$p.value <= 1)
})

test_that("Block.resid.test runs with reduce = 'pca' and with RF / XGB on PLS scores", {
  set.seed(13)
  d <- sim_factor(200)
  r <- Block.resid.test(Y ~ X | V, d, method = "PLS", reduce = "pca", B = 99, seed = 1)
  expect_equal(r$n_reduced, 2L)
  expect_gt(r$p.value, 0.01)
  rp <- Block.resid.test(Y ~ X | V, d, method = "PLS", reduce = "pls", B = 19, seed = 1)
  expect_equal(rp$reduce, "none")                      # ignored for PLS
  skip_if_not_installed("grf")
  rf <- Block.resid.test(Y ~ X | V, d, method = "RF", reduce = "pls", B = 19, nfolds = 3,
                         method_args = list(num.trees = 200), seed = 1)
  expect_equal(rf$reduce, "pls")
  expect_true(rf$p.value > 0 && rf$p.value <= 1)
  rf2 <- Block.resid.test(Y ~ X | V, d, method = "RF", reduce = "pls", reduce_ncomp = 3,
                          linear_first = FALSE, B = 19, nfolds = 3,
                          method_args = list(num.trees = 200), seed = 1)
  expect_true(rf2$p.value > 0 && rf2$p.value <= 1)
  skip_if_not_installed("xgboost")
  xg <- Block.resid.test(Y ~ X | V, d, method = "XGB", reduce = "pls", B = 19, nfolds = 3,
                         method_args = list(nrounds = 20), seed = 1)
  expect_true(xg$p.value > 0 && xg$p.value <= 1)
})
