test_that("rotate_X keeps X-V relation and residual covariance", {
  set.seed(1); n <- 100
  V <- matrix(rnorm(n * 3), n); X <- V %*% matrix(rnorm(3 * 2), 3) + matrix(rnorm(n * 2), n)
  Xr <- rotate_X(X, V)
  expect_equal(dim(Xr), dim(X))
  Q <- qr(cbind(1, V))
  expect_equal(qr.fitted(Q, Xr), qr.fitted(Q, X))                       # same V-explained part
  expect_equal(crossprod(Xr - qr.fitted(Q, Xr)), crossprod(X - qr.fitted(Q, X)))  # same E'E
})

test_that("rotate_X works without V", {
  set.seed(2); X <- matrix(rnorm(60), 30)
  expect_equal(dim(rotate_X(X, NULL)), dim(X))
})

test_that("Block.CCI.test (SOPLS) detects dependence and not independence", {
  set.seed(3); n <- 200
  V <- matrix(rnorm(n * 4), n); X <- V[, 1:2] + matrix(rnorm(n * 2), n)
  Y0 <- cbind(V[, 1] + rnorm(n), rnorm(n))
  Y1 <- Y0; Y1[, 1] <- Y1[, 1] + 1.5 * (X[, 1] - V[, 1])
  r0 <- Block.CCI.test(Y ~ X | V, list(Y = Y0, X = X, V = V), nperm = 39, ncomp = c(4, 2), progress = FALSE, seed = 1)
  r1 <- Block.CCI.test(Y ~ X | V, list(Y = Y1, X = X, V = V), nperm = 39, ncomp = c(4, 2), progress = FALSE, seed = 1)
  expect_s3_class(r1, "BlockCCI")
  expect_gt(r0$p.value, 0.05)
  expect_lt(r1$p.value, 0.05)
})

test_that("permutation option and unconditional test run", {
  set.seed(4); n <- 100
  X <- matrix(rnorm(n * 3), n); Y <- matrix(rnorm(n * 2), n)
  r <- Block.CCI.test(Y ~ X, list(Y = Y, X = X), nperm = 9, resampling = "permutation", progress = FALSE)
  expect_true(r$p.value > 0 && r$p.value <= 1)
})

test_that("Block.CCI.fast runs with SOPLS, PLS, ROSA and RF", {
  set.seed(7)
  d <- sim_1(n = 200, p = 3)
  for (m in c("SOPLS", "PLS", "ROSA")) {
    r <- Block.CCI.fast(Y ~ X | Z1 + Z2, d, method = m, nsplit = 3, seed = 1, progress = FALSE)
    expect_s3_class(r, "BlockCCIfast")
    expect_equal(r$method, m)
    expect_true(is.na(r$p.value) || (r$p.value >= 0 && r$p.value <= 1))
  }
  skip_if_not_installed("grf")
  r <- Block.CCI.fast(Y ~ X | Z1 + Z2, d, method = "RF", nsplit = 3,
                      method_args = list(num.trees = 100), seed = 1, progress = FALSE)
  expect_equal(r$method, "RF")
})

test_that("Block.CCI.test with ROSA detects a false null", {
  set.seed(8)
  d <- sim_1(n = 300, p = 3)
  r <- Block.CCI.test(Y ~ X | Z1, d, method = "ROSA", nperm = 39, progress = FALSE, seed = 1)
  expect_s3_class(r, "BlockCCI")
  expect_lt(r$p.value, 0.05)
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
