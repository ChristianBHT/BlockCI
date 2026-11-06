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

test_that("SOPLS.CCI.test detects dependence and not independence", {
  set.seed(3); n <- 200
  V <- matrix(rnorm(n * 4), n); X <- V[, 1:2] + matrix(rnorm(n * 2), n)
  Y0 <- cbind(V[, 1] + rnorm(n), rnorm(n))
  Y1 <- Y0; Y1[, 1] <- Y1[, 1] + 1.5 * (X[, 1] - V[, 1])
  r0 <- SOPLS.CCI.test(Y ~ X | V, list(Y = Y0, X = X, V = V), nperm = 39, ncomp = c(4, 2), progress = FALSE, seed = 1)
  r1 <- SOPLS.CCI.test(Y ~ X | V, list(Y = Y1, X = X, V = V), nperm = 39, ncomp = c(4, 2), progress = FALSE, seed = 1)
  expect_s3_class(r1, "SOPLSCCI")
  expect_gt(r0$p.value, 0.05)
  expect_lt(r1$p.value, 0.05)
})

test_that("permutation option and unconditional test run", {
  set.seed(4); n <- 100
  X <- matrix(rnorm(n * 3), n); Y <- matrix(rnorm(n * 2), n)
  r <- SOPLS.CCI.test(Y ~ X, list(Y = Y, X = X), nperm = 9, method = "permutation", progress = FALSE)
  expect_true(r$p.value > 0 && r$p.value <= 1)
})
