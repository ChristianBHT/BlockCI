# helpers: resampling of X, losses, p-values, method lookup

test_that("rotate_X keeps the X-V relation and the residual covariance", {
  set.seed(1); n <- 100
  V <- matrix(rnorm(n * 3), n); X <- V %*% matrix(rnorm(3 * 2), 3) + matrix(rnorm(n * 2), n)
  Xr <- rotate_X(X, V)
  expect_equal(dim(Xr), dim(X))
  Q <- qr(cbind(1, V))
  expect_equal(qr.fitted(Q, Xr), qr.fitted(Q, X))                                   # same V part
  expect_equal(crossprod(Xr - qr.fitted(Q, Xr)), crossprod(X - qr.fitted(Q, X)))  # same E'E
})

test_that("rotate_X works without V", {
  set.seed(2); X <- matrix(rnorm(60), 30)
  expect_equal(dim(rotate_X(X, NULL)), dim(X))
})

test_that("permute_X and permute_X_cond permute whole rows", {
  set.seed(3); n <- 100
  V <- matrix(rnorm(n * 2), n); X <- cbind(seq_len(n), seq_len(n) * 10)
  Xp <- permute_X(X)
  expect_setequal(Xp[, 1], X[, 1])
  expect_equal(Xp[, 2], Xp[, 1] * 10)                     # rows kept together
  Xc <- permute_X_cond(X, V, group_size = 10)
  expect_setequal(Xc[, 1], X[, 1])
  expect_equal(Xc[, 2], Xc[, 1] * 10)
  expect_equal(dim(permute_X_cond(X, NULL)), dim(X))
})

test_that("permute_X_cond only swaps rows within V groups", {
  set.seed(4); n <- 120
  V <- matrix(rep(c(-10, 0, 10), each = n / 3) + rnorm(n, sd = 0.01), n)
  X <- matrix(seq_len(n), n)
  Xc <- permute_X_cond(X, V, group_size = 40)
  expect_true(all(ceiling(Xc[, 1] / 40) == ceiling(seq_len(n) / 40)))   # stays in its block
})

test_that(".make_resampler returns X* of the right size for every scheme", {
  set.seed(5); n <- 60
  V <- matrix(rnorm(n * 2), n); X <- matrix(rnorm(n * 3), n)
  for (rs in c("rotation", "permutation", "conditional")) {
    f <- .make_resampler(X, V, rs, group_size = 10)
    expect_equal(dim(f()), dim(X))
  }
  expect_equal(dim(.make_resampler(X, NULL, "conditional")()), dim(X))
})

test_that("rmsep and .parametric_p", {
  expect_equal(rmsep(matrix(1:4, 2), matrix(1:4, 2)), 0)
  expect_equal(rmsep(c(0, 0), c(3, 4)), sqrt(12.5))
  T0 <- c(-1, 0, 1)
  expect_equal(.parametric_p(0, T0), 0.5)
  expect_lt(.parametric_p(-5, T0), 0.01)
  expect_true(is.na(.parametric_p(0, c(1, 1, 1))))      # zero spread
})

test_that(".get_learner finds methods case-insensitively and rejects unknown ones", {
  expect_identical(.get_learner("pls"), .cci_methods$PLS)
  expect_identical(.get_learner("rf", .resid_methods), .resid_methods$RF)
  f <- function(...) NULL
  expect_identical(.get_learner(f), f)
  expect_error(.get_learner("nope"), "unknown method")
  expect_error(.get_learner(c("PLS", "RF")), "must be one of")
  expect_equal(.method_name("xgb"), "XGB")
  expect_equal(.method_name(f), "custom")
})
