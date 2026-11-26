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

test_that("rotate_X draws rotations with no preferred sign", {
  set.seed(14); n <- 10
  V <- matrix(rnorm(n * 2), n); X <- matrix(rnorm(n), n)
  Q <- qr(cbind(1, V))
  E_star <- replicate(4000, qr.resid(Q, rotate_X(X, V))[, 1])      # n x draws
  # under a uniform rotation every rotated residual is symmetric around 0
  expect_true(all(abs(rowMeans(E_star > 0) - 0.5) < 0.04))         # about 5 s.e.
})

test_that("rotate_X works without V", {
  set.seed(2); X <- matrix(rnorm(60), 30)
  expect_equal(dim(rotate_X(X, NULL)), dim(X))
})

test_that("permute_X and permute_X_cond keep rows together", {
  set.seed(3); n <- 100
  V <- matrix(rnorm(n * 2), n); X <- cbind(seq_len(n), seq_len(n) * 10)
  Xp <- permute_X(X)
  expect_setequal(Xp[, 1], X[, 1])
  expect_equal(Xp[, 2], Xp[, 1] * 10)                     # rows kept together
  Xc <- permute_X_cond(X, V)
  expect_equal(dim(Xc), dim(X))
  expect_true(all(Xc[, 1] %in% X[, 1]))                   # only rows of X (may repeat)
  expect_equal(Xc[, 2], Xc[, 1] * 10)
  expect_equal(dim(permute_X_cond(X, NULL)), dim(X))
})

test_that(".knn_index finds the nearest neighbours in the maximum norm of the ranks", {
  set.seed(4)
  for (d in c(1, 2, 6)) {
    n <- 150; k <- 7
    V <- matrix(rnorm(n * d), n)                          # no ties, so ranks are fixed
    D <- as.matrix(stats::dist(apply(V, 2, rank), method = "maximum"))
    diag(D) <- -1                                         # the row itself comes first
    nb <- .knn_index(V, k)
    expect_equal(dim(nb), c(n, k))
    for (i in seq_len(n)) {
      expect_true(i %in% nb[i, ])
      expect_equal(anyDuplicated(nb[i, ]), 0L)
      expect_equal(sort(unname(D[i, nb[i, ]])), sort(unname(D[i, ]))[seq_len(k)])
    }
  }
})

test_that(".local_perm follows the greedy scheme", {
  set.seed(5); n <- 200
  nb <- .knn_index(matrix(rnorm(n * 2), n), 5)
  for (r in 1:20) {
    idx <- .local_perm(nb)
    expect_true(all(rowSums(nb == idx) == 1))             # always one of the row's neighbours
    # an index is used again only by a row whose neighbours were all taken. So among the
    # rows sharing an index, at most one (the first to take it) may have a free neighbour
    blocked <- apply(nb, 1, function(s) all(s %in% idx))
    expect_true(all(tapply(!blocked, idx, sum) <= 1))
    expect_gt(length(unique(idx)), 0.8 * n)
  }
  # neighbourhoods that do not overlap: the result is a strict permutation within them
  blocks <- t(sapply(rep(1:40, each = 5), function(g) (g - 1) * 5 + 1:5))
  idx <- .local_perm(blocks)
  expect_setequal(idx, 1:200)
  expect_true(all(ceiling(idx / 5) == ceiling(1:200 / 5)))
})

test_that("permute_X_cond keeps X* close in V, also with discrete V", {
  set.seed(6); n <- 300
  V <- matrix(sample(1:3, n, TRUE), n)                    # many ties
  X <- matrix(seq_len(n), n)
  Xc <- permute_X_cond(X, V, k_perm = 5)
  expect_gt(mean(V[Xc[, 1]] == V), 0.95)                  # nearly always the same level
  expect_gt(mean(Xc[, 1] != seq_len(n)), 0.5)             # but mostly another row
  V <- matrix(rnorm(n), n)
  Xc <- permute_X_cond(X, V, k_perm = 5)
  expect_lte(max(abs(rank(V)[Xc[, 1]] - rank(V))), 4)     # within 4 ranks
})

test_that("k_perm is checked", {
  X <- matrix(rnorm(20), 10); V <- matrix(rnorm(10), 10)
  for (bad in list(1, 11, 2.5, NA, c(3, 4), "5"))
    expect_error(permute_X_cond(X, V, k_perm = bad), "k_perm")
  expect_equal(dim(permute_X_cond(X, V, k_perm = 10)), dim(X))
})

test_that(".make_resampler returns X* of the right size for every scheme", {
  set.seed(5); n <- 60
  V <- matrix(rnorm(n * 2), n); X <- matrix(rnorm(n * 3), n)
  for (rs in c("rotation", "permutation", "conditional")) {
    f <- .make_resampler(X, V, rs, k_perm = 5)
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
