# Block Computational Conditional Independence (CCI) test --------------------
# Tests Y _||_ X | V for matrices Y, X, V (V optional), mirroring the CCI package:
#   Block.CCI.test() ~ CCI.test()   : user interface (formula + data + method)
#   cci_null_dist()  ~ perm.test()  : builds the null distribution
#   cci_test_gen()   ~ test.gen()   : one MCCV statistic
#   learners.R                      : the regression methods ("SOPLS", "PLS", "XGB", "RF", "KNN")
# Residual test: Block.resid.test() with the residual learners in resid_learners.R
#   ("PLS", "ROSA", "XGB", "RF", "KNN") and .resid_perm_test() below.
# Requires: multiblock (sopls), pls (plsr); xgboost, grf and FNN (KNN) are optional

# ---- helpers ---------------------------------------------------------------

# `a %||% b`: a if it is not NULL, otherwise b
`%||%` <- function(a, b) if (is.null(a)) b else a

# Show a message only once per R session (per distinct text)
.msg_env <- new.env(parent = emptyenv())
.message_once <- function(txt) {
  if (is.null(.msg_env[[txt]])) {
    assign(txt, TRUE, envir = .msg_env)
    message(txt)
  }
  invisible(NULL)
}

.get_block <- function(vars, data) {
  if (length(vars) == 0) return(NULL)
  if (is.data.frame(data)) return(as.matrix(data[, vars, drop = FALSE]))
  do.call(cbind, lapply(vars, function(v) as.matrix(data[[v]])))
}

# Formula "Y ~ X | V1 + V2" -> list(y, x, v) of variable names (v may be empty)
.parse_cci_formula <- function(formula) {
  rhs <- formula[[3]]
  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))) {
    xv <- all.vars(rhs[[2]]); vv <- all.vars(rhs[[3]])
  } else { xv <- all.vars(rhs); vv <- character(0) }
  list(y = all.vars(formula[[2]]), x = xv, v = vv)
}

# Parse `formula` (Y ~ X | V1 + V2) and `data` once. Returns the blocks Y, X and
# V (a list with one matrix per V block, or NULL) with incomplete rows removed, plus
# a canonical `data` list (Y, X, V1, V2, ...) and matching `formula`
# (Y ~ X | V1 + V2) in which X can be swapped for a resampled copy.
.canonical <- function(formula, data) {
  f <- .parse_cci_formula(formula)
  missing_vars <- setdiff(c(f$y, f$x, f$v), names(data))
  if (length(missing_vars)) stop("variables not found in data: ", paste(missing_vars, collapse = ", "))
  Y <- .get_block(f$y, data)
  X <- .get_block(f$x, data)
  V <- if (length(f$v) == 0) NULL else lapply(f$v, .get_block, data = data)
  ok <- stats::complete.cases(cbind(Y, X, if (!is.null(V)) do.call(cbind, V)))
  Y <- Y[ok, , drop = FALSE]
  X <- X[ok, , drop = FALSE]
  if (!is.null(V)) V <- lapply(V, function(m) m[ok, , drop = FALSE])
  vn <- if (is.null(V)) character(0) else paste0("V", seq_along(V))
  dat <- c(list(Y = Y, X = X), if (!is.null(V)) stats::setNames(V, vn))
  fml <- stats::as.formula(paste("Y ~ X", if (length(vn)) paste("|", paste(vn, collapse = " + "))))
  list(Y = Y, X = X, V = V, data = dat, formula = fml)
}

#' Rotation-test resample of X given V (Langsrud-style)
#'
#' Splits X = fitted(X|V) + E. E (rank r) is replaced by E* = A* S B' where A* is a
#' uniformly random orthonormal r-frame in the orthogonal complement of span(1, V).
#' Keeps the X-V relation and the residual covariance E'E exactly; destroys any
#' dependence of the X residual on Y. (Linear V-adjustment.)
#'
#' @param X numeric matrix to resample
#' @param V numeric matrix of conditioning variables, or `NULL`
#' @return matrix of the same dimension as `X`
#' @export
rotate_X <- function(X, V = NULL) {
  n <- nrow(X)
  Q <- qr(cbind(rep(1, n), V))
  fit <- qr.fitted(Q, X); E <- X - fit
  sv <- svd(E); r <- sum(sv$d > 1e-10 * max(sv$d[1], 1e-300))
  if (r == 0) return(X)
  G <- matrix(rnorm(n * r), n, r)
  G <- G - qr.fitted(Q, G)                      # into complement of span(1, V)
  qrG <- qr(G)
  # Householder QR fixes the signs of Q by a convention that depends on G, so Q alone
  # is not uniformly distributed. Flipping the columns so that diag(R) > 0 makes it so
  # (Mezzadri 2007).
  sgn <- ifelse(diag(qr.R(qrG)) < 0, -1, 1)
  Astar <- sweep(qr.Q(qrG), 2, sgn, "*")
  fit + Astar %*% (sv$d[1:r] * t(sv$v[, 1:r, drop = FALSE]))
}

#' Global row permutation of X
#' @param X numeric matrix
#' @return `X` with rows permuted
#' @export
permute_X <- function(X) X[sample.int(nrow(X)), , drop = FALSE]

# A function that returns a fresh resampled copy X* of X each time it is called.
# The neighbours in V for resampling = "conditional" are found once.
.make_resampler <- function(X, Vm, resampling, k_perm = 5) {
  switch(resampling,
         rotation    = function() rotate_X(X, Vm),
         permutation = function() permute_X(X),
         conditional = {
           if (is.null(Vm)) function() permute_X(X) else {
             nb <- .knn_index(Vm, k_perm)
             function() X[.local_perm(nb), , drop = FALSE]
           }
         })
}

#'  Root Mean Square Error of Prediction
#'  @param y numeric vector of observed values
#'  @param yhat numeric vector of predicted values
rmsep <- function(y, yhat) sqrt(mean((y - yhat)^2))

# Left-tail p-value assuming a normal null with the sample mean and sd of T0
.parametric_p <- function(Tobs, T0) {
  s <- stats::sd(T0)
  if (!is.finite(s) || s == 0) return(NA_real_)
  stats::pnorm(Tobs, mean = mean(T0), sd = s)
}



# ---- shared pieces for residual-based tests ----------------------------------

# dependence statistic between two matrices
.resid_stat_fun <- function(stat) {
  switch(stat,
         rv = function(A, Bm) {
           Sab <- crossprod(A, Bm)
           sum(Sab^2) / sqrt(sum(crossprod(A)^2) * sum(crossprod(Bm)^2))
         },
         cc1 = function(A, Bm) stats::cancor(A, Bm)$cor[1],
         dcor = {
           if (!requireNamespace("energy", quietly = TRUE))
             stop("The package 'energy' is required for stat = 'dcor'.")
           function(A, Bm) energy::dcor(A, Bm)
         })
}

# Permutation test of dependence between the residual matrices X_res and Y_res.
# 'conditional': greedy local permutation among the k_perm nearest neighbours in Z
# (Runge 2018; see knn_perm.R).
.resid_perm_test <- function(X_res, Y_res, Z, stat, perm, B, k_perm = 5, k0 = 7) {
  n <- nrow(X_res)
  X_res <- scale(X_res, scale = FALSE)
  Y_res <- scale(Y_res, scale = FALSE)
  # the weighted GCM has its own (Gaussian) null distribution; no permutations
  if (stat == "wgcm") return(.wgcm_test(X_res, Y_res, Z, B, k0))
  stat_fun <- .resid_stat_fun(stat)
  if (perm == "conditional") nb <- .knn_index(Z, k_perm)
  permuted <- function() {
    if (perm == "simple") return(Y_res[sample.int(n), , drop = FALSE])
    # the local scheme can repeat rows, so the columns are centred again
    M <- Y_res[.local_perm(nb), , drop = FALSE]
    sweep(M, 2, colMeans(M))
  }
  obs   <- stat_fun(X_res, Y_res)
  perms <- replicate(B, stat_fun(X_res, permuted()))
  list(obs = obs, perms = perms, p = (1 + sum(perms >= obs)) / (B + 1))
}

# ---- weighted generalised covariance measure (multivariate, fixed weights) --------
# Scheidegger, Hoerrmann and Buehlmann (2022), "The weighted generalised covariance
# measure", JMLR 23; Appendix A.1 (mWGCM.fix) with the default weights of A.1.1.

# Weight functions evaluated at the rows of Z: the constant 1 and, for every column
# d of Z, sign(z_d - a_dk) with a_dk the empirical k / (k0 + 1) quantile of Z_d,
# k = 1, ..., k0 (k0 * ncol(Z) + 1 weights; duplicates, e.g. from discrete Z, are
# removed).
.wgcm_weights <- function(Z, k0 = 7) {
  Z <- as.matrix(Z)
  W <- lapply(seq_len(ncol(Z)), function(d) {
    a <- unique(stats::quantile(Z[, d], probs = seq_len(k0) / (k0 + 1), names = FALSE))
    vapply(a, function(ak) sign(Z[, d] - ak), numeric(nrow(Z)))
  })
  W <- cbind(1, do.call(cbind, W))
  W[, !duplicated(t(W)), drop = FALSE]
}

# mWGCM.fix test. For every pair (X_j, Y_l) of residual columns and every weight
# w_k: R_i = eps_ij * xi_il * w_k(Z_i), T_jlk = sqrt(n) mean(R) / sd(R) (sd with
# divisor n), and S = max |T_jlk|. Under H0 the vector of T_jlk is approximately
# N(0, Sigma), Sigma the correlation matrix of the products. The p-value is
# P(max |N(0, Sigma_hat)| >= S), from B Gaussian multiplier draws
# T* = Rc' g / sqrt(n), g ~ N(0, I_n), Rc the centred and scaled products (this has
# covariance exactly Sigma_hat, without forming Sigma_hat). The products are
# rebuilt per pair and the draws are made in chunks, to limit memory.
.wgcm_test <- function(X_res, Y_res, Z, B, k0 = 7, chunk = 200) {
  n <- nrow(X_res)
  W <- .wgcm_weights(Z, k0)
  pairs <- expand.grid(j = seq_len(ncol(X_res)), l = seq_len(ncol(Y_res)))
  prods <- function(j, l) {
    R <- (X_res[, j] * Y_res[, l]) * W                 # n x (number of weights)
    m <- colMeans(R)
    s <- sqrt(pmax(colMeans(R^2) - m^2, 0))
    ok <- s > 0
    Tstat <- numeric(length(m)); Tstat[ok] <- sqrt(n) * m[ok] / s[ok]
    Rc <- sweep(R, 2, m)
    Rc <- sweep(Rc, 2, ifelse(ok, s, Inf), "/")         # columns with sd 0 become 0
    list(T = Tstat, Rc = Rc)
  }
  S_obs <- 0
  S_sim <- numeric(B)
  done <- 0L
  while (done < B) {
    b <- min(chunk, B - done)
    G <- matrix(stats::rnorm(n * b), n, b)               # the same draws for all pairs
    mx <- numeric(b)
    for (r in seq_len(nrow(pairs))) {
      pr <- prods(pairs$j[r], pairs$l[r])
      if (done == 0L) S_obs <- max(S_obs, abs(pr$T))
      Ts <- crossprod(pr$Rc, G) / sqrt(n)                # (number of weights) x b
      mx <- pmax(mx, apply(abs(Ts), 2, max))
    }
    S_sim[done + seq_len(b)] <- mx
    done <- done + b
  }
  list(obs = S_obs, perms = S_sim, p = (1 + sum(S_sim >= S_obs)) / (B + 1),
       n_stats = ncol(W) * nrow(pairs))
}
