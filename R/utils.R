# SO-PLS based Computational Conditional Independence test ------------------
# Tests Y _||_ X | V for matrices Y, X, V (V optional), mirroring the CCI package:
#   SOPLS.CCI.test()  ~ CCI.test()   : user interface (formula + data)
#   sopls_perm_test() ~ perm.test()  : builds the null distribution
#   sopls_test_gen()  ~ test.gen()   : one MCCV statistic
# Requires: multiblock (sopls)

# ---- helpers ---------------------------------------------------------------

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
  Astar <- qr.Q(qr(G))[, seq_len(r), drop = FALSE]
  fit + Astar %*% (sv$d[1:r] * t(sv$v[, 1:r, drop = FALSE]))
}

#' Global row permutation of X
#' @param X numeric matrix
#' @return `X` with rows permuted
#' @export
permute_X <- function(X) X[sample.int(nrow(X)), , drop = FALSE]

#' Conditional row permutation of X within groups of similar V
#'
#' Rows are clustered on the first (up to 3) principal components of V, with about
#' `group_size` rows per cluster, and the rows of X are permuted only within
#' clusters. This keeps the (possibly nonlinear) dependence between X and V
#' approximately intact while breaking the dependence between X and Y given V.
#'
#' @param X numeric matrix to resample
#' @param V numeric matrix of conditioning variables, or `NULL` (then a global
#'   permutation is used)
#' @param group_size desired number of rows per cluster
#' @return `X` with rows permuted within clusters
#' @export
permute_X_cond <- function(X, V = NULL, group_size = 10) {
  if (is.null(V)) return(permute_X(X))
  .permute_within(X, .v_groups(V, group_size))
}

# cluster labels of the rows of V (kmeans on up to 3 principal components)
.v_groups <- function(V, group_size = 10) {
  n <- nrow(V)
  pc <- stats::prcomp(V, scale. = TRUE)$x[, seq_len(min(ncol(V), 3)), drop = FALSE]
  stats::kmeans(pc, centers = max(2, n %/% group_size), nstart = 5)$cluster
}

.permute_within <- function(X, grp) {
  idx <- seq_len(nrow(X))
  for (g in unique(grp)) {
    i <- which(grp == g)
    idx[i] <- i[sample.int(length(i))]
  }
  X[idx, , drop = FALSE]
}

# A function that returns a fresh resampled copy X* of X each time it is called.
# The V clusters for method = "conditional" are computed once.
.make_resampler <- function(X, Vm, method, group_size = 10) {
  switch(method,
         rotation    = function() rotate_X(X, Vm),
         permutation = function() permute_X(X),
         conditional = {
           if (is.null(Vm)) function() permute_X(X) else {
             grp <- .v_groups(Vm, group_size)
             function() .permute_within(X, grp)
           }
         })
}

#'  Root Mean Square Error of Prediction
#'  @param y numeric vector of observed values
#'  @param yhat numeric vector of predicted values
rmsep <- function(y, yhat) sqrt(mean((y - yhat)^2))


