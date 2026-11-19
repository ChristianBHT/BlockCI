# Dimension reduction of the conditioning blocks V --------------------------------------
# Used by cci_test_gen() / Block.CCI.test() / Block.CCI.fast() (reduce = "pca" or "pls")
# and by Block.resid.test() (reduce = "pca"; "pls" is done inside the XGB and RF
# residual learners, from their linear PLS step). The reductions are fitted on the
# training rows and applied to all rows.

# Standardise the columns of M with the means and sds of the rows `rows`; columns that
# are constant on those rows are dropped.
.std_cols <- function(M, rows) {
  mu <- colMeans(M[rows, , drop = FALSE])
  sdv <- apply(M[rows, , drop = FALSE], 2, stats::sd)
  keep <- is.finite(sdv) & sdv > 0
  sweep(sweep(M[, keep, drop = FALSE], 2, mu[keep]), 2, sdv[keep], "/")
}

# Principal component scores of V (n x k), fitted on the rows `train` (V standardised).
# k = "auto": the number of eigenvalues of the correlation matrix of V above the
# Marchenko-Pastur edge (1 + sqrt(p / n))^2, i.e. the components that stand out from
# noise (factor structure), at most `max_k`. If none does, the smallest number of
# components explaining `var_explained` of the variance (at most `max_k`).
.reduce_pca <- function(Vm, train, k = "auto", max_k = 20, var_explained = 0.9) {
  Vs <- .std_cols(Vm, train)
  auto <- identical(k, "auto")
  nt <- length(train); p <- ncol(Vs)
  kmax <- min(p, nt - 1, if (auto) max_k else k)
  sv <- svd(Vs[train, , drop = FALSE], nu = 0, nv = kmax)
  if (auto) {
    ev <- sv$d^2 / (nt - 1)                       # eigenvalues of the correlation matrix
    k <- sum(ev > (1 + sqrt(p / nt))^2)
    if (k == 0) k <- which(cumsum(ev) / sum(ev) >= var_explained)[1]
    if (is.na(k)) k <- kmax
  }
  k <- max(1L, min(as.integer(k), kmax))
  Vs %*% sv$v[, seq_len(k), drop = FALSE]
}

# PLS scores of V (n x k) from a PLS regression of the target Ts (e.g. the standardised
# Y) on V, fitted on the rows `train` (V standardised). k = "auto": number of components
# by cross-validation on the training rows, at most `max_k`.
.reduce_pls <- function(Vm, Ts, train, k = "auto", max_k = 10, segments = 5) {
  Vs <- .std_cols(Vm, train)
  auto <- identical(k, "auto")
  kmax <- min(ncol(Vs), length(train) - 2, if (auto) max_k else k)
  dtr <- data.frame(Tgt = I(as.matrix(Ts)[train, , drop = FALSE]), V = I(Vs[train, , drop = FALSE]))
  cv_args <- if (auto) list(validation = "CV", segments = min(segments, length(train))) else list()
  fit <- do.call(pls::plsr, c(list(Tgt ~ V, data = dtr, ncomp = kmax, scale = FALSE), cv_args))
  if (auto) {
    cv <- pls::MSEP(fit, estimate = "CV", intercept = FALSE)$val   # [1, response, component]
    k <- which.min(apply(cv[1, , , drop = FALSE], 3, sum))
  }
  k <- max(1L, min(as.integer(k), kmax))
  sc <- stats::predict(fit, newdata = data.frame(V = I(Vs)), comps = seq_len(k), type = "scores")
  matrix(unclass(sc), nrow(Vm), k)
}

# Replace the list of V blocks by one block of scores (CCI tests).
.reduce_blocks <- function(Vl, Ys, train, reduce, k = "auto") {
  if (is.null(Vl) || reduce == "none") return(Vl)
  Vm <- do.call(cbind, Vl)
  list(switch(reduce,
              pca = .reduce_pca(Vm, train, k),
              pls = .reduce_pls(Vm, Ys, train, k)))
}
