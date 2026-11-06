# SO-PLS based Computational Conditional Independence test ------------------
# Tests Y _||_ X | V for matrices Y, X, V (V optional), mirroring the CCI package:
#   SOPLS.CCI.test()  ~ CCI.test()   : user interface (formula + data)
#   sopls_perm_test() ~ perm.test()  : builds the null distribution
#   sopls_test_gen()  ~ test.gen()   : one MCCV statistic
# Requires: multiblock (sopls)

# ---- helpers ---------------------------------------------------------------

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

rmsep <- function(y, yhat) sqrt(mean((y - yhat)^2))

# ---- one MCCV statistic (test.gen analogue) --------------------------------

#' One Monte Carlo CV statistic for the SO-PLS CCI test
#'
#' M0: SO-PLS of Y on V (a comps). M1: after adding b comps of X (sequentially, on
#' the V-orthogonalised Y and X). Both come from one `multiblock::sopls` fit.
#' With `ncomp = "auto"`, (a, b) are chosen by inner CV on the TRAINING part only
#' (test rows are never used for selection; identical for observed and null
#' statistics): a = argmin over `1:max_ncomp[1]` of RMSECV(a,0); then b = argmin
#' over `0:max_ncomp[2]` of RMSECV(a,b). b = 0 means X is not selected, so M1 = M0
#' and T = 0.
#'
#' @param Y,X,V numeric matrices (`V` may be `NULL`)
#' @param p training fraction
#' @param ncomp `c(a, b)` or `"auto"`
#' @param scale scale the blocks
#' @param max_ncomp maximum components searched when `ncomp = "auto"`
#' @param cv_segments CV segments for the inner selection
#' @return named vector `M0, M1, T = M1 - M0, a, b` (RMSEP metric)
#' @export
sopls_test_gen <- function(Y, X, V = NULL, p = 0.5, ncomp = "auto", scale = TRUE,
                           max_ncomp = c(6, 10), cv_segments = 3) {
  n <- nrow(Y) 
  train <- sample.int(n, floor(p * n)) 
  test <- setdiff(seq_len(n), train)
  mu <- colMeans(Y[train, , drop = FALSE]); 
  sdv <- if (scale) apply(Y[train, , drop = FALSE], 2, sd) else rep(1, ncol(Y))
  sdv[sdv == 0] <- 1
  Ys <- sweep(sweep(Y, 2, mu), 2, sdv, "/")
  auto <- identical(ncomp, "auto"); 
  lim <- length(train) - 2
  if (is.null(V)) {   # unconditional test: M0 = training-mean predictor
    nb <- min(if (auto) max_ncomp[2] else ncomp[2], ncol(X), lim)
    dtr <- list(Y = Ys[train, , drop = FALSE], X = X[train, , drop = FALSE])
    fit <- multiblock::sopls(Y ~ X, ncomp = nb, data = dtr, scale = scale, progress = FALSE,
                             validation = if (auto) "CV" else "none", segments = cv_segments)
    if (auto) { 
      cvr <- fit$validation$RMSECV; 
      nb <- as.integer(sub(".*,", "", names(cvr)[which.min(cvr)])) 
      }
    
    M0 <- rmsep(Ys[test, ], 0)
    if (nb == 0) M1 <- M0 else {
      pr <- predict(fit, newdata = list(Y = Ys[test, , drop = FALSE], X = X[test, , drop = FALSE]), ncomp = nb)
      M1 <- rmsep(Ys[test, ], pr[, , dim(pr)[3]])
    }
    na <- 0
  } else {
    # V may be one matrix or a list of blocks; each V block gets ONE component
    # (ncomp = "auto") or the number given in ncomp = c(a_1, ..., a_k, b).
    Vl <- if (is.list(V)) V else list(V)
    k <- length(Vl)

    vn <- paste0("V", seq_len(k))
    if (auto) {
      av <- rep(1L, k)
      nb <- max_ncomp[2]
    } else {
      if (length(ncomp) != k + 1) stop("ncomp must have length (number of V blocks) + 1")
      av <- ncomp[seq_len(k)]
      nb <- ncomp[k + 1]
    }
    av <- pmin(av, vapply(Vl, ncol, 1L), lim)
    nb <- min(nb, ncol(X), lim)
    dtr <- c(list(Y = Ys[train, , drop = FALSE]),
             setNames(lapply(Vl, function(m) m[train, , drop = FALSE]), vn),
             list(X = X[train, , drop = FALSE]))
    dte <- c(list(Y = Ys[test, , drop = FALSE]),
             setNames(lapply(Vl, function(m) m[test, , drop = FALSE]), vn),
             list(X = X[test, , drop = FALSE]))
    fml <- stats::as.formula(paste("Y ~", paste(c(vn, "X"), collapse = " + ")))
    fit <- multiblock::sopls(fml, ncomp = c(av, nb), data = dtr, scale = scale, progress = FALSE,
                             validation = if (auto) "CV" else "none", segments = cv_segments)
    if (auto) {   # V components fixed; choose b by inner CV
      cvr <- fit$validation$RMSECV
      ab <- do.call(rbind, lapply(strsplit(names(cvr), ","), as.integer))
      cand <- which(rowSums(abs(ab[, seq_len(k), drop = FALSE] -
                                matrix(av, nrow(ab), k, byrow = TRUE))) == 0)
      nb <- ab[cand[which.min(cvr[cand])], k + 1]
    }
    pr <- predict(fit, newdata = dte, ncomp = c(av, nb))
    cn <- dimnames(pr)[[3]]
    M0 <- rmsep(Ys[test, ], pr[, , which(cn == paste(c(av, 0), collapse = ","))])
    M1 <- rmsep(Ys[test, ], pr[, , which(cn == paste(c(av, nb), collapse = ","))])
    na <- sum(av)
  }
  c(M0 = M0, M1 = M1, T = M1 - M0, a = na, b = nb)
}

# ---- null distribution (perm.test analogue) --------------------------------

#' Null distribution of T for the SO-PLS CCI test
#'
#' R MCCV instances, each on a resampled X*, using `method = "rotation"`
#' (default) or `"permutation"`.
#'
#' @inheritParams sopls_test_gen
#' @param nperm number of null samples R
#' @param method `"rotation"` or `"permutation"`
#' @param progress show a progress bar
#' @return numeric vector of `nperm` null statistics
#' @export
sopls_perm_test <- function(Y, X, V = NULL, nperm = 160, method = c("rotation", "permutation"),
                            p = 0.5, ncomp = "auto", scale = TRUE, max_ncomp = c(10, 5),
                            cv_segments = 5, progress = FALSE) {
  method <- match.arg(method)
  if (progress) pb <- txtProgressBar(0, nperm, style = 3)
  T0 <- numeric(nperm)
  Vm <- if (is.list(V)) do.call(cbind, V) else V   # rotation adjusts for all of V
  for (r in seq_len(nperm)) {
    Xs <- if (method == "rotation") rotate_X(X, Vm) else permute_X(X)
    T0[r] <- sopls_test_gen(Y, Xs, V, p, ncomp, scale, max_ncomp, cv_segments)[["T"]]
    if (progress) setTxtProgressBar(pb, r)
  }
  if (progress) close(pb)
  T0
}

# ---- user interface (CCI.test analogue) ------------------------------------

#' SO-PLS Computational Conditional Independence test
#'
#' Tests Y independent of X given V using sequential orthogonal PLS: the test
#' statistic is T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V) on held-out data; the null
#' distribution comes from rotation (or permutation) of X. The p-value is the
#' share of null statistics <= the observed one.
#'
#' @param formula `Y ~ X | V` (V optional; names refer to matrices/columns in
#'   `data`; several blocks allowed, e.g. `Y ~ X | V1 + V2`, they are cbind-ed)
#' @param data named list of matrices (or a data.frame)
#' @param p training fraction of each MCCV split
#' @param nperm number of null samples R
#' @param ncomp `"auto"` (one component per V block, number of X components chosen
#'   by inner CV) or `c(a_1, ..., a_k, b)`: components for each of the k V blocks,
#'   then components added from X
#' @param method `"rotation"` (default) or `"permutation"` resampling of X
#' @param scale scale the blocks
#' @param max_ncomp,cv_segments settings for `ncomp = "auto"`
#' @param seed optional random seed
#' @param progress show a progress bar
#' @return object of class `c("SOPLSCCI", "htest")`
#' @export
SOPLS.CCI.test <- function(formula, data, p = 0.5, nperm = 160, ncomp = "auto",
                           method = c("rotation", "permutation"), scale = TRUE,
                           max_ncomp = c(10, 5), cv_segments = 5,
                           seed = NA, progress = TRUE) {
  method <- match.arg(method)
  if (!is.na(seed)) set.seed(seed)
  f <- .parse_cci_formula(formula)
  Y <- .get_block(f$y, data); X <- .get_block(f$x, data)
  # each name after "|" is its own V block (one latent variable per block)
  V <- if (length(f$v) == 0) NULL else lapply(f$v, .get_block, data = data)
  ok <- complete.cases(cbind(Y, X, if (!is.null(V)) do.call(cbind, V)))
  Y <- Y[ok, , drop = FALSE]
  X <- X[ok, , drop = FALSE]
  if (!is.null(V)) V <- lapply(V, function(m) m[ok, , drop = FALSE])
  obs <- sopls_test_gen(Y, X, V, p, ncomp, scale, max_ncomp, cv_segments)
  T0 <- sopls_perm_test(Y, X, V, nperm, method, p, ncomp, scale, max_ncomp, cv_segments, progress)
  pval <- (sum(T0 <= obs[["T"]]) + 1) / (nperm + 1)   # RMSEP is a loss: small T = evidence against H0
  structure(list(statistic = c(T = obs[["T"]]), p.value = pval, M0 = obs[["M0"]], M1 = obs[["M1"]],
                 null.distribution = T0, method = method, nperm = nperm, ncomp = ncomp,
                 ncomp_selected = c(V = obs[["a"]], X = obs[["b"]]), p = p,
                 formula = formula, tail = "left"),
            class = c("SOPLSCCI", "htest"))
}

#' @param x a `SOPLSCCI` object
#' @param ... unused
#' @rdname SOPLS.CCI.test
#' @export
print.SOPLSCCI <- function(x, ...) {
  cat("SO-PLS Computational Conditional Independence test\n")
  cat("Formula:", deparse(x$formula), 
      "| resampling:", x$method, 
      "| R =", x$nperm,
      "| ncomp =", paste(x$ncomp, collapse = ","),
      "(observed split uses V,X =", paste(x$ncomp_selected, collapse = ","), ")\n")
  cat(sprintf("RMSEP given V: %.4f   after adding X: %.4f   T = M1 - M0 = %.4f\n", x$M0, x$M1, x$statistic))
  cat(sprintf("p-value = %.4f (share of T* <= T_obs)\n", x$p.value))
  invisible(x)
}

