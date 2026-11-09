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
#' @param split Define the training rows (optional; if `NULL` a random split is drawn) 
#' @return named vector `M0, M1, T = M1 - M0, a, b` (RMSEP metric)
#' @export
sopls_test_gen <- function(Y, X, V = NULL, p = 0.75, ncomp = "auto", scale = TRUE,
                           max_ncomp = c(6, 10), cv_segments = 5,
                           split = NULL, losses = FALSE) {
  n <- nrow(Y)
  train <- if (is.null(split)) sample.int(n, floor(p * n)) else split
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
    l0 <- rowSums(Ys[test, , drop = FALSE]^2)
    if (nb == 0) { M1 <- M0; l1 <- l0 } else {
      pr <- predict(fit, newdata = list(Y = Ys[test, , drop = FALSE], X = X[test, , drop = FALSE]), ncomp = nb)
      M1 <- rmsep(Ys[test, ], pr[, , dim(pr)[3]])
      l1 <- rowSums((Ys[test, , drop = FALSE] - matrix(pr[, , dim(pr)[3]], length(test)))^2)
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
    yh0 <- pr[, , which(cn == paste(c(av, 0), collapse = ","))]
    yh1 <- pr[, , which(cn == paste(c(av, nb), collapse = ","))]
    M0 <- rmsep(Ys[test, ], yh0)
    M1 <- rmsep(Ys[test, ], yh1)
    l0 <- rowSums((Ys[test, , drop = FALSE] - matrix(yh0, length(test)))^2)
    l1 <- rowSums((Ys[test, , drop = FALSE] - matrix(yh1, length(test)))^2)
    na <- sum(av)
  }
  out <- c(M0 = M0, M1 = M1, T = M1 - M0, a = na, b = nb)
  if (!losses) return(out)
  # per-test-row squared errors (summed over Y columns) for the V-only and V + X models
  list(stats = out, test = test, l0 = l0, l1 = l1)
}
