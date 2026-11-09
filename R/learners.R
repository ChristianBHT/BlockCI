# Learner wrappers for the CCI-type test -------------------------------------
# Every wrapper has the same interface and is chosen in sopls_test_gen() through
# its `learner` argument. Inputs: the standardised response `Ys` (all rows), the
# X block, the list of V blocks (or NULL), and the train / test row indices.
# Output: list(yh0, yh1, a, b) where
#   yh0 = test-row predictions of Y from V only (zeros if there is no V),
#   yh1 = test-row predictions of Y from V and X,
#   a, b = complexity chosen for the V-only and the V + X model (for SO-PLS:
#          number of V components and of added X components; for XGBoost: the
#          number of boosting rounds of the two models).
# All tuning uses TRAINING rows only; the test rows are never touched.

.learner_sopls <- function(Ys, X, Vl, train, test, ncomp, scale, max_ncomp, cv_segments) {
  q <- ncol(Ys); nt <- length(test)
  auto <- identical(ncomp, "auto")
  lim <- length(train) - 2
  if (is.null(Vl)) {   # unconditional test: V-only predictor is the training mean (0 on the Ys scale)
    nb <- min(if (auto) max_ncomp[2] else ncomp[2], ncol(X), lim)
    dtr <- list(Y = Ys[train, , drop = FALSE], X = X[train, , drop = FALSE])
    fit <- multiblock::sopls(Y ~ X, ncomp = nb, data = dtr, scale = scale, progress = FALSE,
                             validation = if (auto) "CV" else "none", segments = cv_segments)
    if (auto) {
      cvr <- fit$validation$RMSECV
      nb <- as.integer(sub(".*,", "", names(cvr)[which.min(cvr)]))
    }
    yh0 <- matrix(0, nt, q)
    if (nb == 0) yh1 <- yh0 else {
      pr <- predict(fit, newdata = list(Y = Ys[test, , drop = FALSE], X = X[test, , drop = FALSE]), ncomp = nb)
      yh1 <- matrix(pr[, , dim(pr)[3]], nt)
    }
    return(list(yh0 = yh0, yh1 = yh1, a = 0, b = nb))
  }
  # V is a list of blocks; each V block gets ONE component (ncomp = "auto") or the
  # number given in ncomp = c(a_1, ..., a_k, b).
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
           stats::setNames(lapply(Vl, function(m) m[train, , drop = FALSE]), vn),
           list(X = X[train, , drop = FALSE]))
  dte <- c(list(Y = Ys[test, , drop = FALSE]),
           stats::setNames(lapply(Vl, function(m) m[test, , drop = FALSE]), vn),
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
  list(yh0 = matrix(pr[, , which(cn == paste(c(av, 0), collapse = ","))], nt),
       yh1 = matrix(pr[, , which(cn == paste(c(av, nb), collapse = ","))], nt),
       a = sum(av), b = nb)
}

# XGBoost with multi-target regression (xgboost >= 2.0). V-only model: Y ~ V; full
# model: Y ~ V + X. Number of rounds chosen by early stopping on an inner
# validation part of the TRAINING rows.
# learner_args (all optional): nrounds (max rounds, default 300),
#   early_stopping_rounds (10), val_fraction (0.2), params (list merged into the
#   defaults: max_depth = 3, eta = 0.1, tree_method = "hist",
#   multi_strategy = "multi_output_tree", nthread = 1).
.learner_xgboost <- function(Ys, X, Vl, train, test, learner_args = list()) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("The package 'xgboost' (>= 2.0) is required for learner = \"xgboost\".")
  q <- ncol(Ys); nt <- length(test)
  la <- utils::modifyList(list(nrounds = 300, early_stopping_rounds = 10, val_fraction = 0.2, params = list()),
                          learner_args)
  params <- utils::modifyList(
    list(objective = "reg:squarederror", tree_method = "hist", multi_strategy = "multi_output_tree",
         max_depth = 3, eta = 0.1, nthread = 1),
    la$params)

  fit_pred <- function(F) {
    nv <- max(2L, ceiling(la$val_fraction * length(train)))
    va <- sample(train, nv)                       # inner validation rows (training part only)
    tr2 <- setdiff(train, va)
    dtr <- xgboost::xgb.DMatrix(F[tr2, , drop = FALSE], label = Ys[tr2, , drop = FALSE])
    dva <- xgboost::xgb.DMatrix(F[va, , drop = FALSE],  label = Ys[va, , drop = FALSE])
    fit <- xgboost::xgb.train(params = params, data = dtr, nrounds = la$nrounds,
                              evals = list(val = dva),
                              early_stopping_rounds = la$early_stopping_rounds, verbose = 0)
    best <- xgboost::xgb.attr(fit, "best_iteration")
    iter <- if (is.null(best)) xgboost::xgb.get.num.boosted.rounds(fit) else as.integer(best) + 1L
    list(pred = matrix(stats::predict(fit, F[test, , drop = FALSE]), nt, q), iter = iter)
  }

  if (is.null(Vl)) {
    f1 <- fit_pred(X)
    return(list(yh0 = matrix(0, nt, q), yh1 = f1$pred, a = 0, b = f1$iter))
  }
  Vm <- do.call(cbind, Vl)
  f0 <- fit_pred(Vm)
  f1 <- fit_pred(cbind(Vm, X))
  list(yh0 = f0$pred, yh1 = f1$pred, a = f0$iter, b = f1$iter)
}
