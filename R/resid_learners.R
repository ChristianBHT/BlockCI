# Residual learners for the block residual test ---------------------------------
# A residual learner removes the conditioning blocks Z from one target block M (X or
# Y) by regressing M on Z, and returns the residuals M - prediction.
#
# Interface:   function(M, Zl, folds, opts)
#   M       target matrix (n x q; the X or the Y block, on its original scale)
#   Zl      list of Z blocks (each n x p_k)
#   folds   integer vector of length n: the cross-fitting fold of each row. The same
#           folds are used for the X and the Y regression.
#   opts    list(ncomp, max_ncomp, scale, resid, linear_first, reduce, reduce_ncomp,
#           args): `resid` is "cv" (residuals
#           from predictions of held-out rows) or "fitted" (in-sample residuals);
#           `args` is the user's `method_args`. A learner uses the fields it needs.
# Returns list(res, complexity):
#   res         n x q residual matrix
#   complexity  number of components, boosting rounds, trees, ...
#   lin_ncomp   (optional) number of PLS components of a linear first step
#
# To add a method: write a function with this interface and add it to
# `.resid_methods` at the bottom of this file. A function can also be given
# directly as `method` in Block.resid.test().

# list of row indices per fold, for the `segments` argument of pls / multiblock
.fold_segments <- function(folds) unname(split(seq_along(folds), folds))

# PLS regression of M on cbind(Z) (pls::plsr). With ncomp = "auto" the number of
# components is the one with the lowest CV error over the folds; resid = "cv" uses
# the CV predictions of that model.
.resid_pls <- function(M, Zl, folds, opts) {
  Z <- do.call(cbind, Zl)
  n <- nrow(M)
  auto <- identical(opts$ncomp, "auto")
  k_max <- min(if (auto) opts$max_ncomp else opts$ncomp, ncol(Z), n - 2)
  need_cv <- auto || opts$resid == "cv"
  dat <- data.frame(M = I(M), Z = I(Z))
  cv_args <- if (need_cv) list(validation = "CV", segments = .fold_segments(folds)) else list()
  fit <- do.call(pls::plsr, c(list(M ~ Z, data = dat, ncomp = k_max, scale = opts$scale), cv_args))
  nc <- k_max
  if (auto) {
    cv <- pls::MSEP(fit, estimate = "CV", intercept = FALSE)$val   # [1, response, component]
    nc <- as.integer(which.min(apply(cv[1, , , drop = FALSE], 3, sum)))
  }
  pred <- if (opts$resid == "cv") fit$validation$pred[, , nc] else stats::predict(fit, newdata = dat, ncomp = nc)
  list(res = M - matrix(pred, nrow = n), complexity = as.integer(nc))
}

# ROSA (multiblock::rosa) of M on the Z blocks, kept as separate blocks. Components
# chosen and CV predictions made as for PLS.
.resid_rosa <- function(M, Zl, folds, opts) {
  n <- nrow(M)
  vn <- paste0("V", seq_along(Zl))
  fml <- stats::as.formula(paste("M ~", paste(vn, collapse = " + ")))
  dat <- c(list(M = M), stats::setNames(Zl, vn))
  auto <- identical(opts$ncomp, "auto")
  k_max <- min(if (auto) opts$max_ncomp else opts$ncomp, sum(vapply(Zl, ncol, 1L)), n - 2)
  need_cv <- auto || opts$resid == "cv"
  cv_args <- if (need_cv) list(validation = "CV", segments = .fold_segments(folds)) else list()
  fit <- do.call(multiblock::rosa, c(list(fml, data = dat, ncomp = k_max, scale = opts$scale), cv_args))
  nc <- k_max
  if (need_cv) {
    cvp <- fit$validation$pred                                    # n x q x ncomp
    nc <- min(nc, dim(cvp)[3])
    if (auto)
      nc <- which.min(vapply(seq_len(dim(cvp)[3]),
                             function(a) sum(colMeans((M - matrix(cvp[, , a], n))^2)), 0))
  }
  pred <- if (opts$resid == "cv") cvp[, , nc] else stats::fitted(fit)[, , nc]
  list(res = M - matrix(pred, nrow = n), complexity = as.integer(nc))
}

# Standardise M with the means and sds of the training rows `tr`, call `fit_pred(Ms)`
# (which fits on the training rows and returns predictions for the rows `te`), and
# return those predictions on the original scale of M.
.std_predict <- function(M, tr, te, fit_pred) {
  mu <- colMeans(M[tr, , drop = FALSE])
  sdv <- apply(M[tr, , drop = FALSE], 2, stats::sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
  Ms <- sweep(sweep(M, 2, mu), 2, sdv, "/")
  pr <- fit_pred(Ms)
  sweep(sweep(pr, 2, sdv, "*"), 2, mu, "+")
}

# ---- PLS first, then a nonlinear learner ----------------------------------------
# Linear step for XGB and RF: PLS regression of M on Z fitted on the training rows
# `tr`. With ncomp = "auto" the number of components is chosen by 5-fold CV inside
# the training rows (at most max_ncomp). Returns the fitted values on `tr`, the
# predictions for `te` and the number of components. With scores = TRUE also the
# PLS scores of the rows `tr` (F_tr) and `te` (F_te): the first nc scores, or the
# first opts$reduce_ncomp if that is a number (reduce = "pls").
.pls_fit_predict <- function(M, Z, tr, te, opts, scores = FALSE) {
  auto <- identical(opts$ncomp, "auto")
  lim <- min(ncol(Z), length(tr) - 2)
  k_max <- min(if (auto) opts$max_ncomp else opts$ncomp, lim)
  k_sc <- if (scores && is.numeric(opts$reduce_ncomp)) min(as.integer(opts$reduce_ncomp), lim) else NA
  k_fit <- max(k_max, k_sc, na.rm = TRUE)
  dtr <- data.frame(M = I(M[tr, , drop = FALSE]), Z = I(Z[tr, , drop = FALSE]))
  dte <- data.frame(Z = I(Z[te, , drop = FALSE]))
  cv_args <- if (auto) list(validation = "CV", segments = min(5, length(tr))) else list()
  fit <- do.call(pls::plsr, c(list(M ~ Z, data = dtr, ncomp = k_fit, scale = opts$scale), cv_args))
  nc <- k_max
  if (auto) {
    cv <- pls::MSEP(fit, estimate = "CV", intercept = FALSE)$val   # [1, response, component]
    nc <- as.integer(which.min(apply(cv[1, , seq_len(k_max), drop = FALSE], 3, sum)))
  }
  out <- list(fit_tr  = matrix(stats::fitted(fit)[, , nc], length(tr)),
              pred_te = matrix(stats::predict(fit, newdata = dte, ncomp = nc)[, , 1], length(te)),
              nc = nc)
  if (scores) {
    ks <- if (is.na(k_sc)) nc else k_sc
    out$F_tr <- matrix(unclass(fit$scores)[, seq_len(ks)], length(tr), ks)
    out$F_te <- matrix(unclass(stats::predict(fit, newdata = dte, comps = seq_len(ks),
                                              type = "scores")), length(te), ks)
  }
  out
}

# Cross-fitted predictions of M from Z with a nonlinear learner, optionally after a
# linear PLS step (opts$linear_first). `nl_fit(Tm, tr, te, Fm)` trains the nonlinear
# learner on the rows `tr` of the target Tm and returns predictions for the rows
# `te` (original scale of Tm); its inputs are the columns of Fm, or of Z when Fm is
# NULL. With the linear step, in each fold:
#   1. PLS of M on Z on the training rows; 2. the nonlinear learner on the training
#   residuals of that fit; 3. prediction of the held-out rows = PLS prediction +
#   nonlinear prediction.
# With opts$reduce = "pls" the nonlinear learner gets the PLS scores of that fold's
# PLS fit (a few supervised directions of Z) instead of Z itself; the PLS fit is
# then made also when linear_first = FALSE (its prediction is not used).
# resid = "cv": every fold is predicted by models trained on the other folds.
# resid = "fitted": one fit on all rows, in-sample predictions.
# Returns the predictions (n x q) and the mean number of PLS components.
.two_stage <- function(M, Z, folds, opts, nl_fit) {
  n <- nrow(M); q <- ncol(M)
  lin <- isTRUE(opts$linear_first)
  use_sc <- identical(opts$reduce, "pls")
  fit_rows <- function(tr, te) {
    if (!lin && !use_sc) return(list(pred = nl_fit(M, tr, te, NULL), nc = NA_real_))
    l <- .pls_fit_predict(M, Z, tr, te, opts, scores = use_sc)
    Fm <- NULL
    if (use_sc) {                                      # features for the tree learner
      Fm <- matrix(NA_real_, n, ncol(l$F_tr))
      Fm[te, ] <- l$F_te
      Fm[tr, ] <- l$F_tr
    }
    if (!lin) return(list(pred = nl_fit(M, tr, te, Fm), nc = l$nc))
    R <- M
    R[tr, ] <- M[tr, , drop = FALSE] - l$fit_tr       # training residuals of the linear fit
    list(pred = l$pred_te + nl_fit(R, tr, te, Fm), nc = l$nc)
  }
  pred <- matrix(NA_real_, n, q)
  if (opts$resid == "fitted") {
    r <- fit_rows(seq_len(n), seq_len(n))
    pred[] <- r$pred
    ncs <- r$nc
  } else {
    ncs <- numeric(0)
    for (k in unique(folds)) {
      te <- which(folds == k)
      r <- fit_rows(which(folds != k), te)
      pred[te, ] <- r$pred
      ncs <- c(ncs, r$nc)
    }
  }
  list(pred = pred, lin_ncomp = mean(ncs))
}

# Multi-target XGBoost (xgboost >= 2.0, multi_strategy = "multi_output_tree") of M
# on cbind(Z), by default after a linear PLS step (see .two_stage()). Each XGBoost
# model is trained with a fixed number of rounds (no early stopping) on the
# training rows only; targets are standardised with the training rows.
# opts$args (optional): nrounds (default 1000) and params (list merged into
#   max_depth = 3, eta = 0.1, tree_method = "hist",
#   multi_strategy = "multi_output_tree", nthread = 1).
.resid_xgb <- function(M, Zl, folds, opts) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("The package 'xgboost' (>= 2.0) is required for method = \"XGB\".")
  Z <- do.call(cbind, Zl)
  q <- ncol(M)
  la <- utils::modifyList(list(nrounds = 1000, params = list()), opts$args %||% list())
  prm <- utils::modifyList(
    list(objective = "reg:squarederror", tree_method = "hist", multi_strategy = "multi_output_tree",
         max_depth = 3, eta = 0.1, nthread = 1),
    la$params)
  nl_fit <- function(Tm, tr, te, Fm) .std_predict(Tm, tr, te, function(Ts) {
    Xf <- Fm %||% Z
    dtr <- xgboost::xgb.DMatrix(Xf[tr, , drop = FALSE], label = Ts[tr, , drop = FALSE])
    fit <- xgboost::xgb.train(params = prm, data = dtr, nrounds = la$nrounds, verbose = 0)
    matrix(stats::predict(fit, Xf[te, , drop = FALSE]), length(te), q)
  })
  r <- .two_stage(M, Z, folds, opts, nl_fit)
  list(res = M - r$pred, complexity = la$nrounds, lin_ncomp = r$lin_ncomp)
}

# Multi-output random forest (grf::multi_regression_forest) of M on cbind(Z), with
# splits that reduce the error summed over the columns of M; by default after a
# linear PLS step (see .two_stage()), cross-fitted over the folds. Targets are
# standardised with the training rows.
# Without the linear step (linear_first = FALSE) and with resid = "cv", the default
# is out-of-bag predictions of one forest grown on all rows (each row is predicted
# only by trees that did not use it); `oob = FALSE` cross-fits over the folds.
# opts$args (optional): oob (TRUE) and any arguments of
#   grf::multi_regression_forest(); defaults here num.trees = 2000 and
#   num.threads = 1 (NULL uses all cores).
.resid_rf <- function(M, Zl, folds, opts) {
  if (!requireNamespace("grf", quietly = TRUE))
    stop("The package 'grf' is required for method = \"RF\".")
  Z <- as.matrix(do.call(cbind, Zl))
  n <- nrow(M); q <- ncol(M)
  args <- opts$args %||% list()
  oob <- args$oob %||% TRUE
  args$oob <- NULL
  la <- utils::modifyList(list(num.trees = 2000, num.threads = 1), args)
  grow <- function(rows, Ts, oob_pred = FALSE, Xf = Z)
    do.call(grf::multi_regression_forest,
            c(list(X = Xf[rows, , drop = FALSE], Y = Ts[rows, , drop = FALSE],
                   compute.oob.predictions = oob_pred), la))
  if (opts$resid == "cv" && oob && !isTRUE(opts$linear_first) && !identical(opts$reduce, "pls")) {
    all_rows <- seq_len(n)
    pred <- .std_predict(M, all_rows, all_rows, function(Ts)
      matrix(stats::predict(grow(all_rows, Ts, oob_pred = TRUE))$predictions, n, q))
    return(list(res = M - pred, complexity = la$num.trees, lin_ncomp = NA_real_))
  }
  nl_fit <- function(Tm, tr, te, Fm) .std_predict(Tm, tr, te, function(Ts) {
    Xf <- Fm %||% Z
    matrix(stats::predict(grow(tr, Ts, Xf = Xf), newdata = Xf[te, , drop = FALSE])$predictions,
           length(te), q)
  })
  r <- .two_stage(M, Z, folds, opts, nl_fit)
  list(res = M - r$pred, complexity = la$num.trees, lin_ncomp = r$lin_ncomp)
}

# ---- registry -----------------------------------------------------------------
# Methods available by name in Block.resid.test(method = ...). Add new ones here.
.resid_methods <- list(
  PLS  = .resid_pls,
  ROSA = .resid_rosa,
  XGB  = .resid_xgb,
  RF   = .resid_rf
)
