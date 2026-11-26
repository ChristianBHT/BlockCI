# Learners (regression methods) for the CCI-type tests ---------------------------
# A learner fits the two models of the test on the TRAINING rows and predicts the
# TEST rows:   M0: Y ~ V   and   M1: Y ~ V + X.
#
# Interface:   function(Ys, X, Vl, train, test, opts)
#   Ys           standardised response matrix (all rows)
#   X            the X block (all rows; the observed X or a resampled X*)
#   Vl           list of V blocks (all rows), or NULL for an unconditional test
#   train, test  row indices
#   opts         list(ncomp, scale, max_ncomp, cv_segments, args, fit_m0); `args`
#                is the user's `method_args` (learner-specific options). If
#                `fit_m0` is FALSE the caller already has the M0 predictions for
#                this split, and the learner may skip M0 (return yh0 = NULL,
#                a = NA). A learner uses the fields it needs and ignores the rest.
# Returns list(yh0, yh1, a, b):
#   yh0, yh1     test-row predictions of Y from M0 (zeros if there is no V) and M1
#   a, b         complexity of M0 and M1 (components, boosting rounds, trees, ...)
# All tuning must use the training rows only; the test rows are only predicted.
#
# To add a method: write a function with this interface and add it to
# `.cci_methods` at the bottom of this file under a short upper-case name. A
# learner function can also be given directly as `method` in Block.CCI.test(),
# Block.CCI.fast(), cci_test_gen() and cci_null_dist().

# SO-PLS (multiblock::sopls). One fit gives both models: M0 = the V components,
# M1 = after adding b components of X (sequential, on the V-orthogonalised Y and X).
# Each V block is its own SO-PLS block.
# opts: ncomp ("auto": one component per V block and b chosen by inner CV, or
#   c(a_1, ..., a_k, b)), scale, max_ncomp (maximum number of X components when
#   ncomp = "auto"; default 5), cv_segments.
.learner_sopls <- function(Ys, X, Vl, train, test, opts) {
  ncomp <- opts$ncomp %||% "auto"
  scale <- opts$scale %||% TRUE
  max_x <- opts$max_ncomp %||% 5
  cv_segments <- opts$cv_segments %||% 5
  q <- ncol(Ys); nt <- length(test)
  auto <- identical(ncomp, "auto")
  lim <- length(train) - 2
  if (is.null(Vl)) {   # unconditional test: M0 is the training mean (0 on the Ys scale)
    nb <- min(if (auto) max_x else ncomp[length(ncomp)], ncol(X), lim)
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
  k <- length(Vl)
  vn <- paste0("V", seq_len(k))
  if (auto) {
    av <- rep(1L, k)
    nb <- max_x
  } else {
    if (length(ncomp) != k + 1) stop("for SOPLS, ncomp must have length (number of V blocks) + 1")
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

# Regular PLS regression (pls::plsr). M0: Y ~ V; M1: Y ~ V + X, where all V blocks
# and X are cbind-ed (no block structure, unlike SO-PLS).
# opts: ncomp ("auto": chosen by inner CV for each model, or c(a, b) for M0 and M1),
#   scale, max_ncomp (maximum components of each model when ncomp = "auto";
#   default 10), cv_segments.
.learner_pls <- function(Ys, X, Vl, train, test, opts) {
  ncomp <- opts$ncomp %||% "auto"
  scale <- opts$scale %||% TRUE
  max_ncomp <- opts$max_ncomp %||% 10
  cv_segments <- opts$cv_segments %||% 5
  q <- ncol(Ys); nt <- length(test)
  auto <- identical(ncomp, "auto")
  if (!auto && length(ncomp) != 2) stop("for PLS, ncomp must be \"auto\" or c(a, b): components for Y ~ V and for Y ~ V + X")

  fit_pred <- function(Fm, k_fixed) {
    kmax <- min(if (auto) max_ncomp[1] else k_fixed, ncol(Fm), length(train) - 2)
    dtr <- data.frame(Y = I(Ys[train, , drop = FALSE]), Fm = I(Fm[train, , drop = FALSE]))
    dte <- data.frame(Fm = I(Fm[test, , drop = FALSE]))
    fit <- pls::plsr(Y ~ Fm, ncomp = kmax, data = dtr, scale = scale,
                     validation = if (auto) "CV" else "none",
                     segments = min(cv_segments, length(train)))
    nc <- kmax
    if (auto) {
      cv <- pls::MSEP(fit, estimate = "CV", intercept = FALSE)$val   # [1, response, component]
      nc <- as.integer(which.min(apply(cv[1, , , drop = FALSE], 3, sum)))
    }
    list(pred = matrix(stats::predict(fit, newdata = dte, ncomp = nc), nt, q), nc = as.integer(nc))
  }

  if (is.null(Vl)) {   # unconditional test: M0 is the training mean (0 on the Ys scale)
    f1 <- fit_pred(X, if (auto) NULL else ncomp[2])
    return(list(yh0 = matrix(0, nt, q), yh1 = f1$pred, a = 0, b = f1$nc))
  }
  Vm <- do.call(cbind, Vl)
  f0 <- if (isFALSE(opts$fit_m0)) NULL else fit_pred(Vm, if (auto) NULL else ncomp[1])
  f1 <- fit_pred(cbind(Vm, X), if (auto) NULL else ncomp[2])
  list(yh0 = f0$pred, yh1 = f1$pred, a = f0$nc %||% NA_integer_, b = f1$nc)
}

# ROSA, response-oriented sequential alternation (multiblock::rosa). M0: Y ~ V1 + ...
# + Vk; M1: Y ~ V1 + ... + Vk + X, with every block kept as a separate ROSA block (each
# component is taken from the block that best explains the current Y residual).
# opts: ncomp ("auto": chosen by inner CV for each model, or c(a, b) for M0 and M1),
#   scale, max_ncomp (maximum components of each model when ncomp = "auto";
#   default 10), cv_segments.
.learner_rosa <- function(Ys, X, Vl, train, test, opts) {
  ncomp <- opts$ncomp %||% "auto"
  scale <- opts$scale %||% TRUE
  max_ncomp <- opts$max_ncomp %||% 10
  cv_segments <- opts$cv_segments %||% 5
  q <- ncol(Ys); nt <- length(test); ntr <- length(train)
  auto <- identical(ncomp, "auto")
  if (!auto && length(ncomp) != 2) stop("for ROSA, ncomp must be \"auto\" or c(a, b): components for Y ~ V and for Y ~ V + X")
  Ytr <- Ys[train, , drop = FALSE]

  fit_pred <- function(blocks, k_fixed) {
    bn <- paste0("B", seq_along(blocks))
    # standardise every block column with the TRAINING means (and sds if scale), and
    # fit rosa without its own scaling; predictions are then computed directly from
    # the coefficients (not with predict()), so the centring and scaling of the test
    # rows are under our control
    Xall <- do.call(cbind, blocks)
    blk <- rep(seq_along(blocks), vapply(blocks, ncol, 1L))
    mu <- colMeans(Xall[train, , drop = FALSE])
    sdv <- if (scale) apply(Xall[train, , drop = FALSE], 2, stats::sd) else rep(1, ncol(Xall))
    sdv[!is.finite(sdv) | sdv == 0] <- 1
    Xs <- sweep(sweep(Xall, 2, mu), 2, sdv, "/")
    kmax <- min(if (auto) max_ncomp[1] else k_fixed, ncol(Xall), ntr - 2)
    dtr <- c(list(Y = Ytr), stats::setNames(lapply(seq_along(blocks), function(b)
      Xs[train, blk == b, drop = FALSE]), bn))
    fml <- stats::as.formula(paste("Y ~", paste(bn, collapse = " + ")))
    cv_args <- if (auto) list(validation = "CV", segments = min(cv_segments, ntr)) else list()
    fit <- do.call(multiblock::rosa, c(list(fml, data = dtr, ncomp = kmax, scale = FALSE), cv_args))
    nc <- min(kmax, fit$ncomp)
    if (auto) {   # component count with the lowest CV error (summed over Y columns)
      cvp <- fit$validation$pred                                   # ntr x q x ncomp
      nc <- which.min(vapply(seq_len(dim(cvp)[3]),
                             function(a) sum((Ytr - matrix(cvp[, , a], ntr))^2), 0))
    }
    # prediction = training mean of Y + (centred test rows) %*% coefficients
    B <- matrix(stats::coef(fit, ncomp = nc, intercept = FALSE), ncol(Xall), q)
    xm <- colMeans(Xs[train, , drop = FALSE])                      # 0 up to rounding
    pred <- sweep(Xs[test, , drop = FALSE], 2, xm) %*% B
    pred <- sweep(pred, 2, colMeans(Ytr), "+")
    list(pred = matrix(pred, nt, q), nc = as.integer(nc))
  }

  if (is.null(Vl)) {   # unconditional test: M0 is the training mean (0 on the Ys scale)
    f1 <- fit_pred(list(X), if (auto) NULL else ncomp[2])
    return(list(yh0 = matrix(0, nt, q), yh1 = f1$pred, a = 0, b = f1$nc))
  }
  f0 <- if (isFALSE(opts$fit_m0)) NULL else fit_pred(Vl, if (auto) NULL else ncomp[1])
  f1 <- fit_pred(c(Vl, list(X)), if (auto) NULL else ncomp[2])
  list(yh0 = f0$pred, yh1 = f1$pred, a = f0$nc %||% NA_integer_, b = f1$nc)
}

# Multi-target XGBoost (xgboost >= 2.0, multi_strategy = "multi_output_tree").
# M0: Y ~ V; M1: Y ~ V + X (cbind-ed). Each model is trained with a fixed number of
# rounds on ALL the training rows (no early stopping) and only predicts the test rows.
# opts$args (all optional): nrounds (default 1000) and params (list merged into the
#   defaults: max_depth = 3, eta = 0.1, tree_method = "hist",
#   multi_strategy = "multi_output_tree", nthread = 1).
.learner_xgboost <- function(Ys, X, Vl, train, test, opts) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("The package 'xgboost' (>= 2.0) is required for method = \"XGB\".")
  q <- ncol(Ys)
  nt <- length(test)
  la <- utils::modifyList(list(nrounds = 2000, params = list()), opts$args %||% list())
  params <- utils::modifyList(
    list(objective = "reg:squarederror", tree_method = "hist", multi_strategy = "multi_output_tree",
         max_depth = 3, eta = 0.1, nthread = 1),
    la$params)

  fit_pred <- function(F) {
    dtr <- xgboost::xgb.DMatrix(F[train, , drop = FALSE], label = Ys[train, , drop = FALSE])
    fit <- xgboost::xgb.train(params = params, data = dtr, nrounds = la$nrounds, verbose = 0)
    matrix(stats::predict(fit, F[test, , drop = FALSE]), nt, q)
  }

  if (is.null(Vl)) {
    return(list(yh0 = matrix(0, nt, q), yh1 = fit_pred(X), a = 0, b = la$nrounds))
  }
  Vm <- do.call(cbind, Vl)
  skip0 <- isFALSE(opts$fit_m0)
  list(yh0 = if (skip0) NULL else fit_pred(Vm), yh1 = fit_pred(cbind(Vm, X)),
       a = if (skip0) NA else la$nrounds, b = la$nrounds)
}

# Multi-output random forest (grf::multi_regression_forest). M0: Y ~ V; M1: Y ~ V + X
# (cbind-ed). Each forest is grown on ALL the training rows, with splits chosen to
# reduce the error summed over the Y columns, and only predicts the test rows.
# opts$args (all optional) are passed to grf::multi_regression_forest(); defaults here
#   are num.trees = 2000, num.threads = 1 (so parallel simulations do not
#   oversubscribe the cores; NULL uses all cores) and compute.oob.predictions = FALSE.
#   The other grf defaults apply (honesty = TRUE, sample.fraction = 0.5,
#   min.node.size = 5, ...). The forest seed is drawn from R's RNG, so set.seed()
#   makes the result reproducible.
.learner_rf <- function(Ys, X, Vl, train, test, opts) {
  if (!requireNamespace("grf", quietly = TRUE))
    stop("The package 'grf' is required for method = \"RF\".")
  q <- ncol(Ys)
  nt <- length(test)
  la <- utils::modifyList(list(num.trees = 2000, num.threads = 1, compute.oob.predictions = FALSE),
                          opts$args %||% list())

  fit_pred <- function(F) {
    F <- as.matrix(F)
    fit <- do.call(grf::multi_regression_forest,
                   c(list(X = F[train, , drop = FALSE], Y = Ys[train, , drop = FALSE]), la))
    matrix(stats::predict(fit, newdata = F[test, , drop = FALSE])$predictions, nt, q)
  }

  if (is.null(Vl)) {
    return(list(yh0 = matrix(0, nt, q), yh1 = fit_pred(X), a = 0, b = la$num.trees))
  }
  Vm <- do.call(cbind, Vl)
  skip0 <- isFALSE(opts$fit_m0)
  list(yh0 = if (skip0) NULL else fit_pred(Vm), yh1 = fit_pred(cbind(Vm, X)),
       a = if (skip0) NA else la$num.trees, b = la$num.trees)
}

# Multi-output k nearest neighbour regression (see knn_regression.R). M0: Y ~ V; M1:
# Y ~ V + X (cbind-ed, features standardised with the training rows if scale = TRUE).
# opts$args (optional): k ("auto": chosen for each model by cross-validation within
#   the training rows with cv_segments segments; or one number used for both models)
#   and k_grid (candidates for "auto").
.learner_knn <- function(Ys, X, Vl, train, test, opts) {
  la <- .knn_args(opts$args)
  fit_pred <- function(Fm)
    .knn_fit_predict(Fm, Ys, train, test, k = la$k, grid = la$grid,
                     scale = opts$scale %||% TRUE, segments = opts$cv_segments %||% 5)
  if (is.null(Vl)) {   # unconditional test: M0 is the training mean (0 on the Ys scale)
    f1 <- fit_pred(X)
    return(list(yh0 = matrix(0, length(test), ncol(Ys)), yh1 = f1$pred, a = 0, b = f1$k))
  }
  Vm <- do.call(cbind, Vl)
  f0 <- if (isFALSE(opts$fit_m0)) NULL else fit_pred(Vm)
  f1 <- fit_pred(cbind(Vm, X))
  list(yh0 = f0$pred, yh1 = f1$pred, a = f0$k %||% NA_integer_, b = f1$k)
}

# ---- registry -----------------------------------------------------------------
# Methods available by name in `method`. Add new learners here.
.cci_methods <- list(
  SOPLS = .learner_sopls,
  PLS   = .learner_pls,
  ROSA  = .learner_rosa,
  XGB   = .learner_xgboost,
  RF    = .learner_rf,
  KNN   = .learner_knn
)

# Look up a learner: `method` is a name in `registry` (case-insensitive) or a
# learner function with the interface of that registry.
.get_learner <- function(method, registry = .cci_methods) {
  if (is.function(method)) return(method)
  if (!is.character(method) || length(method) != 1)
    stop("`method` must be one of ", paste(names(registry), collapse = ", "),
         ", or a learner function")
  key <- toupper(method)
  if (!key %in% names(registry))
    stop("unknown method '", method, "'; available: ",
         paste(names(registry), collapse = ", "), ", or a learner function")
  registry[[key]]
}

# Display name of a method (the registry name, or "custom" for a function)
.method_name <- function(method) if (is.function(method)) "custom" else toupper(method)
