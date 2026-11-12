#' Test of conditional independence Y _||_ X | Z via cross-fitted XGBoost residuals
#' (experimental, needs more work)
#'
#' **Status: needs more work.** In simulations (`inst/benchmarks/multixg_power.R`,
#' see `LOG.md`) this test rejects a true null far too often (all of the `sim_1` and
#' `sim_5` runs at n = 500): the boosted trees do not remove the effect of Z from X
#' and Y well enough, so the residuals stay dependent. Do not use it for real
#' decisions yet.
#'
#' Same logic as [pls_ci_test()] and [rosa_ci_test()]: remove Z from X and from Y
#' with a regression, then test with a permutation test whether the two residual
#' matrices are still dependent. Here the regression is multi-target XGBoost
#' (`multi_strategy = "multi_output_tree"`, xgboost >= 2.0), so nonlinear effects
#' of Z are removed. The predictions are made by *cross-fitting*: the rows are split
#' in `nfolds` folds, every fold is predicted by a model trained on the other folds,
#' so no row is predicted by a model that has seen it. The residuals are then
#' `M - cross-fitted prediction` for M = X and M = Y.
#'
#' In each training part the targets are standardised (training mean and sd), and
#' the number of boosting rounds is chosen by early stopping on a validation part
#' (`val_fraction`) of that training part, so the held-out fold is never used for
#' tuning.
#'
#' @param formula `Y ~ X | Z` (the conditioning set is required). Several blocks
#'   may be given after the bar, e.g. `Y ~ X | Z1 + Z2`; they are combined
#'   column-wise.
#' @param data named list of matrices (or a data.frame) holding the variables in
#'   `formula`. Rows with missing values are dropped.
#' @param nfolds number of cross-fitting folds.
#' @param nrounds maximum number of boosting rounds.
#' @param early_stopping_rounds early stopping patience.
#' @param val_fraction share of each training part used as early-stopping validation.
#' @param params list merged into the xgboost parameters (defaults:
#'   `max_depth = 3, eta = 0.1, tree_method = "hist",
#'   multi_strategy = "multi_output_tree", nthread = 1`).
#' @param stat `"rv"`, `"cc1"` or `"dcor"` (the last needs the `energy` package).
#' @param perm `"conditional"` (permute residual rows within groups of
#'   observations with similar Z) or `"simple"` (permute all rows).
#' @param B number of permutations.
#' @param group_size desired group size for conditional permutation.
#' @param seed optional seed.
#' @return an object of class `"htest"` with the mean number of boosting rounds in
#'   `rounds` and the permutation statistics in `perm_stats`.
#' @export
multixg_ci_test <- function(formula, data,
                            nfolds = 5, nrounds = 300, early_stopping_rounds = 10,
                            val_fraction = 0.2, params = list(),
                            stat = c("rv", "cc1", "dcor"),
                            perm = c("conditional", "simple"),
                            B = 999, group_size = 10, seed = NULL) {
  stat <- match.arg(stat)
  perm <- match.arg(perm)
  if (!is.null(seed)) set.seed(seed)
  .message_once("multixg_ci_test() is experimental and needs more work: it can reject true nulls far too often (see ?multixg_ci_test).")

  cc <- .canonical(formula, data)            # parses formula, drops incomplete rows
  if (is.null(cc$V)) stop("formula must have a conditioning set: Y ~ X | Z")
  Z <- do.call(cbind, cc$V)
  rs <- .multixg_residuals(cc$X, cc$Y, Z, nfolds, nrounds, early_stopping_rounds,
                           val_fraction, params)
  res <- .resid_perm_test(rs$res_x, rs$res_y, Z, stat, perm, B, group_size)

  structure(list(
    statistic   = c(stat = res$obs),
    p.value     = res$p,
    method      = sprintf("XGBoost cross-fitted residual test for X _||_ Y | Z (%s, %s permutation, %d folds)",
                          stat, perm, nfolds),
    data.name   = paste(deparse(formula), collapse = ""),
    null.value  = NULL,
    alternative = "X and Y are dependent given Z",
    rounds      = c(X = rs$rounds_x, Y = rs$rounds_y),
    perm_stats  = res$perms
  ), class = "htest")
}

# Cross-fitted XGBoost residuals of X and of Y given Z.
.multixg_residuals <- function(X, Y, Z, nfolds, nrounds, early_stopping_rounds,
                               val_fraction, params) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("The package 'xgboost' (>= 2.0) is required for multixg_ci_test().")
  n <- nrow(X)
  if (nfolds < 2 || nfolds > n) stop("nfolds must be between 2 and the number of rows")
  prm <- utils::modifyList(
    list(objective = "reg:squarederror", tree_method = "hist", multi_strategy = "multi_output_tree",
         max_depth = 3, eta = 0.1, nthread = 1),
    params)
  folds <- sample(rep_len(seq_len(nfolds), n))      # shared by the X and the Y regression

  crossfit <- function(M) {
    q <- ncol(M)
    res <- matrix(NA_real_, n, q)
    rounds <- numeric(nfolds)
    for (k in seq_len(nfolds)) {
      te <- which(folds == k)
      tr <- which(folds != k)
      mu <- colMeans(M[tr, , drop = FALSE])
      sdv <- apply(M[tr, , drop = FALSE], 2, stats::sd); sdv[sdv == 0] <- 1
      Ms <- sweep(sweep(M, 2, mu), 2, sdv, "/")      # targets standardised with training stats
      nv <- max(2L, ceiling(val_fraction * length(tr)))
      va <- sample(tr, nv)                           # early-stopping rows (inside the training part)
      tr2 <- setdiff(tr, va)
      dtr <- xgboost::xgb.DMatrix(Z[tr2, , drop = FALSE], label = Ms[tr2, , drop = FALSE])
      dva <- xgboost::xgb.DMatrix(Z[va, , drop = FALSE],  label = Ms[va, , drop = FALSE])
      fit <- xgboost::xgb.train(params = prm, data = dtr, nrounds = nrounds,
                                evals = list(val = dva),
                                early_stopping_rounds = early_stopping_rounds, verbose = 0)
      best <- xgboost::xgb.attr(fit, "best_iteration")
      rounds[k] <- if (is.null(best)) xgboost::xgb.get.num.boosted.rounds(fit) else as.integer(best) + 1L
      pr <- matrix(stats::predict(fit, Z[te, , drop = FALSE]), length(te), q)
      res[te, ] <- M[te, , drop = FALSE] - sweep(sweep(pr, 2, sdv, "*"), 2, mu, "+")
    }
    list(res = res, rounds = mean(rounds))
  }

  rx <- crossfit(X)
  ry <- crossfit(Y)
  list(res_x = rx$res, res_y = ry$res, rounds_x = rx$rounds, rounds_y = ry$rounds)
}
