#' One Monte Carlo CV statistic for the CCI-type test
#'
#' M0: prediction error of Y from V on held-out rows. M1: prediction error of Y
#' from V and X. `T = M1 - M0` (RMSEP metric). The regression is done by a
#' *learner* wrapper chosen with `learner`:
#'
#' * `"sopls"` (default): both come from one `multiblock::sopls` fit (M0: V
#'   components; M1: after adding b components of X, sequentially on the
#'   V-orthogonalised Y and X). With `ncomp = "auto"` the number of X components b
#'   is chosen by inner CV on the TRAINING part only (test rows are never used for
#'   selection; identical for observed and null statistics); b = 0 means X is not
#'   selected, so M1 = M0 and T = 0.
#' * `"xgboost"`: two multi-target XGBoost models (Y ~ V and Y ~ V + X, needs
#'   xgboost >= 2.0); boosting rounds are chosen by early stopping on a validation
#'   part of the training rows. Options go in `learner_args`.
#'
#' Learners are defined in `R/learners.R`; a new learner is a function with the
#' same interface that is added to the `switch()` below.
#'
#' @param formula `Y ~ X | V` (V optional; several blocks allowed, e.g.
#'   `Y ~ X | V1 + V2`). Names refer to matrices or columns in `data`; each name
#'   after the bar is its own V block.
#' @param data named list of matrices (or a data.frame). Rows with missing values
#'   in any used variable are dropped.
#' @param p training fraction
#' @param ncomp (`learner = "sopls"`) `"auto"` (one component per V block, X
#'   components chosen by inner CV) or `c(a_1, ..., a_k, b)`: components for each of
#'   the k V blocks, then the components added from X (for a test without V,
#'   `c(0, b)`)
#' @param scale (`learner = "sopls"`) scale the blocks
#' @param max_ncomp (`learner = "sopls"`) maximum components searched when
#'   `ncomp = "auto"`
#' @param cv_segments (`learner = "sopls"`) CV segments for the inner selection
#' @param split Define the training rows (optional; if `NULL` a random split is drawn)
#' @param losses if `TRUE`, return per-test-row squared errors as well
#' @param learner `"sopls"` (default) or `"xgboost"`
#' @param learner_args named list of options for the learner. For `"xgboost"`:
#'   `nrounds` (maximum rounds, 300), `early_stopping_rounds` (10),
#'   `val_fraction` (inner validation share of the training rows, 0.2) and
#'   `params` (a list merged into the xgboost parameters, defaults
#'   `max_depth = 3, eta = 0.1, tree_method = "hist",
#'   multi_strategy = "multi_output_tree", nthread = 1`).
#' @return named vector `M0, M1, T = M1 - M0, a, b` (RMSEP metric; `a`, `b` are the
#'   complexities chosen by the learner: SO-PLS components, or XGBoost rounds of
#'   the V-only and V + X models); with `losses = TRUE` a list with `stats` (that
#'   vector), `test` (test rows), `l0`, `l1`
#' @export
sopls_test_gen <- function(formula, data, p = 0.75, ncomp = "auto", scale = TRUE,
                           max_ncomp = c(6, 10), cv_segments = 5,
                           split = NULL, losses = FALSE,
                           learner = c("sopls", "xgboost"), learner_args = list()) {
  learner <- match.arg(learner)
  cc <- .canonical(formula, data)
  Y <- cc$Y; X <- cc$X; V <- cc$V
  n <- nrow(Y)
  train <- if (is.null(split)) sample.int(n, floor(p * n)) else split
  test <- setdiff(seq_len(n), train)

  # standardise Y with training means and sds only
  mu <- colMeans(Y[train, , drop = FALSE])
  sdv <- if (scale) apply(Y[train, , drop = FALSE], 2, stats::sd) else rep(1, ncol(Y))
  sdv[sdv == 0] <- 1
  Ys <- sweep(sweep(Y, 2, mu), 2, sdv, "/")

  fit <- switch(learner,
                sopls   = .learner_sopls(Ys, X, V, train, test, ncomp, scale, max_ncomp, cv_segments),
                xgboost = .learner_xgboost(Ys, X, V, train, test, learner_args))

  Yt <- Ys[test, , drop = FALSE]
  M0 <- rmsep(Yt, fit$yh0)
  M1 <- rmsep(Yt, fit$yh1)
  out <- c(M0 = M0, M1 = M1, T = M1 - M0, a = fit$a, b = fit$b)
  if (!losses) return(out)
  # per-test-row squared errors (summed over Y columns) for the V-only and V + X models
  list(stats = out, test = test,
       l0 = rowSums((Yt - fit$yh0)^2), l1 = rowSums((Yt - fit$yh1)^2))
}
