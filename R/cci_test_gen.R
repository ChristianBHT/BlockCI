#' One Monte Carlo CV statistic for the block CCI test
#'
#' Splits the rows into a training and a test part, fits two regressions on the
#' training rows, M0: Y ~ V and M1: Y ~ V + X, and compares their prediction errors
#' on the test rows: `T = M1 - M0` (RMSEP, with Y standardised by the training means
#' and sds). Negative `T` means that X improves the prediction of Y beyond V. This
#' is one draw of the statistic: [cci_null_dist()] repeats it on resampled copies of
#' X to build the null distribution, and [Block.CCI.test()] combines the two into a
#' p-value.
#'
#' The regression is chosen with `method`:
#'
#' * `"SOPLS"` (default): both models come from one `multiblock::sopls` fit (M0: the
#'   V components; M1: after adding b components of X, fitted sequentially on the
#'   V-orthogonalised Y and X). Each V block is a separate SO-PLS block. With
#'   `ncomp = "auto"` each V block gets one component and the number of X
#'   components b is chosen by inner CV; b = 0 means X is not selected, so
#'   M1 = M0 and T = 0.
#' * `"PLS"`: ordinary PLS regression (`pls::plsr`) of Y on cbind(V) and on
#'   cbind(V, X), so the block structure is not used. With `ncomp = "auto"` the
#'   number of components of each model is chosen by inner CV.
#' * `"XGB"` (experimental, needs more work: type I error is not controlled on
#'   nonlinear scenarios, see `LOG.md`): two multi-target XGBoost models (needs
#'   xgboost >= 2.0), each trained with a fixed number of boosting rounds on all
#'   the training rows (no early stopping).
#' * `"RF"`: two multi-output random forests (`grf::multi_regression_forest`, needs
#'   the grf package), each grown on all the training rows with splits that reduce
#'   the error summed over the Y columns.
#' * a learner function: any function with the interface described in
#'   `R/learners.R`, `function(Ys, X, Vl, train, test, opts)` returning
#'   `list(yh0, yh1, a, b)`. To make a new method available by name, add it to
#'   `.cci_methods` in that file.
#'
#' All tuning uses the training rows only, so the test rows are never used for
#' model selection, and the same procedure is applied to the observed X and to the
#' resampled X*.
#'
#' @param formula `Y ~ X | V` (V optional; several blocks allowed, e.g.
#'   `Y ~ X | V1 + V2`). Names refer to matrices or columns in `data`; each name
#'   after the bar is its own V block.
#' @param data named list of matrices (or a data.frame). Rows with missing values
#'   in any used variable are dropped.
#' @param method regression method: `"SOPLS"` (default), `"PLS"`, `"XGB"`, `"RF"`
#'   (case-insensitive) or a learner function (see Details).
#' @param p training fraction
#' @param ncomp (`"SOPLS"`, `"PLS"`) `"auto"` (chosen by inner CV as described
#'   above) or fixed numbers of components. For `"SOPLS"`: `c(a_1, ..., a_k, b)`,
#'   the components for each of the k V blocks, then the components added from X
#'   (for a test without V, `c(0, b)`). For `"PLS"`: `c(a, b)`, the components of
#'   Y ~ V and of Y ~ V + X.
#' @param scale scale the blocks (`"SOPLS"`, `"PLS"`) and standardise Y
#' @param max_ncomp (`"SOPLS"`, `"PLS"`) maximum number of components searched when
#'   `ncomp = "auto"`: for `"SOPLS"` the maximum number of X components (default
#'   5), for `"PLS"` the maximum for each of the two models (default 10). `NULL`
#'   uses the method's default.
#' @param cv_segments (`"SOPLS"`, `"PLS"`) number of CV segments for the inner
#'   selection
#' @param method_args named list of method-specific options. For `"XGB"`:
#'   `nrounds` (number of boosting rounds, 1000) and `params` (a list merged into
#'   the xgboost parameters, defaults `max_depth = 3, eta = 0.1,
#'   tree_method = "hist", multi_strategy = "multi_output_tree", nthread = 1`).
#'   For `"RF"`: any arguments of [grf::multi_regression_forest()], e.g.
#'   `num.trees` (2000), `num.threads` (1; `NULL` uses all cores),
#'   `min.node.size`, `mtry`, `sample.fraction` or `honesty`. A learner function
#'   receives the list as `opts$args`.
#' @param split training rows (optional; if `NULL` a random split of `floor(p * n)`
#'   rows is drawn). [Block.CCI.fast()] uses it to fit X and X* on the same split.
#' @param losses if `TRUE`, return per-test-row squared errors as well
#' @return named vector `M0, M1, T = M1 - M0, a, b` (RMSEP scale). `a` and `b` are
#'   the complexities of M0 and M1: for `"SOPLS"` the total number of V components
#'   and the number of added X components; for `"PLS"` the components of Y ~ V and
#'   of Y ~ V + X; for `"XGB"` the number of boosting rounds and for `"RF"` the
#'   number of trees (`a = 0` without V). With `losses = TRUE`: a list with `stats`
#'   (that vector), `test` (the test rows), `l0` and `l1` (per-test-row squared
#'   errors, summed over the Y columns, of M0 and M1).
#' @export
cci_test_gen <- function(formula, data, method = "SOPLS", p = 0.75, ncomp = "auto",
                         scale = TRUE, max_ncomp = NULL, cv_segments = 5,
                         method_args = list(), split = NULL, losses = FALSE) {
  learner <- .get_learner(method)
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

  opts <- list(ncomp = ncomp, scale = scale, max_ncomp = max_ncomp,
               cv_segments = cv_segments, args = method_args)
  fit <- learner(Ys, X, V, train, test, opts)

  Yt <- Ys[test, , drop = FALSE]
  M0 <- rmsep(Yt, fit$yh0)
  M1 <- rmsep(Yt, fit$yh1)
  out <- c(M0 = M0, M1 = M1, T = M1 - M0, a = fit$a, b = fit$b)
  if (!losses) return(out)
  # per-test-row squared errors (summed over Y columns) for M0 and M1
  list(stats = out, test = test,
       l0 = rowSums((Yt - fit$yh0)^2), l1 = rowSums((Yt - fit$yh1)^2))
}
