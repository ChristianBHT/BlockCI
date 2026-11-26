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
#' * `"PLS"` (default): ordinary PLS regression (`pls::plsr`) of Y on cbind(V) and
#'   on cbind(V, X), so the block structure is not used. With `ncomp = "auto"` the
#'   number of components of each model is chosen by inner CV.
#' * `"SOPLS"`: both models come from one `multiblock::sopls` fit (M0: the
#'   V components; M1: after adding b components of X, fitted sequentially on the
#'   V-orthogonalised Y and X). Each V block is a separate SO-PLS block. With
#'   `ncomp = "auto"` each V block gets one component and the number of X
#'   components b is chosen by inner CV; b = 0 means X is not selected, so
#'   M1 = M0 and T = 0.
#' * `"ROSA"`: response-oriented sequential alternation (`multiblock::rosa`) of Y on
#'   the V blocks and on the V blocks plus X, every block kept separate (each
#'   component comes from the block that best explains the current Y residual).
#'   With `ncomp = "auto"` the number of components of each model is chosen by
#'   inner CV.
#' * `"XGB"` (experimental, needs more work: type I error is not controlled on
#'   nonlinear scenarios, see `LOG.md`): two multi-target XGBoost models (needs
#'   xgboost >= 2.0), each trained with a fixed number of boosting rounds on all
#'   the training rows (no early stopping).
#' * `"RF"`: two multi-output random forests (`grf::multi_regression_forest`, needs
#'   the grf package), each grown on all the training rows with splits that reduce
#'   the error summed over the Y columns.
#' * `"KNN"`: multi-output k nearest neighbour regression (needs FNN) of Y on
#'   cbind(V) and on cbind(V, X): the prediction is the mean Y of the k nearest
#'   training rows (Euclidean distance, features standardised with the training rows
#'   if `scale = TRUE`). k is chosen for each model by cross-validation within the
#'   training rows (`cv_segments` segments), or fixed with `method_args`. All columns
#'   count equally in the distance, so many uninformative columns weaken it.
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
#' @param method regression method: `"PLS"` (default), `"SOPLS"`, `"ROSA"`, `"XGB"`,
#'   `"RF"`, `"KNN"` (case-insensitive) or a learner function (see Details).
#' @param p training fraction (default 0.5)
#' @param ncomp (`"SOPLS"`, `"PLS"`, `"ROSA"`) `"auto"` (chosen by inner CV as described
#'   above) or fixed numbers of components. For `"SOPLS"`: `c(a_1, ..., a_k, b)`,
#'   the components for each of the k V blocks, then the components added from X
#'   (for a test without V, `c(0, b)`). For `"PLS"` and `"ROSA"`: `c(a, b)`, the
#'   components of Y ~ V and of Y ~ V + X.
#' @param scale scale the blocks (`"SOPLS"`, `"PLS"`, `"ROSA"`) and standardise Y
#' @param max_ncomp (`"SOPLS"`, `"PLS"`, `"ROSA"`) maximum number of components
#'   searched when `ncomp = "auto"`: for `"SOPLS"` the maximum number of X
#'   components (default 5), for `"PLS"` and `"ROSA"` the maximum for each of the
#'   two models (default 10). `NULL` uses the method's default.
#' @param cv_segments (`"SOPLS"`, `"PLS"`, `"ROSA"`) number of CV segments for the inner
#'   selection
#' @param method_args named list of method-specific options. For `"XGB"`:
#'   `nrounds` (number of boosting rounds, 1000) and `params` (a list merged into
#'   the xgboost parameters, defaults `max_depth = 3, eta = 0.1,
#'   tree_method = "hist", multi_strategy = "multi_output_tree", nthread = 1`).
#'   For `"RF"`: any arguments of [grf::multi_regression_forest()], e.g.
#'   `num.trees` (2000), `num.threads` (1; `NULL` uses all cores),
#'   `min.node.size`, `mtry`, `sample.fraction` or `honesty`. For `"KNN"`: `k`
#'   (`"auto"`, default, or one number of neighbours for both models) and `k_grid`
#'   (candidates for `"auto"`; default 3, 5, 7, 10, 15, 20, 30, 50, 75, 100). A
#'   learner function receives the list as `opts$args`.
#' @param split training rows (optional; if `NULL` a random split of `floor(p * n)`
#'   rows is drawn). [Block.CCI.fast()] uses it to fit X and X* on the same split.
#' @param reduce dimension reduction of the conditioning blocks before the
#'   regression: `"none"` (default), `"pca"` (principal component scores of all V
#'   blocks together, standardised) or `"pls"` (scores of a PLS regression of Y on V).
#'   The reduction is fitted on the training rows only, and the V blocks are replaced
#'   by one block of scores for both M0 and M1 (X is not reduced). Useful when V has
#'   many columns and a nonlinear learner ("RF", "XGB") has to learn functions of a few
#'   directions of V.
#' @param reduce_ncomp number of scores: `"auto"` or a number. For `"pca"`, `"auto"`
#'   keeps the components whose eigenvalue (correlation matrix of V) exceeds the
#'   Marchenko-Pastur noise edge `(1 + sqrt(p/n))^2` (at most 20; if none does, the
#'   components explaining 90 % of the variance); for `"pls"` the number is chosen by
#'   cross-validation on the training rows (at most 10).
#' @param losses if `TRUE`, return per-test-row squared errors as well
#' @param yh0 (advanced) test-row predictions of M0 (Y ~ V) from an earlier call
#'   with the same `data` (apart from X), `split` and settings, as returned with
#'   `losses = TRUE`. M0 does not depend on X, so it is then not refitted: only M1
#'   is fitted. [Block.CCI.fast()] uses this for the resampled X*.
#' @return named vector `M0, M1, T = M1 - M0, a, b` (RMSEP scale). `a` and `b` are
#'   the complexities of M0 and M1: for `"SOPLS"` the total number of V components
#'   and the number of added X components; for `"PLS"` and `"ROSA"` the components
#'   of Y ~ V and of Y ~ V + X; for `"XGB"` the number of boosting rounds and for `"RF"` the
#'   number of trees, and for `"KNN"` the number of neighbours k (`a = 0` without V). With `losses = TRUE`: a list with `stats`
#'   (that vector), `test` (the test rows), `l0` and `l1` (per-test-row squared
#'   errors, summed over the Y columns, of M0 and M1) and `yh0` (the test-row
#'   predictions of M0, standardised Y scale).
#' @export
cci_test_gen <- function(formula, data, method = "PLS", p = 0.5, ncomp = "auto",
                         scale = TRUE, max_ncomp = NULL, cv_segments = 5,
                         method_args = list(), split = NULL, losses = FALSE,
                         yh0 = NULL, reduce = c("none", "pca", "pls"), reduce_ncomp = "auto") {
  learner <- .get_learner(method)
  reduce <- match.arg(reduce)
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

  # optional dimension reduction of V (fitted on the training rows)
  V <- .reduce_blocks(V, Ys, train, reduce, reduce_ncomp)

  opts <- list(ncomp = ncomp, scale = scale, max_ncomp = max_ncomp,
               cv_segments = cv_segments, args = method_args, fit_m0 = is.null(yh0))
  fit <- learner(Ys, X, V, train, test, opts)
  if (!is.null(yh0)) {   # reuse the M0 predictions of an earlier call on the same split
    yh0 <- as.matrix(yh0)
    if (nrow(yh0) != length(test) || ncol(yh0) != ncol(Ys))
      stop("`yh0` does not match the test rows of this split")
    fit$yh0 <- yh0
    fit$a <- NA
  }

  Yt <- Ys[test, , drop = FALSE]
  M0 <- rmsep(Yt, fit$yh0)
  M1 <- rmsep(Yt, fit$yh1)
  out <- c(M0 = M0, M1 = M1, T = M1 - M0, a = fit$a, b = fit$b)
  if (!losses) return(out)
  # per-test-row squared errors (summed over Y columns) for M0 and M1
  list(stats = out, test = test,
       l0 = rowSums((Yt - fit$yh0)^2), l1 = rowSums((Yt - fit$yh1)^2),
       yh0 = fit$yh0)
}

# Two-fold cross-fitted version of the statistic (crossfit = TRUE in Block.CCI.test()
# and cci_null_dist()). The rows are split at random into two halves; M0 and M1 are
# fitted on each half and predict the other (everything, tuning and standardisation
# included, uses the fitting half only, as in cci_test_gen()). The squared errors of
# all rows are pooled into ONE RMSEP for each model. With an odd number of rows the
# halves differ by one row. Takes the arguments of cci_test_gen() through `...`
# (`p` is accepted so that both can be called alike; the halves are 50/50). Returns the
# same vector as cci_test_gen(): M0, M1, T, and a, b averaged over the two folds.
.cci_crossfit <- function(formula, data, ..., p = 0.5) {
  cc <- .canonical(formula, data)            # incomplete rows are dropped before halving
  n <- nrow(cc$Y)
  half <- sample.int(n, n %/% 2)
  folds <- lapply(list(half, setdiff(seq_len(n), half)), function(train)
    cci_test_gen(cc$formula, cc$data, ..., split = train, losses = TRUE))
  pooled <- function(loss)
    sqrt(sum(vapply(folds, function(f) sum(f[[loss]]), 0)) / (n * ncol(cc$Y)))
  M0 <- pooled("l0"); M1 <- pooled("l1")
  complexity <- rowMeans(vapply(folds, function(f) f$stats[c("a", "b")], numeric(2)))
  c(M0 = M0, M1 = M1, T = M1 - M0, complexity)
}

.check_crossfit <- function(crossfit, p) {
  if (!isTRUE(crossfit) && !isFALSE(crossfit)) stop("`crossfit` must be TRUE or FALSE")
  if (crossfit && !isTRUE(all.equal(p, 0.5)))
    stop("crossfit = TRUE uses two halves of the rows: `p` must be 0.5")
}
