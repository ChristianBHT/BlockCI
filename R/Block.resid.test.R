#' Block residual test of conditional independence
#'
#' Tests whether the block Y is independent of the block X given the blocks Z
#' (`Y ~ X | Z1 + Z2 + ...`) with residuals: Z is removed from X and from Y with a
#' regression, and a permutation test checks whether the two residual matrices are
#' still dependent. Under the null, and if the regression removes the effect of Z
#' well enough, the residuals of X and Y are independent.
#'
#' The regression is chosen with `method`:
#'
#' * `"PLS"` (default): PLS regression (`pls::plsr`) on the `cbind`-ed Z blocks.
#' * `"ROSA"`: response-oriented sequential alternation (`multiblock::rosa`), with
#'   the Z blocks kept as separate blocks.
#' * `"XGB"`: multi-target XGBoost (needs xgboost >= 2.0).
#' * `"RF"`: multi-output random forest (`grf::multi_regression_forest`, needs grf).
#' * `"KNN"`: multi-output k nearest neighbour regression (needs FNN), k chosen by
#'   cross-validation within the training rows of every fold.
#'
#' For `"XGB"`, `"RF"` and `"KNN"` the linear part of the effect of Z is removed first
#' (`linear_first = TRUE`, default): in each fold a PLS regression of the target on
#' Z is fitted on the training rows, the tree learner is fitted on the training
#' residuals of that PLS fit, and the held-out rows are predicted by the sum of the
#' two. Tree ensembles approximate linear effects poorly (piecewise constant, pulled
#' towards the mean), so without this step a linear effect of Z is left partly in
#' both residual matrices and makes them dependent even when the null is true. The
#' PLS step uses `ncomp` and `max_ncomp` (with `ncomp = "auto"`, components chosen
#' by 5-fold CV within the training rows). With `linear_first = FALSE` the tree
#' learner is used alone; for `"RF"` the residuals then come from out-of-bag
#' predictions unless `method_args = list(oob = FALSE)`.
#' * a residual learner function with the interface described in
#'   `R/resid_learners.R`. New methods can be made available by name by adding them
#'   to `.resid_methods` in that file.
#'
#' **Dimension reduction (`reduce`).** With many conditioning variables a tree
#' learner may not find the few directions of Z that matter. `reduce = "pca"`
#' replaces the Z blocks by their principal component scores (one block; fitted on
#' all rows, which is allowed since PCA does not use X or Y) for the regression, the
#' conditional permutation and the wGCM weights. `reduce = "pls"` (`"XGB"`, `"RF"` and
#' `"KNN"` only) feeds the tree learner with the PLS scores of the linear step of each fold
#' (fitted on the training rows, separately for X and Y) instead of Z itself; the
#' linear step then removes the linear part as before. For `"PLS"` and `"ROSA"`,
#' which already work on a few components, `reduce = "pls"` is ignored.
#' A nonlinear effect of Z that is shared by X and Y and lies outside the retained
#' directions is not removed, and the test then rejects also under the null.
#'
#' **Weighted GCM (`stat = "wgcm"`).** The multivariate weighted generalised
#' covariance measure with fixed weights (mWGCM.fix; Scheidegger, Hoerrmann and
#' Buehlmann 2022, JMLR 23, Appendix A.1). For every pair of a residual column of X
#' (eps_j) and of Y (xi_l) and every weight function w_k(Z) it computes
#' `T_jlk = sqrt(n) mean(eps_j xi_l w_k) / sd(eps_j xi_l w_k)`, and the statistic is
#' `max |T_jlk|`. The weights (Appendix A.1.1) are the constant 1 and
#' `sign(Z_d - a_dk)` for every column d of Z, with `a_dk` the empirical
#' `k / (k0 + 1)` quantiles, `k = 1, ..., k0`. Weights that depend on Z let the test
#' detect conditional dependence that changes sign across Z and averages to zero,
#' which an unweighted covariance misses. Under the null the vector of `T_jlk` is
#' approximately normal with the correlation matrix of the weighted products; the
#' p-value is the share of `B` draws of `max |N(0, Sigma_hat)|` that reach the
#' observed maximum. `perm` and `k_perm` are not used. There are
#' `ncol(X) * ncol(Y) * (k0 * ncol(Z) + 1)` statistics; time and memory grow with
#' that number times n.
#'
#' With `resid = "cv"` (recommended) every row's residual comes from a model that
#' did not use that row: the rows are split into `nfolds` folds, shared by the X and
#' the Y regression. For PLS and ROSA the number of components is chosen by the CV
#' error over the same folds (`ncomp = "auto"`).
#'
#' @param formula `Y ~ X | Z` (the conditioning set is required). Several blocks may
#'   be given after the bar, e.g. `Y ~ X | Z1 + Z2`.
#' @param data named list of matrices (or a data.frame) holding the variables in
#'   `formula`. Rows with missing values are dropped.
#' @param method regression used to remove Z: `"PLS"` (default), `"ROSA"`, `"XGB"`,
#'   `"RF"`, `"KNN"` (case-insensitive) or a residual learner function (see Details).
#' @param stat dependence statistic between the residual matrices: `"rv"` (RV
#'   coefficient), `"cc1"` (first canonical correlation), `"dcor"` (distance
#'   correlation; needs the energy package) or `"wgcm"` (weighted generalised
#'   covariance measure, see Details). `"rv"`, `"cc1"` and `"dcor"` are calibrated
#'   by permutation, `"wgcm"` by a Gaussian approximation.
#' @param resid `"cv"` (cross-validated / cross-fitted residuals, recommended) or
#'   `"fitted"` (in-sample residuals).
#' @param perm `"conditional"` (greedy local permutation of the residual rows: each
#'   row gets the Y residual of a row with similar Z, see [permute_X_cond()]) or `"simple"` (permute all rows).
#' @param B number of permutations (for `stat = "wgcm"`: number of draws from the
#'   Gaussian approximation).
#' @param k_perm neighbourhood size (the row itself included) for
#'   `perm = "conditional"` (default 5).
#' @param k0 (`stat = "wgcm"`) number of quantile-based weight functions per column
#'   of Z (default 7, as in Scheidegger et al. 2022).
#' @param nfolds number of folds for `resid = "cv"` and for choosing the number of
#'   components.
#' @param ncomp number of PLS or ROSA components, or `"auto"` (chosen by CV,
#'   separately for X and Y). Also used by the linear step of `"XGB"`, `"RF"` and `"KNN"`.
#' @param max_ncomp upper limit for `ncomp` when `ncomp = "auto"`.
#' @param scale scale the Z variables in the PLS / ROSA regression.
#' @param linear_first (`"XGB"`, `"RF"`, `"KNN"`) remove the linear part of Z with PLS before
#'   the tree learner (default `TRUE`); ignored for `"PLS"` and `"ROSA"`.
#' @param reduce dimension reduction of Z before the regression: `"none"` (default),
#'   `"pca"` or `"pls"` (see Details).
#' @param reduce_ncomp number of scores: `"auto"` or a number. For `"pca"`, `"auto"`
#'   keeps the components above the Marchenko-Pastur noise edge (see
#'   [cci_test_gen()]); for `"pls"`, the number of components of the linear step
#'   (chosen by CV with `ncomp = "auto"`).
#' @param method_args named list of method-specific options. For `"KNN"`: `k`
#'   (`"auto"`, default, or a number of neighbours) and `k_grid` (candidates for
#'   `"auto"`). For `"XGB"`: `nrounds`
#'   (number of boosting rounds, 1000) and `params` (list merged into
#'   `max_depth = 3, eta = 0.1, tree_method = "hist",
#'   multi_strategy = "multi_output_tree", nthread = 1`). For `"RF"`: `oob` (`TRUE`)
#'   and any arguments of [grf::multi_regression_forest()], e.g. `num.trees` (2000)
#'   or `num.threads` (1; `NULL` uses all cores). A learner function receives the
#'   list as `opts$args`.
#' @param seed optional random seed.
#' @return an object of class `c("BlockResid", "htest")` with the statistic, the
#'   p-value, the regression `regression`, the complexity of the X and the Y
#'   regression in `complexity` (components, rounds or trees), the mean number of
#'   PLS components of the linear step in `linear_ncomp` (XGB and RF with
#'   `linear_first = TRUE`, otherwise `NA`), the permutation statistics (for
#'   `"wgcm"`: the simulated maxima) in `perm_stats`, and for `"wgcm"` the number of
#'   statistics in the maximum in `n_stats`.
#' @examples
#' \dontrun{
#' data <- sim_1(n = 500, p = 4)
#' Block.resid.test(Y ~ X | Z1 + Z2, data)                      # H0 true, PLS
#' Block.resid.test(Y ~ X | Z1, data, method = "ROSA")          # H0 false
#' Block.resid.test(Y ~ X | Z1 + Z2, data, method = "RF", B = 199,   # PLS, then RF
#'                  method_args = list(num.trees = 1000, num.threads = NULL))
#' Block.resid.test(Y ~ X | Z1 + Z2, data, method = "XGB", B = 199,
#'                  method_args = list(nrounds = 300))
#' Block.resid.test(Y ~ X | Z1 + Z2, data, method = "RF", linear_first = FALSE)  # RF alone
#' Block.resid.test(Y ~ X | Z1 + Z2, data, method = "KNN", B = 199)         # PLS, then KNN
#' Block.resid.test(Y ~ X | Z1 + Z2, data, method = "KNN", stat = "wgcm",
#'                  method_args = list(k = 20))
#' Block.resid.test(Y ~ X | Z1 + Z2, data, stat = "wgcm")                # weighted GCM
#' # many conditioning variables: RF on the PLS scores, or PLS on PCA scores
#' Block.resid.test(Y ~ X | Z1 + Z2, data, method = "RF", reduce = "pls")
#' Block.resid.test(Y ~ X | Z1 + Z2, data, method = "PLS", reduce = "pca")
#' }
#' @export
Block.resid.test <- function(formula, data, method = "PLS",
                             stat = c("rv", "cc1", "dcor", "wgcm"),
                             resid = c("cv", "fitted"),
                             perm = c("conditional", "simple"),
                             B = 999, k_perm = 5, nfolds = 10,
                             ncomp = "auto", max_ncomp = 10, scale = TRUE,
                             linear_first = TRUE, k0 = 7, method_args = list(), seed = NA,
                             reduce = c("none", "pca", "pls"), reduce_ncomp = "auto") {
  stat  <- match.arg(stat)
  reduce <- match.arg(reduce)
  resid <- match.arg(resid)
  perm  <- match.arg(perm)
  learner <- .get_learner(method, .resid_methods)
  method_name <- .method_name(method)
  if (!is.na(seed)) set.seed(seed)
  if (stat == "dcor" && !requireNamespace("energy", quietly = TRUE))
    stop("The package 'energy' is required for stat = 'dcor'.")

  cc <- .canonical(formula, data)          # parses the formula, drops incomplete rows
  if (is.null(cc$V)) stop("formula must have a conditioning set: Y ~ X | Z")
  n <- nrow(cc$X)
  if (nfolds < 2 || nfolds > n) stop("nfolds must be between 2 and the number of rows")
  folds <- sample(rep_len(seq_len(nfolds), n))   # shared by the X and the Y regression
  if (reduce == "pls" && method_name %in% c("PLS", "ROSA")) reduce <- "none"
  n_reduced <- NA_integer_
  if (reduce == "pca") {                    # unsupervised: fitted on all rows
    cc$V <- list(.reduce_pca(do.call(cbind, cc$V), seq_len(n), reduce_ncomp))
    n_reduced <- ncol(cc$V[[1]])
  }
  opts <- list(ncomp = ncomp, max_ncomp = max_ncomp, scale = scale, resid = resid,
               linear_first = linear_first, reduce = reduce, reduce_ncomp = reduce_ncomp,
               args = method_args)

  rx <- learner(cc$X, cc$V, folds, opts)
  ry <- learner(cc$Y, cc$V, folds, opts)
  Z <- do.call(cbind, cc$V)
  res <- .resid_perm_test(rx$res, ry$res, Z, stat, perm, B, k_perm, k0 = k0)
  calib <- if (stat == "wgcm") "Gaussian approximation" else paste(perm, "permutation")

  two_step <- linear_first && method_name %in% c("XGB", "RF", "KNN")
  reg_label <- if (two_step) paste0("PLS + ", method_name) else method_name
  if (reduce == "pca") reg_label <- paste0(reg_label, " on ", n_reduced, " PCA scores of Z")
  if (reduce == "pls") reg_label <- paste0(reg_label, " on PLS scores of Z")
  structure(list(
    statistic   = stats::setNames(res$obs, toupper(stat)),
    p.value     = res$p,
    method      = sprintf("Block residual test for Y _||_ X | Z (%s regression, %s, %s, %s residuals)",
                          reg_label, stat, calib, resid),
    data.name   = paste(deparse(formula), collapse = ""),
    alternative = "X and Y are dependent given Z",
    regression  = method_name,
    complexity  = c(X = rx$complexity, Y = ry$complexity),
    linear_ncomp = c(X = rx$lin_ncomp %||% NA_real_, Y = ry$lin_ncomp %||% NA_real_),
    linear_first = two_step,
    reduce = reduce, n_reduced = n_reduced,
    perm_stats  = res$perms,
    n_stats = res$n_stats,
    B = B, nfolds = nfolds, resid = resid, perm = if (stat == "wgcm") NA else perm, stat = stat
  ), class = c("BlockResid", "htest"))
}
