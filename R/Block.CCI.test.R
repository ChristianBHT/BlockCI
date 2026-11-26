#' Block Computational Conditional Independence test
#'
#' Tests whether the block Y is independent of the block X given the blocks V
#' (`Y ~ X | V1 + V2 + ...`), with the Computational Conditional Independence (CCI)
#' approach: compare how well Y is predicted on held-out rows with and without X,
#' and calibrate that comparison against resampled copies of X.
#'
#' * **Statistic:** `T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V)` on the test rows of a
#'   random train/test split ([cci_test_gen()]). Negative `T` means that X improves
#'   the prediction of Y beyond V.
#' * **Null distribution:** `nperm` repetitions of the same procedure on a resampled
#'   X* ([cci_null_dist()]): rotation (default), permutation or conditional
#'   permutation of X, chosen with `resampling`.
#' * **p-value:** left tail, the share of null statistics at or below the observed
#'   one (`pvalue = "empirical"`), or from a normal fit to the null statistics
#'   (`pvalue = "parametric"`). By default the observed statistic comes from one
#'   random train/test split. With `npval > 1` (experimental) it is computed on
#'   `npval` splits, each is compared with the same null statistics and the p-values
#'   are averaged.
#'
#' The regression behind M0 (Y ~ V) and M1 (Y ~ V + X) is chosen with `method`:
#' 
#'* `"PLS"`(default): ordinary PLS (`pls::plsr`) on the `cbind`-ed blocks; the number of
#'   components of each model is chosen by inner CV.
#' * `"SOPLS"` : sequential and orthogonalised PLS (`multiblock::sopls`);
#'   each V block is a separate block. With `ncomp = "auto"` each V block gets one
#'   component and the number of X components is chosen by inner CV.
#' * `"ROSA"`: response-oriented sequential alternation (`multiblock::rosa`); every
#'   block is kept separate and each component comes from the block that best
#'   explains the current Y residual. Components chosen by inner CV.
#' * `"XGB"`: multi-target XGBoost (experimental; needs xgboost >= 2.0).
#' * `"RF"`: multi-output random forest (`grf::multi_regression_forest`; needs grf).
#' * `"KNN"`: multi-output k nearest neighbour regression (needs FNN); k chosen by
#'   inner CV.
#' * a learner function with the interface described in `R/learners.R`. New methods
#'   can be made available by name by adding them to `.cci_methods` in that file.
#'
#' With many conditioning variables, `reduce = "pca"` or `"pls"` replaces the V blocks
#' by a few scores before the regression (fitted on the training rows of every split),
#' so that a nonlinear learner works on a low-dimensional input; see [cci_test_gen()].
#' The resampling of X still adjusts for the full V.
#'
#' All tuning (components, ...) uses the training rows only, identically for the
#' observed and the resampled X. See [cci_test_gen()] for the details of each
#' method and its options.
#'
#' @inheritParams cci_null_dist
#' @param pvalue `"empirical"` (default): `(#{T* <= T_obs} + 1) / (R + 1)`.
#'   `"parametric"`: assumes the null distribution is normal with the mean and sd
#'   of the R null statistics and returns `pnorm(T_obs, mean, sd)`. This is not
#'   limited to a resolution of `1/(R + 1)`, but relies on the normal assumption.
#'   Both p-values are stored in the result as `p.empirical` and `p.parametric`.
#' @param npval number of observed statistics, each on its own random train/test
#'   split of the observed data; the p-value is the average of their `npval`
#'   p-values. A whole number from 1 to `nperm`; the default 1 is the test with a
#'   single observed split. **Warning:** values above 1 are experimental. An average
#'   of p-values computed from the same data is not uniformly distributed under the
#'   null. In linear Gaussian simulations (PLS, rotation, `nperm = 99`, level 0.05)
#'   the type I error fell from about 0.04 with `npval = 1` to about 0.01 with 2 and
#'   0.001 or less with 5, and the power at that level fell by more than half
#'   already with 2. Nonlinear data and other learners have not been tested. A
#'   larger `npval` does make the p-value depend less on the random split (but not
#'   on the randomness of the null draws). Each extra one costs one more statistic,
#'   so the number of statistics grows by the factor
#'   `(nperm + npval) / (nperm + 1)`.
#' @param crossfit (experimental, default `FALSE`) two-fold cross-fitting of the
#'   observed statistic **and** of every null draw, see [cci_null_dist()]: the models
#'   are fitted on each half of the rows and predict the other, and one statistic is
#'   computed from the errors of all rows. Needs `p = 0.5` and `npval = 1`, and about
#'   twice the model fits. In simulations with linear Gaussian data (PLS, n = 500,
#'   rotation, `nperm = 99`, level 0.05) the type I error was 0.039 (single 50/50
#'   split: 0.045) and the power at effect 0.10 / 0.15 was 0.30 / 0.68 (single split:
#'   0.24 / 0.57). A second run on new data (effect 0 / 0.10) gave 0.038 and 0.235
#'   (single split: 0.047 and 0.208), so the power gain at 0.10 may be smaller than
#'   in the first run. **Warning:** with a random forest (200 trees, n = 500) it rejected
#'   true nulls too often (0.13 with nonlinear relations, 0.10 with skewed noise;
#'   single split: 0.08 and 0.045). Other learners and settings have not been tested.
#' @param seed optional random seed
#' @return object of class `c("BlockCCI", "htest")`: a list with `statistic` (T),
#'   `p.value`, `M0`, `M1`, `p.empirical`, `p.parametric`, `null.distribution`,
#'   `method`, `resampling`, `complexity` (of M0 and M1, see [cci_test_gen()]) and
#'   the settings used, including `npval` and `crossfit` (with `crossfit = TRUE`,
#'   `statistic`, `M0` and `M1` are the pooled values and `complexity` the mean over
#'   the two folds). `statistic`, `M0`, `M1`, `complexity`,
#'   `p.empirical` and `p.parametric` are averages over the `npval` observed
#'   splits; the values of the single splits are in `T.splits`,
#'   `p.empirical.splits` and `p.parametric.splits`. [plot()] shows the null
#'   distribution.
#' @examples
#' \dontrun{
#' data <- sim_1(n = 500, p = 4)
#' Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99)                  # H0 true, PLS
#' Block.CCI.test(Y ~ X | Z1, data, nperm = 99, method = "PLS")      # H0 false
#' Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, crossfit = TRUE) # experimental
#' Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, method = "RF",
#'                method_args = list(num.trees = 500, num.threads = NULL))
#' Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, method = "KNN")
#' Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, method = "KNN",
#'                method_args = list(k = 10), resampling = "conditional")
#' Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 49, method = "XGB",
#'                method_args = list(nrounds = 300, params = list(max_depth = 4)))
#' }
#' @export
Block.CCI.test <- function(formula, data, method = "PLS",
                           resampling = c("rotation", "permutation", "conditional"),
                           nperm = 160, p = 0.5, ncomp = "auto", scale = TRUE,
                           max_ncomp = NULL, cv_segments = 5, method_args = list(),
                           pvalue = c("empirical", "parametric"),
                           k_perm = 5, seed = NA, progress = TRUE,
                           reduce = c("none", "pca", "pls"), reduce_ncomp = "auto",
                           npval = 1, crossfit = FALSE) {
  .check_crossfit(crossfit, p)
  resampling <- match.arg(resampling)
  reduce <- match.arg(reduce)
  pvalue <- match.arg(pvalue)
  if (!is.numeric(npval) || length(npval) != 1 ||
      !isTRUE(npval >= 1 && npval <= nperm && npval == round(npval)))
    stop("`npval` must be a whole number from 1 to `nperm` (", nperm, ")")
  if (crossfit && npval != 1) stop("crossfit = TRUE gives one observed statistic: `npval` must be 1")
  statistic <- if (crossfit) .cci_crossfit else cci_test_gen
  learner <- .get_learner(method)
  method_name <- .method_name(method)
  if (!is.na(seed)) set.seed(seed)
  # each name after "|" is its own V block; incomplete rows are dropped identically
  # for the observed statistic and the null draws (both use .canonical())
  observed <- function() {   # the statistic on the observed X, new random split
    statistic(formula, data, method = learner, p = p, ncomp = ncomp, scale = scale,
              max_ncomp = max_ncomp, cv_segments = cv_segments,
              method_args = method_args, reduce = reduce, reduce_ncomp = reduce_ncomp)
  }
  first <- observed()
  T0 <- cci_null_dist(formula, data, method = learner, nperm = nperm, resampling = resampling,
                      p = p, ncomp = ncomp, scale = scale, max_ncomp = max_ncomp,
                      cv_segments = cv_segments, method_args = method_args,
                      k_perm = k_perm, progress = progress,
                      reduce = reduce, reduce_ncomp = reduce_ncomp, crossfit = crossfit)
  # the other observed splits are drawn after the null draws, so that for a given
  # seed the first split and the null draws do not depend on npval
  obs <- do.call(rbind, c(list(first), lapply(seq_len(npval - 1), function(i) observed())))
  T_obs <- unname(obs[, "T"])
  # RMSEP is a loss: small T = evidence against H0 (left tail). One p-value per
  # observed split, all against the same null draws, then averaged
  p_emp_splits <- vapply(T_obs, function(t) (sum(T0 <= t) + 1) / (nperm + 1), 0)
  p_par_splits <- vapply(T_obs, .parametric_p, 0, T0 = T0)
  p_emp <- mean(p_emp_splits)
  p_par <- mean(p_par_splits)
  avg <- colMeans(obs)
  pval <- if (pvalue == "parametric") p_par else p_emp
  structure(list(statistic = c(T = avg[["T"]]), p.value = pval, M0 = avg[["M0"]], M1 = avg[["M1"]],
                 pvalue = pvalue, p.empirical = p_emp, p.parametric = p_par,
                 null.distribution = T0, method = method_name, resampling = resampling,
                 nperm = nperm, ncomp = ncomp, method_args = method_args, reduce = reduce,
                 complexity = c(M0 = avg[["a"]], M1 = avg[["b"]]), p = p,
                 npval = npval, T.splits = T_obs, p.empirical.splits = p_emp_splits,
                 p.parametric.splits = p_par_splits, crossfit = crossfit,
                 formula = formula, tail = "left"),
            class = c("BlockCCI", "htest"))
}

#' @param x a `BlockCCI` object
#' @param ... unused
#' @rdname Block.CCI.test
#' @export
print.BlockCCI <- function(x, ...) {
  red <- if (!is.null(x$reduce) && x$reduce != "none") paste("| V reduced by", toupper(x$reduce)) else ""
  cat("Block CCI test | method:", x$method, "| resampling:", x$resampling, "| R =", x$nperm, red, "\n")
  cat("Formula:", paste(deparse(x$formula), collapse = ""), "\n")
  unit <- switch(x$method,
                 SOPLS = "V components, added X components",
                 PLS   = "components of Y ~ V, Y ~ V + X",
                 ROSA  = "components of Y ~ V, Y ~ V + X",
                 XGB   = "boosting rounds of M0, M1",
                 RF    = "trees of M0, M1",
                 KNN   = "neighbours (k) of M0, M1",
                 "complexity of M0, M1")
  if (x$method %in% c("SOPLS", "PLS", "ROSA"))
    cat("ncomp =", paste(x$ncomp, collapse = ","), "| ")
  multi <- isTRUE(x$npval > 1)
  cf <- isTRUE(x$crossfit)
  cat(unit, if (cf) "(mean of the 2 folds):"
            else if (multi) paste0("(mean of ", x$npval, " observed splits):") else "(observed split):",
      paste(signif(x$complexity, 4), collapse = ", "), "\n")
  if (cf) cat("2-fold cross-fitted statistic and null draws (experimental)\n")
  cat(sprintf("RMSEP given V: %.4f   after adding X: %.4f   T = M1 - M0 = %.4f\n",
              x$M0, x$M1, x$statistic))
  if (identical(x$pvalue, "parametric")) {
    cat(sprintf("p-value = %.4f (parametric: normal null, mean %.4f, sd %.4f; empirical %.4f)\n",
                x$p.value, mean(x$null.distribution), stats::sd(x$null.distribution), x$p.empirical))
  } else {
    cat(sprintf("p-value = %.4f (share of T* <= T_obs; parametric %.4f)\n", x$p.value, x$p.parametric))
  }
  if (multi) {
    ps <- if (identical(x$pvalue, "parametric")) x$p.parametric.splits else x$p.empirical.splits
    cat(sprintf("average of %d p-values, one per observed split (range %.4f to %.4f)\n",
                x$npval, min(ps), max(ps)))
  }
  invisible(x)
}
