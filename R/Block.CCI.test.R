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
#'   (`pvalue = "parametric"`).
#'
#' The regression behind M0 (Y ~ V) and M1 (Y ~ V + X) is chosen with `method`:
#'
#' * `"SOPLS"` (default): sequential and orthogonalised PLS (`multiblock::sopls`);
#'   each V block is a separate block. With `ncomp = "auto"` each V block gets one
#'   component and the number of X components is chosen by inner CV.
#' * `"PLS"`: ordinary PLS (`pls::plsr`) on the `cbind`-ed blocks; the number of
#'   components of each model is chosen by inner CV.
#' * `"XGB"`: multi-target XGBoost (experimental; needs xgboost >= 2.0).
#' * `"RF"`: multi-output random forest (`grf::multi_regression_forest`; needs grf).
#' * a learner function with the interface described in `R/learners.R`. New methods
#'   can be made available by name by adding them to `.cci_methods` in that file.
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
#' @param seed optional random seed
#' @return object of class `c("BlockCCI", "htest")`: a list with `statistic` (T),
#'   `p.value`, `M0`, `M1`, `p.empirical`, `p.parametric`, `null.distribution`,
#'   `method`, `resampling`, `complexity` (of M0 and M1 on the observed split, see
#'   [cci_test_gen()]) and the settings used. [plot()] shows the null distribution.
#' @examples
#' \dontrun{
#' data <- sim_1(n = 500, p = 4)
#' Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99)                  # H0 true, SO-PLS
#' Block.CCI.test(Y ~ X | Z1, data, nperm = 99, method = "PLS")      # H0 false
#' Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99, method = "RF",
#'                method_args = list(num.trees = 500, num.threads = NULL))
#' Block.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 49, method = "XGB",
#'                method_args = list(nrounds = 300, params = list(max_depth = 4)))
#' }
#' @export
Block.CCI.test <- function(formula, data, method = "SOPLS",
                           resampling = c("rotation", "permutation", "conditional"),
                           nperm = 160, p = 0.75, ncomp = "auto", scale = TRUE,
                           max_ncomp = NULL, cv_segments = 5, method_args = list(),
                           pvalue = c("empirical", "parametric"),
                           group_size = 10, seed = NA, progress = TRUE) {
  resampling <- match.arg(resampling)
  pvalue <- match.arg(pvalue)
  learner <- .get_learner(method)
  method_name <- .method_name(method)
  if (!is.na(seed)) set.seed(seed)
  # each name after "|" is its own V block; incomplete rows are dropped identically
  # for the observed statistic and the null draws (both use .canonical())
  obs <- cci_test_gen(formula, data, method = learner, p = p, ncomp = ncomp, scale = scale,
                      max_ncomp = max_ncomp, cv_segments = cv_segments,
                      method_args = method_args)
  T0 <- cci_null_dist(formula, data, method = learner, nperm = nperm, resampling = resampling,
                      p = p, ncomp = ncomp, scale = scale, max_ncomp = max_ncomp,
                      cv_segments = cv_segments, method_args = method_args,
                      group_size = group_size, progress = progress)
  # RMSEP is a loss: small T = evidence against H0 (left tail)
  p_emp <- (sum(T0 <= obs[["T"]]) + 1) / (nperm + 1)
  p_par <- .parametric_p(obs[["T"]], T0)
  pval <- if (pvalue == "parametric") p_par else p_emp
  structure(list(statistic = c(T = obs[["T"]]), p.value = pval, M0 = obs[["M0"]], M1 = obs[["M1"]],
                 pvalue = pvalue, p.empirical = p_emp, p.parametric = p_par,
                 null.distribution = T0, method = method_name, resampling = resampling,
                 nperm = nperm, ncomp = ncomp, method_args = method_args,
                 complexity = c(M0 = obs[["a"]], M1 = obs[["b"]]), p = p,
                 formula = formula, tail = "left"),
            class = c("BlockCCI", "htest"))
}

#' @param x a `BlockCCI` object
#' @param ... unused
#' @rdname Block.CCI.test
#' @export
print.BlockCCI <- function(x, ...) {
  cat("Block CCI test | method:", x$method, "| resampling:", x$resampling, "| R =", x$nperm, "\n")
  cat("Formula:", paste(deparse(x$formula), collapse = ""), "\n")
  unit <- switch(x$method,
                 SOPLS = "V components, added X components",
                 PLS   = "components of Y ~ V, Y ~ V + X",
                 XGB   = "boosting rounds of M0, M1",
                 RF    = "trees of M0, M1",
                 "complexity of M0, M1")
  if (x$method %in% c("SOPLS", "PLS"))
    cat("ncomp =", paste(x$ncomp, collapse = ","), "| ")
  cat(unit, "(observed split):", paste(x$complexity, collapse = ", "), "\n")
  cat(sprintf("RMSEP given V: %.4f   after adding X: %.4f   T = M1 - M0 = %.4f\n",
              x$M0, x$M1, x$statistic))
  if (identical(x$pvalue, "parametric")) {
    cat(sprintf("p-value = %.4f (parametric: normal null, mean %.4f, sd %.4f; empirical %.4f)\n",
                x$p.value, mean(x$null.distribution), stats::sd(x$null.distribution), x$p.empirical))
  } else {
    cat(sprintf("p-value = %.4f (share of T* <= T_obs; parametric %.4f)\n", x$p.value, x$p.parametric))
  }
  invisible(x)
}
