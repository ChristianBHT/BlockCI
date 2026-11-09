
#' SO-PLS Computational Conditional Independence test
#'
#' Tests Y independent of X given V using sequential orthogonal PLS: the test
#' statistic is T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V) on held-out data; the null
#' distribution comes from rotation (or permutation) of X. The p-value is the
#' share of null statistics <= the observed one.
#'
#' @param formula `Y ~ X | V` (V optional; names refer to matrices/columns in
#'   `data`; several blocks allowed, e.g. `Y ~ X | V1 + V2`, they are cbind-ed)
#' @param data named list of matrices (or a data.frame)
#' @param p training fraction of each MCCV split
#' @param nperm number of null samples R
#' @param ncomp `"auto"` (one component per V block, number of X components chosen
#'   by inner CV) or `c(a_1, ..., a_k, b)`: components for each of the k V blocks,
#'   then components added from X
#' @param method resampling of X for the null distribution: `"rotation"` (default,
#'   linear adjustment for V), `"permutation"` (global row permutation) or
#'   `"conditional"` (permute rows of X within groups of rows with similar V, see
#'   [permute_X_cond()]; keeps nonlinear X-V dependence)
#' @param group_size rows per cluster for `method = "conditional"`
#' @param scale scale the blocks
#' @param max_ncomp,cv_segments settings for `ncomp = "auto"`
#' @param pvalue `"empirical"` (default): `(#{T* <= T_obs} + 1) / (R + 1)`.
#'   `"parametric"`: assumes the null distribution is normal with the mean and sd
#'   of the R null statistics and returns `pnorm(T_obs, mean, sd)`. This is not
#'   limited to a resolution of `1/(R + 1)`, but relies on the normal assumption.
#'   Both p-values are stored in the result as `p.empirical` and `p.parametric`.
#' @param learner regression used for M0 and M1: `"sopls"` (default) or
#'   `"xgboost"` (multi-target XGBoost, needs xgboost >= 2.0). See
#'   [sopls_test_gen()].
#' @param learner_args named list of options for the learner (for `"xgboost"`:
#'   `nrounds`, `early_stopping_rounds`, `val_fraction`, `params`)
#' @param seed optional random seed
#' @param progress show a progress bar
#' @return object of class `c("SOPLSCCI", "htest")`
#' @export
SOPLS.CCI.test <- function(formula, data, p = 0.75, nperm = 160, ncomp = "auto",
                           method = c("rotation", "permutation", "conditional"), scale = TRUE,
                           max_ncomp = c(10, 5), cv_segments = 5,
                           pvalue = c("empirical", "parametric"),
                           learner = c("sopls", "xgboost"), learner_args = list(),
                           group_size = 10, seed = NA, progress = TRUE) {
  method <- match.arg(method)
  pvalue <- match.arg(pvalue)
  learner <- match.arg(learner)
  if (!is.na(seed)) set.seed(seed)
  # each name after "|" is its own V block; incomplete rows are dropped identically
  # for the observed statistic and the null draws (both use .canonical())
  obs <- sopls_test_gen(formula, data, p, ncomp, scale, max_ncomp, cv_segments,
                        learner = learner, learner_args = learner_args)
  T0 <- sopls_perm_test(formula, data, nperm, method, p, ncomp, scale, max_ncomp, cv_segments, progress,
                        learner = learner, learner_args = learner_args, group_size = group_size)
  # RMSEP is a loss: small T = evidence against H0 (left tail)
  p_emp <- (sum(T0 <= obs[["T"]]) + 1) / (nperm + 1)
  p_par <- .parametric_p(obs[["T"]], T0)
  pval <- if (pvalue == "parametric") p_par else p_emp
  structure(list(statistic = c(T = obs[["T"]]), p.value = pval, M0 = obs[["M0"]], M1 = obs[["M1"]],
                 pvalue = pvalue, p.empirical = p_emp, p.parametric = p_par,
                 null.distribution = T0, method = method, learner = learner, nperm = nperm, ncomp = ncomp,
                 ncomp_selected = c(V = obs[["a"]], X = obs[["b"]]), p = p,
                 formula = formula, tail = "left"),
            class = c("SOPLSCCI", "htest"))
}

#' @param x a `SOPLSCCI` object
#' @param ... unused
#' @rdname SOPLS.CCI.test
#' @export
print.SOPLSCCI <- function(x, ...) {
  learner <- if (is.null(x$learner)) "sopls" else x$learner
  cat("Computational Conditional Independence test (learner:", learner, ")\n")
  if (learner == "sopls") {
    cat("Formula:", deparse(x$formula),
        "| resampling:", x$method,
        "| R =", x$nperm,
        "| ncomp =", paste(x$ncomp, collapse = ","),
        "(observed split uses V,X =", paste(x$ncomp_selected, collapse = ","), ")\n")
  } else {
    cat("Formula:", deparse(x$formula), "| resampling:", x$method, "| R =", x$nperm,
        "| rounds on V, V+X (observed split) =", paste(x$ncomp_selected, collapse = ","), "\n")
  }
  cat(sprintf("RMSEP given V: %.4f   after adding X: %.4f   T = M1 - M0 = %.4f\n", x$M0, x$M1, x$statistic))
  if (identical(x$pvalue, "parametric")) {
    cat(sprintf("p-value = %.4f (parametric: normal null, mean %.4f, sd %.4f; empirical %.4f)\n",
                x$p.value, mean(x$null.distribution), stats::sd(x$null.distribution), x$p.empirical))
  } else {
    cat(sprintf("p-value = %.4f (share of T* <= T_obs; parametric %.4f)\n", x$p.value, x$p.parametric))
  }
  invisible(x)
}

