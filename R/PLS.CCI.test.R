#' PLS Computational Conditional Independence test
#'
#' Same design as [SOPLS.CCI.test()], but the regression is ordinary PLS
#' (`pls::plsr`) instead of SO-PLS. Tests Y independent of X given V with
#' `T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V)` on held-out rows, where V (all blocks
#' after the bar) and X are simply `cbind`-ed into one predictor matrix, so the
#' block structure is not used. The null distribution of T comes from rotating or
#' permuting X (R repeats of the same train/test procedure on X*). The p-value is
#' the share of null statistics at or below the observed one.
#'
#' The number of PLS components of each of the two models is chosen by cross-
#' validation on the TRAINING rows only (`ncomp = "auto"`), identically for the
#' observed and the null statistics. Under H0 the extra X columns in `Y ~ V + X`
#' can only add noise, so T is typically slightly positive; the null distribution
#' accounts for this.
#'
#' @inheritParams SOPLS.CCI.test
#' @param ncomp `"auto"` (cross-validated on the training rows) or `c(a, b)`: the
#'   number of PLS components for `Y ~ V` and for `Y ~ V + X`
#' @param max_ncomp maximum number of components searched when `ncomp = "auto"`
#' @return object of class `c("PLSCCI", "SOPLSCCI", "htest")`; [plot()] works as
#'   for [SOPLS.CCI.test()]
#' @examples
#' \dontrun{
#' data <- sim_1(n = 500, p = 4)
#' PLS.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 99)   # H0 true
#' PLS.CCI.test(Y ~ X | Z1, data, nperm = 99)       # H0 false
#' }
#' @export
PLS.CCI.test <- function(formula, data, p = 0.75, nperm = 160, ncomp = "auto",
                         method = c("rotation", "permutation", "conditional"), scale = TRUE,
                         max_ncomp = 10, cv_segments = 5,
                         pvalue = c("empirical", "parametric"),
                         group_size = 10, seed = NA, progress = TRUE) {
  method <- match.arg(method)
  pvalue <- match.arg(pvalue)
  if (!is.na(seed)) set.seed(seed)
  obs <- sopls_test_gen(formula, data, p, ncomp, scale, max_ncomp, cv_segments, learner = "pls")
  T0 <- sopls_perm_test(formula, data, nperm, method, p, ncomp, scale, max_ncomp, cv_segments,
                        progress, learner = "pls", group_size = group_size)
  # RMSEP is a loss: small T = evidence against H0 (left tail)
  p_emp <- (sum(T0 <= obs[["T"]]) + 1) / (nperm + 1)
  p_par <- .parametric_p(obs[["T"]], T0)
  pval <- if (pvalue == "parametric") p_par else p_emp
  structure(list(statistic = c(T = obs[["T"]]), p.value = pval, M0 = obs[["M0"]], M1 = obs[["M1"]],
                 pvalue = pvalue, p.empirical = p_emp, p.parametric = p_par,
                 null.distribution = T0, method = method, learner = "pls", nperm = nperm, ncomp = ncomp,
                 ncomp_selected = c(V = obs[["a"]], VX = obs[["b"]]), p = p,
                 formula = formula, tail = "left"),
            class = c("PLSCCI", "SOPLSCCI", "htest"))
}

#' @param x a `PLSCCI` object
#' @param ... unused
#' @rdname PLS.CCI.test
#' @export
print.PLSCCI <- function(x, ...) {
  cat("PLS Computational Conditional Independence test\n")
  cat("Formula:", deparse(x$formula),
      "| resampling:", x$method,
      "| R =", x$nperm,
      "| components on V, V+X (observed split) =", paste(x$ncomp_selected, collapse = ","), "\n")
  cat(sprintf("RMSEP given V: %.4f   after adding X: %.4f   T = M1 - M0 = %.4f\n", x$M0, x$M1, x$statistic))
  if (identical(x$pvalue, "parametric")) {
    cat(sprintf("p-value = %.4f (parametric: normal null, mean %.4f, sd %.4f; empirical %.4f)\n",
                x$p.value, mean(x$null.distribution), stats::sd(x$null.distribution), x$p.empirical))
  } else {
    cat(sprintf("p-value = %.4f (share of T* <= T_obs; parametric %.4f)\n", x$p.value, x$p.parametric))
  }
  invisible(x)
}
