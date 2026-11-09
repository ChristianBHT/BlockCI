
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
#' @param method `"rotation"` (default) or `"permutation"` resampling of X
#' @param scale scale the blocks
#' @param max_ncomp,cv_segments settings for `ncomp = "auto"`
#' @param pvalue `"empirical"` (default): `(#{T* <= T_obs} + 1) / (R + 1)`.
#'   `"parametric"`: assumes the null distribution is normal with the mean and sd
#'   of the R null statistics and returns `pnorm(T_obs, mean, sd)`. This is not
#'   limited to a resolution of `1/(R + 1)`, but relies on the normal assumption.
#'   Both p-values are stored in the result as `p.empirical` and `p.parametric`.
#' @param seed optional random seed
#' @param progress show a progress bar
#' @return object of class `c("SOPLSCCI", "htest")`
#' @export
SOPLS.CCI.test <- function(formula, data, p = 0.75, nperm = 160, ncomp = "auto",
                           method = c("rotation", "permutation"), scale = TRUE,
                           max_ncomp = c(10, 5), cv_segments = 5,
                           pvalue = c("empirical", "parametric"),
                           seed = NA, progress = TRUE) {
  method <- match.arg(method)
  pvalue <- match.arg(pvalue)
  if (!is.na(seed)) set.seed(seed)
  f <- .parse_cci_formula(formula)
  Y <- .get_block(f$y, data); X <- .get_block(f$x, data)
  # each name after "|" is its own V block (one latent variable per block)
  V <- if (length(f$v) == 0) NULL else lapply(f$v, .get_block, data = data)
  ok <- complete.cases(cbind(Y, X, if (!is.null(V)) do.call(cbind, V)))
  Y <- Y[ok, , drop = FALSE]
  X <- X[ok, , drop = FALSE]
  if (!is.null(V)) V <- lapply(V, function(m) m[ok, , drop = FALSE])
  obs <- sopls_test_gen(Y, X, V, p, ncomp, scale, max_ncomp, cv_segments)
  T0 <- sopls_perm_test(Y, X, V, nperm, method, p, ncomp, scale, max_ncomp, cv_segments, progress)
  # RMSEP is a loss: small T = evidence against H0 (left tail)
  p_emp <- (sum(T0 <= obs[["T"]]) + 1) / (nperm + 1)
  p_par <- .parametric_p(obs[["T"]], T0)
  pval <- if (pvalue == "parametric") p_par else p_emp
  structure(list(statistic = c(T = obs[["T"]]), p.value = pval, M0 = obs[["M0"]], M1 = obs[["M1"]],
                 pvalue = pvalue, p.empirical = p_emp, p.parametric = p_par,
                 null.distribution = T0, method = method, nperm = nperm, ncomp = ncomp,
                 ncomp_selected = c(V = obs[["a"]], X = obs[["b"]]), p = p,
                 formula = formula, tail = "left"),
            class = c("SOPLSCCI", "htest"))
}

#' @param x a `SOPLSCCI` object
#' @param ... unused
#' @rdname SOPLS.CCI.test
#' @export
print.SOPLSCCI <- function(x, ...) {
  cat("SO-PLS Computational Conditional Independence test\n")
  cat("Formula:", deparse(x$formula), 
      "| resampling:", x$method, 
      "| R =", x$nperm,
      "| ncomp =", paste(x$ncomp, collapse = ","),
      "(observed split uses V,X =", paste(x$ncomp_selected, collapse = ","), ")\n")
  cat(sprintf("RMSEP given V: %.4f   after adding X: %.4f   T = M1 - M0 = %.4f\n", x$M0, x$M1, x$statistic))
  if (identical(x$pvalue, "parametric")) {
    cat(sprintf("p-value = %.4f (parametric: normal null, mean %.4f, sd %.4f; empirical %.4f)\n",
                x$p.value, mean(x$null.distribution), stats::sd(x$null.distribution), x$p.empirical))
  } else {
    cat(sprintf("p-value = %.4f (share of T* <= T_obs; parametric %.4f)\n", x$p.value, x$p.parametric))
  }
  invisible(x)
}

