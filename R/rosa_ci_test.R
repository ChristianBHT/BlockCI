#' Test of conditional independence Y _||_ X | Z via ROSA residuals and permutation
#'
#' Removes the conditioning blocks Z from X and from Y with ROSA (response-oriented
#' sequential alternation, `multiblock::rosa`), keeping the Z blocks separate,
#' and then tests whether the two residual matrices are still dependent
#' (RV coefficient, first canonical correlation or distance correlation) with a
#' permutation test. Same idea as [pls_ci_test()], with ROSA instead of PLS.
#'
#' @param formula `Y ~ X | Z` (the conditioning set is required). Several blocks
#'   may be given after the bar, e.g. `Y ~ X | Z1 + Z2`; each is a separate ROSA block.
#' @param data named list of matrices (or a data.frame) holding the variables in
#'   `formula`. Rows with missing values are dropped.
#' @param ncomp number of ROSA components used to remove Z, or `"auto"`
#'   (chosen by cross-validation, separately for X and Y).
#' @param max_ncomp upper limit for `ncomp` when `ncomp = "auto"`.
#' @param stat `"rv"`, `"cc1"` or `"dcor"` (the last needs the `energy` package).
#' @param resid `"cv"` (cross-validated residuals, recommended) or `"fitted"`
#'   (in-sample residuals).
#' @param perm `"conditional"` (permute residual rows within groups of
#'   observations with similar Z) or `"simple"` (permute all rows).
#' @param B number of permutations.
#' @param group_size desired group size for conditional permutation.
#' @param segments number of CV segments.
#' @param seed optional seed.
#' @return an object of class `"htest"`, with the chosen numbers of components in
#'   `ncomp` and the permutation statistics in `perm_stats`.
#' @export
rosa_ci_test <- function(formula, data,
                         ncomp = "auto", max_ncomp = 10,
                         stat = c("rv", "cc1", "dcor"),
                         resid = c("cv", "fitted"),
                         perm = c("conditional", "simple"),
                         B = 999, group_size = 10,
                         segments = 10, seed = NULL) {
  stat  <- match.arg(stat)
  resid <- match.arg(resid)
  perm  <- match.arg(perm)
  if (!is.null(seed)) set.seed(seed)

  cc <- .canonical(formula, data)            # parses formula, drops incomplete rows
  if (is.null(cc$V)) stop("formula must have a conditioning set: Y ~ X | Z")
  X <- cc$X; Y <- cc$Y
  Zl <- cc$V                                 # list of Z blocks, kept separate for ROSA
  Z <- do.call(cbind, Zl)
  n <- nrow(X)
  vn <- paste0("V", seq_along(Zl))
  rosa_formula <- stats::as.formula(paste("M ~", paste(vn, collapse = " + ")))

  # regress a block M on all Z blocks with ROSA and return its residuals
  remove_Z <- function(M) {
    dat <- c(list(M = M), stats::setNames(Zl, vn))
    k_max <- min(if (identical(ncomp, "auto")) max_ncomp else ncomp, ncol(Z), n - 2)
    fit <- multiblock::rosa(rosa_formula, data = dat, ncomp = k_max, scale = TRUE,
                            validation = "CV", segments = segments)
    cvp <- fit$validation$pred                                   # n x q x ncomp
    nc <- if (identical(ncomp, "auto")) {
      which.min(vapply(seq_len(dim(cvp)[3]),
                       function(a) sum(colMeans((M - matrix(cvp[, , a], n))^2)), 0))
    } else k_max
    pred <- if (resid == "cv") cvp[, , nc] else stats::fitted(fit)[, , nc]
    list(res = M - matrix(pred, nrow = n), ncomp = nc)
  }

  rx <- remove_Z(X)
  ry <- remove_Z(Y)
  res <- .resid_perm_test(rx$res, ry$res, Z, stat, perm, B, group_size)

  structure(list(
    statistic   = c(stat = res$obs),
    p.value     = res$p,
    method      = sprintf("ROSA residual test for X _||_ Y | Z (%s, %s permutation, %s residuals)",
                          stat, perm, resid),
    data.name   = paste(deparse(formula), collapse = ""),
    null.value  = NULL,
    alternative = "X and Y are dependent given Z",
    ncomp       = c(X = rx$ncomp, Y = ry$ncomp),
    perm_stats  = res$perms
  ), class = "htest")
}
