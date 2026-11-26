#' Null distribution of T for the block CCI test
#'
#' `nperm` Monte Carlo CV statistics ([cci_test_gen()]), each computed on a resampled
#' copy X* of X. The resampling keeps the relation between X and V (exactly or
#' approximately) but breaks any dependence between X and Y given V.
#'
#' @inheritParams cci_test_gen
#' @param nperm number of null samples R
#' @param resampling resampling of X: `"rotation"` (default; linear adjustment for
#'   V, see [rotate_X()]), `"permutation"` (global row permutation) or
#'   `"conditional"` (greedy local permutation: each row of X is replaced by the X
#'   of a row with similar V, see [permute_X_cond()])
#' @param k_perm neighbourhood size (the row itself included) for
#'   `resampling = "conditional"`, see [permute_X_cond()] (default 5)
#' @param progress show a progress bar
#' @param crossfit (experimental) if `TRUE`, every statistic is cross-fitted over two
#'   random halves of the rows instead of being computed on one train/test split: the
#'   models are fitted on each half and predict the other, the squared errors of all
#'   rows are pooled, and `T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V)` is computed once from
#'   them. Every null draw gets its own resampled X* and its own random halves. Needs
#'   two statistics' worth of model fits per draw, and `p = 0.5`.
#' @return numeric vector of `nperm` null statistics
#' @export
cci_null_dist <- function(formula, data, method = "PLS", nperm = 160,
                          resampling = c("rotation", "permutation", "conditional"),
                          p = 0.5, ncomp = "auto", scale = TRUE, max_ncomp = NULL,
                          cv_segments = 5, method_args = list(), k_perm = 5,
                          progress = FALSE, reduce = c("none", "pca", "pls"),
                          reduce_ncomp = "auto", crossfit = FALSE) {
  .check_crossfit(crossfit, p)
  statistic <- if (crossfit) .cci_crossfit else cci_test_gen
  reduce <- match.arg(reduce)
  learner <- .get_learner(method)
  resampling <- match.arg(resampling)
  if (progress) pb <- txtProgressBar(0, nperm, style = 3)
  cc <- .canonical(formula, data)
  dat <- cc$data
  T0 <- numeric(nperm)
  Vm <- if (is.null(cc$V)) NULL else do.call(cbind, cc$V)   # all of V together
  resample <- .make_resampler(cc$X, Vm, resampling, k_perm)
  for (r in seq_len(nperm)) {
    dat$X <- resample()
    T0[r] <- statistic(cc$formula, dat, method = learner, p = p, ncomp = ncomp,
                       scale = scale, max_ncomp = max_ncomp, cv_segments = cv_segments,
                       method_args = method_args, reduce = reduce,
                       reduce_ncomp = reduce_ncomp)[["T"]]
    if (progress) setTxtProgressBar(pb, r)
  }
  if (progress) close(pb)
  T0
}
