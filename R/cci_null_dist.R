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
#'   `"conditional"` (permute rows of X within groups of rows with similar V, see
#'   [permute_X_cond()])
#' @param group_size rows per cluster for `resampling = "conditional"`
#' @param progress show a progress bar
#' @return numeric vector of `nperm` null statistics
#' @export
cci_null_dist <- function(formula, data, method = "SOPLS", nperm = 160,
                          resampling = c("rotation", "permutation", "conditional"),
                          p = 0.75, ncomp = "auto", scale = TRUE, max_ncomp = NULL,
                          cv_segments = 5, method_args = list(), group_size = 10,
                          progress = FALSE) {
  learner <- .get_learner(method)
  resampling <- match.arg(resampling)
  if (progress) pb <- txtProgressBar(0, nperm, style = 3)
  cc <- .canonical(formula, data)
  dat <- cc$data
  T0 <- numeric(nperm)
  Vm <- if (is.null(cc$V)) NULL else do.call(cbind, cc$V)   # all of V together
  resample <- .make_resampler(cc$X, Vm, resampling, group_size)
  for (r in seq_len(nperm)) {
    dat$X <- resample()
    T0[r] <- cci_test_gen(cc$formula, dat, method = learner, p = p, ncomp = ncomp,
                          scale = scale, max_ncomp = max_ncomp, cv_segments = cv_segments,
                          method_args = method_args)[["T"]]
    if (progress) setTxtProgressBar(pb, r)
  }
  if (progress) close(pb)
  T0
}
