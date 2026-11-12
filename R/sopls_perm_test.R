#' Null distribution of T for the SO-PLS CCI test
#'
#' R MCCV instances, each on a resampled X*, using `method = "rotation"`
#' (default), `"permutation"` or `"conditional"`.
#'
#' @inheritParams sopls_test_gen
#' @param nperm number of null samples R
#' @param method resampling of X: `"rotation"` (linear adjustment for V),
#'   `"permutation"` (global row permutation) or `"conditional"` (permute rows of X
#'   within groups of rows with similar V, see [permute_X_cond()]; suited to
#'   nonlinear dependence on V)
#' @param progress show a progress bar
#' @param group_size rows per cluster for `method = "conditional"`
#' @return numeric vector of `nperm` null statistics
#' @export
sopls_perm_test <- function(formula, data, nperm = 160,
                            method = c("rotation", "permutation", "conditional"),
                            p = 0.75, ncomp = "auto", scale = TRUE, max_ncomp = c(10, 5),
                            cv_segments = 5, progress = FALSE,
                            learner = c("sopls", "xgboost", "pls"), learner_args = list(),
                            group_size = 10) {
  learner <- match.arg(learner)
  method <- match.arg(method)
  if (progress) pb <- txtProgressBar(0, nperm, style = 3)
  cc <- .canonical(formula, data)
  dat <- cc$data
  T0 <- numeric(nperm)
  Vm <- if (is.null(cc$V)) NULL else do.call(cbind, cc$V)   # all of V together
  resample <- .make_resampler(cc$X, Vm, method, group_size)
  for (r in seq_len(nperm)) {
    dat$X <- resample()
    T0[r] <- sopls_test_gen(cc$formula, dat, p, ncomp, scale, max_ncomp, cv_segments,
                            learner = learner, learner_args = learner_args)[["T"]]
    if (progress) setTxtProgressBar(pb, r)
  }
  if (progress) close(pb)
  T0
}

# ---- user interface (CCI.test analogue) ------------------------------------

# Left-tail p-value assuming a normal null with the sample mean and sd of T0
.parametric_p <- function(Tobs, T0) {
  s <- stats::sd(T0)
  if (!is.finite(s) || s == 0) return(NA_real_)
  stats::pnorm(Tobs, mean = mean(T0), sd = s)
}

