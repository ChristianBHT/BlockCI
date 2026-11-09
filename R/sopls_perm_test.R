#' Null distribution of T for the SO-PLS CCI test
#'
#' R MCCV instances, each on a resampled X*, using `method = "rotation"`
#' (default) or `"permutation"`.
#'
#' @inheritParams sopls_test_gen
#' @param nperm number of null samples R
#' @param method `"rotation"` or `"permutation"`
#' @param progress show a progress bar
#' @return numeric vector of `nperm` null statistics
#' @export
sopls_perm_test <- function(formula, data, nperm = 160, method = c("rotation", "permutation"),
                            p = 0.75, ncomp = "auto", scale = TRUE, max_ncomp = c(10, 5),
                            cv_segments = 5, progress = FALSE) {
  method <- match.arg(method)
  if (progress) pb <- txtProgressBar(0, nperm, style = 3)
  cc <- .canonical(formula, data)
  dat <- cc$data
  T0 <- numeric(nperm)
  Vm <- if (is.null(cc$V)) NULL else do.call(cbind, cc$V)   # rotation adjusts for all of V
  for (r in seq_len(nperm)) {
    dat$X <- if (method == "rotation") rotate_X(cc$X, Vm) else permute_X(cc$X)
    T0[r] <- sopls_test_gen(cc$formula, dat, p, ncomp, scale, max_ncomp, cv_segments)[["T"]]
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

