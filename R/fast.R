#' Fast SO-PLS conditional independence test (no Monte Carlo null distribution)
#'
#' Instead of building a null distribution from R refits, this compares SO-PLS
#' fitted on the real X with SO-PLS fitted on a rotated (or permuted) X* using the
#' SAME train/test split, and repeats this for `nsplit` splits.
#'
#' For split `k` and each test row `i`, `d_ki = l1_ki - l0_ki` is the change in
#' squared error from adding X (negative = X helps), and `d*_ki` is the same for
#' X*. The split-level statistic is `delta_k = mean_i(d_ki - d*_ki)`. A single
#' split cannot be tested on its own: the standard error from the test rows ignores
#' the variability of the training fit, which made the one-split version reject
#' true nulls far too often. Instead the `delta_k` are combined with the
#' Nadeau-Bengio corrected resampled t statistic,
#' `t = mean(delta) / sqrt((1/K + n_test/n_train) * var(delta))`, with `K - 1`
#' degrees of freedom. The p-value is the left tail `pt(t, K - 1)`: small `t` means
#' X improves prediction more than X*, i.e. evidence against conditional
#' independence.
#'
#' The number of components is chosen as in [sopls_test_gen()], on the training
#' rows only, separately for X and X*.
#'
#' @inheritParams SOPLS.CCI.test
#' @param nsplit number of train/test splits K (at least 2)
#' @param nrot number of rotated copies X* per split, averaged (default 1)
#' @param correction `"nadeau-bengio"` (default, valid but conservative with low
#'   power) or `"none"`: plain t-test across the `nsplit` split statistics, i.e.
#'   variance `var(delta) / K`. Splits overlap, so `"none"` ignores their
#'   dependence and can be anti-conservative (about 8% rejection at the nominal
#'   5% level in a small simulation), but it has much more power.
#' @return object of class `SOPLSCCIfast` (a list with `statistic` (t),
#'   `p.value`, `delta` (per split), `T_obs` and `T_null` (mean over splits,
#'   RMSEP scale), ...)
#' @export
SOPLS.CCI.fast <- function(formula, data, p = 0.5, nsplit = 10, nrot = 1, ncomp = "auto",
                           method = c("rotation", "permutation"), scale = TRUE,
                           max_ncomp = c(10, 5), cv_segments = 5,
                           correction = c("nadeau-bengio", "none"), seed = NA) {
  method <- match.arg(method)
  correction <- match.arg(correction)
  if (nsplit < 2) stop("nsplit must be at least 2")
  if (!is.na(seed)) set.seed(seed)
  f <- .parse_cci_formula(formula)
  Y <- .get_block(f$y, data); X <- .get_block(f$x, data)
  V <- if (length(f$v) == 0) NULL else lapply(f$v, .get_block, data = data)
  ok <- complete.cases(cbind(Y, X, if (!is.null(V)) do.call(cbind, V)))
  Y <- Y[ok, , drop = FALSE]
  X <- X[ok, , drop = FALSE]
  if (!is.null(V)) V <- lapply(V, function(m) m[ok, , drop = FALSE])
  Vm <- if (is.null(V)) NULL else do.call(cbind, V)
  n <- nrow(Y); n_tr <- floor(p * n); n_te <- n - n_tr

  per_split <- lapply(seq_len(nsplit), function(k) {
    split <- sample.int(n, n_tr)
    obs <- sopls_test_gen(Y, X, V, p, ncomp, scale, max_ncomp, cv_segments,
                          split = split, losses = TRUE)
    rots <- lapply(seq_len(nrot), function(i) {
      Xs <- if (method == "rotation") rotate_X(X, Vm) else permute_X(X)
      sopls_test_gen(Y, Xs, V, p, ncomp, scale, max_ncomp, cv_segments,
                     split = split, losses = TRUE)
    })
    dstar <- rowMeans(vapply(rots, function(r) r$l1 - r$l0, numeric(n_te)))
    list(delta = mean((obs$l1 - obs$l0) - dstar),
         T_obs = obs$stats[["T"]],
         T_null = mean(vapply(rots, function(r) r$stats[["T"]], 0)),
         a = obs$stats[["a"]], b = obs$stats[["b"]])
  })
  delta <- vapply(per_split, `[[`, 0, "delta")
  s2 <- stats::var(delta)
  se <- sqrt((1 / nsplit + if (correction == "nadeau-bengio") n_te / n_tr else 0) * s2)
  tstat <- if (is.finite(se) && se > 0) mean(delta) / se else NA_real_
  pval <- if (is.na(tstat)) NA_real_ else stats::pt(tstat, df = nsplit - 1)
  structure(list(statistic = c(t = tstat), parameter = c(df = nsplit - 1), p.value = pval,
                 delta = delta, mean_delta = mean(delta), se = se,
                 correction = correction, n_train = n_tr, n_test = n_te, nsplit = nsplit, nrot = nrot, method = method,
                 T_obs = mean(vapply(per_split, `[[`, 0, "T_obs")),
                 T_null = mean(vapply(per_split, `[[`, 0, "T_null")),
                 ncomp_selected = c(V = mean(vapply(per_split, `[[`, 0, "a")),
                                    X = mean(vapply(per_split, `[[`, 0, "b"))),
                 formula = formula, tail = "left"),
            class = c("SOPLSCCIfast", "htest"))
}

#' @param x a `SOPLSCCIfast` object
#' @param ... unused
#' @rdname SOPLS.CCI.fast
#' @export
print.SOPLSCCIfast <- function(x, ...) {
  cat("SO-PLS CCI test, fast version (resampled t-test, correction:", x$correction, ")\n")
  cat("Formula:", deparse(x$formula), "| resampling:", x$method, "| splits:", x$nsplit,
      "| rotations per split:", x$nrot, "| train/test rows:", x$n_train, "/", x$n_test, "\n")
  cat(sprintf("T_obs = %.4f   T_null = %.4f   (mean over splits, RMSEP scale, M1 - M0)\n",
              x$T_obs, x$T_null))
  cat(sprintf("t = %.3f (df = %d)   p-value = %.4f (left tail)\n",
              x$statistic, x$parameter, x$p.value))
  invisible(x)
}
