#' Fast block CCI test (no Monte Carlo null distribution)
#'
#' Instead of building a null distribution from R refits, this compares the model
#' fitted on the real X with the same model fitted on a resampled X* (rotation by
#' default) using the SAME train/test split, and repeats this for `nsplit` splits.
#' The regression is chosen with `method`, as in [Block.CCI.test()]: `"SOPLS"`
#' (default), `"PLS"`, `"ROSA"`, `"RF"` (multi-output random forest), `"XGB"`, `"KNN"` or a
#' learner function.
#'
#' For split `k` and each test row `i`, `d_ki = l1_ki - l0_ki` is the change in
#' squared error from adding X (negative = X helps), and `d*_ki` is the same for
#' X*. The split-level statistic is `delta_k = mean_i(d_ki - d*_ki)`. A single
#' split cannot be tested on its own: the standard error from the test rows ignores
#' the variability of the training fit, which made the one-split version reject
#' true nulls far too often. Instead the `delta_k` are combined with a t statistic
#' across splits, `t = mean(delta) / se` with `K - 1` degrees of freedom (see
#' `correction`). The p-value is the left tail `pt(t, K - 1)`: small `t` means X
#' improves prediction more than X*, i.e. evidence against conditional
#' independence.
#'
#' Tuning (components, ...) is done as in [cci_test_gen()], on the training rows
#' only, separately for X and X*. M0 (Y ~ V) does not depend on X, so it is fitted
#' once per split and reused for X*; only M1 (Y ~ V + X*) is refitted. Hence
#' `d_ki - d*_ki = l1_ki - l1*_ki`: the test compares the model with X directly
#' with the same model with X* on the same rows.
#'
#' @inheritParams Block.CCI.test
#' @param method regression method: `"SOPLS"` (default), `"PLS"`, `"ROSA"`, `"XGB"`,
#'   `"RF"`, `"KNN"` (case-insensitive) or a learner function (see [cci_test_gen()]).
#' @examples
#' \dontrun{
#' data <- sim_1(n = 500, p = 4)
#' Block.CCI.fast(Y ~ X | Z1 + Z2, data)                      # SO-PLS
#' Block.CCI.fast(Y ~ X | Z1, data, method = "PLS")
#' Block.CCI.fast(Y ~ X | Z1, data, method = "ROSA", nsplit = 20)
#' Block.CCI.fast(Y ~ X | Z1 + Z2, data, method = "RF",
#'                method_args = list(num.trees = 500, num.threads = NULL))
#' }
#' @param p training fraction of each split
#' @param nsplit number of train/test splits K (at least 2)
#' @param nrot number of resampled copies X* per split, averaged (default 1)
#' @param correction `"none"` (default): plain t-test across the `nsplit` split
#'   statistics, i.e. variance `var(delta) / K`. Splits overlap, so this ignores
#'   their dependence and can be anti-conservative (0.05-0.08 rejection at the
#'   nominal 5% level in small linear simulations), but it has much more power.
#'   `"nadeau-bengio"`: corrected resampled t,
#'   `se = sqrt((1/K + n_test/n_train) * var(delta))`; valid but conservative with
#'   low power.
#' @param progress show a progress bar over the `nsplit` splits
#' @return object of class `BlockCCIfast` (a list with `statistic` (t),
#'   `p.value`, `delta` (per split), `T_obs` and `T_null` (mean over splits,
#'   RMSEP scale), `complexity` (mean over splits), ...)
#' @export
Block.CCI.fast <- function(formula, data, method = "SOPLS",
                           resampling = c("rotation", "permutation", "conditional"),
                           p = 0.5, nsplit = 10, nrot = 1, ncomp = "auto", scale = TRUE,
                           max_ncomp = NULL, cv_segments = 5, method_args = list(),
                           correction = c("none", "nadeau-bengio"),
                           k_perm = 5, seed = NA, progress = TRUE,
                           reduce = c("none", "pca", "pls"), reduce_ncomp = "auto") {
  resampling <- match.arg(resampling)
  reduce <- match.arg(reduce)
  correction <- match.arg(correction)
  learner <- .get_learner(method)
  method_name <- .method_name(method)
  if (nsplit < 2) stop("nsplit must be at least 2")
  if (!is.na(seed)) set.seed(seed)
  cc <- .canonical(formula, data)
  Vm <- if (is.null(cc$V)) NULL else do.call(cbind, cc$V)
  n <- nrow(cc$Y); n_tr <- floor(p * n); n_te <- n - n_tr
  resample <- .make_resampler(cc$X, Vm, resampling, k_perm)
  gen <- function(dat, split, yh0 = NULL)
    cci_test_gen(cc$formula, dat, method = learner, p = p, ncomp = ncomp, scale = scale,
                 max_ncomp = max_ncomp, cv_segments = cv_segments,
                 method_args = method_args, split = split, losses = TRUE, yh0 = yh0,
                 reduce = reduce, reduce_ncomp = reduce_ncomp)

  if (progress) {
    pb <- txtProgressBar(0, nsplit, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  per_split <- lapply(seq_len(nsplit), function(k) {
    split <- sample.int(n, n_tr)
    obs <- gen(cc$data, split)                 # M0 (Y ~ V) and M1 (Y ~ V + X)
    rots <- lapply(seq_len(nrot), function(i) {
      dat <- cc$data
      dat$X <- resample()
      gen(dat, split, yh0 = obs$yh0)           # M0 does not depend on X: reuse it, fit M1 only
    })
    dstar <- rowMeans(vapply(rots, function(r) r$l1 - r$l0, numeric(n_te)))
    if (progress) setTxtProgressBar(pb, k)
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
                 correction = correction, method = method_name, resampling = resampling, reduce = reduce,
                 n_train = n_tr, n_test = n_te, nsplit = nsplit, nrot = nrot,
                 T_obs = mean(vapply(per_split, `[[`, 0, "T_obs")),
                 T_null = mean(vapply(per_split, `[[`, 0, "T_null")),
                 complexity = c(M0 = mean(vapply(per_split, `[[`, 0, "a")),
                                M1 = mean(vapply(per_split, `[[`, 0, "b"))),
                 formula = formula, tail = "left"),
            class = c("BlockCCIfast", "htest"))
}

#' @param x a `BlockCCIfast` object
#' @param ... unused
#' @rdname Block.CCI.fast
#' @export
print.BlockCCIfast <- function(x, ...) {
  cat("Block CCI test, fast version (resampled t-test, correction:", x$correction, ")\n")
  cat("Formula:", paste(deparse(x$formula), collapse = ""), "| method:", x$method,
      "| resampling:", x$resampling, "\n")
  cat("Splits:", x$nsplit, "| resampled copies per split:", x$nrot,
      "| train/test rows:", x$n_train, "/", x$n_test, "\n")
  cat(sprintf("T_obs = %.4f   T_null = %.4f   (mean over splits, RMSEP scale, M1 - M0)\n",
              x$T_obs, x$T_null))
  cat(sprintf("t = %.3f (df = %d)   p-value = %.4f (left tail)\n",
              x$statistic, x$parameter, x$p.value))
  invisible(x)
}
