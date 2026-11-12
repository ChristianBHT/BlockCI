#' Multi-target XGBoost Computational Conditional Independence test
#'
#' Same design as [SOPLS.CCI.test()], but the regression is multi-target XGBoost
#' (`learner = "xgboost"`, see `.learner_xgboost()` in `R/learners.R`; needs
#' xgboost >= 2.0 with `multi_strategy = "multi_output_tree"`). Tests Y independent
#' of X given V with `T = RMSEP(Y ~ V + X) - RMSEP(Y ~ V)` on held-out rows; V (all
#' blocks after the bar) and X are `cbind`-ed. The null distribution of T comes from
#' resampling X (rotation, permutation or conditional) `nperm` times. The p-value is
#' the share of null statistics at or below the observed one.
#'
#' The number of boosting rounds of each model is chosen by early stopping on a
#' validation part of the TRAINING rows only.
#'
#' **Experimental, needs more work:** in our simulations the XGBoost-based tests
#' reject true nulls too often on nonlinear data, because rotation / permutation of
#' X does not condition on V nonlinearly (see `LOG.md`). It prints a message when used.
#'
#' @inheritParams SOPLS.CCI.test
#' @param nrounds maximum number of boosting rounds
#' @param early_stopping_rounds stop after this many rounds without validation improvement
#' @param val_fraction share of the training rows used as inner validation set
#' @param params named list merged into the xgboost parameters (defaults
#'   `max_depth = 3, eta = 0.1, tree_method = "hist",
#'   multi_strategy = "multi_output_tree", nthread = 1`)
#' @return object of class `c("XGBCCI", "SOPLSCCI", "htest")`; `print()` and
#'   [plot()] work as for [SOPLS.CCI.test()]
#' @examples
#' \dontrun{
#' data <- sim_1(n = 500, p = 4)
#' multixgb.CCI.test(Y ~ X | Z1 + Z2, data, nperm = 49)   # H0 true
#' multixgb.CCI.test(Y ~ X | Z1, data, nperm = 49)        # H0 false
#' }
#' @export
multixgb.CCI.test <- function(formula, data, p = 0.75, nperm = 160,
                              method = c("rotation", "permutation", "conditional"),
                              nrounds = 300, early_stopping_rounds = 10, val_fraction = 0.2,
                              params = list(),
                              pvalue = c("empirical", "parametric"),
                              group_size = 10, seed = NA, progress = TRUE) {
  method <- match.arg(method)
  pvalue <- match.arg(pvalue)
  if (!is.na(seed)) set.seed(seed)
  la <- list(nrounds = nrounds, early_stopping_rounds = early_stopping_rounds,
             val_fraction = val_fraction, params = params)
  obs <- sopls_test_gen(formula, data, p, learner = "xgboost", learner_args = la)
  T0 <- sopls_perm_test(formula, data, nperm, method, p, progress = progress,
                        learner = "xgboost", learner_args = la, group_size = group_size)
  # RMSEP is a loss: small T = evidence against H0 (left tail)
  p_emp <- (sum(T0 <= obs[["T"]]) + 1) / (nperm + 1)
  p_par <- .parametric_p(obs[["T"]], T0)
  pval <- if (pvalue == "parametric") p_par else p_emp
  structure(list(statistic = c(T = obs[["T"]]), p.value = pval, M0 = obs[["M0"]], M1 = obs[["M1"]],
                 pvalue = pvalue, p.empirical = p_emp, p.parametric = p_par,
                 null.distribution = T0, method = method, learner = "xgboost", nperm = nperm,
                 ncomp = NA, ncomp_selected = c(V = obs[["a"]], VX = obs[["b"]]), p = p,
                 formula = formula, tail = "left"),
            class = c("XGBCCI", "SOPLSCCI", "htest"))
}
