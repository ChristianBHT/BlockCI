# Multi-block GCM statistics (signed sum, sum of squares, max, Wald) on cross-fitted residuals from
# three regressions: OLS ("lin"), multi-target XGBoost ("xgb"), OLS then XGBoost on its residuals ("linxgb").
# Rejection rates at 0.05 on sim_1, sim_3, sim_5 (null true and false). Set `pkg` to the package folder.
# Needs: xgboost, parallel. Takes about 4 minutes with 7 workers (40 reps, n = 500).
pkg <- "C:/Users/chtho2119/Documents/GitHub/BlockCI"
lib <- file.path(tempdir(), "lib13"); dir.create(lib)
install.packages(pkg, repos = NULL, type = "source", lib = lib, quiet = TRUE)
.libPaths(c(lib, .libPaths())); Sys.setenv(R_LIBS = lib)
suppressMessages({library(CausalBlocks); library(parallel)})

# ---- multi-block GCM statistics from residual matrices ----
gcm_multi <- function(EX, EY, nsim = 20000) {
  n <- nrow(EX); p <- ncol(EX); q <- ncol(EY)
  R <- do.call(cbind, lapply(seq_len(q), function(k) EX * EY[, k]))   # n x (p*q) products
  m <- colMeans(R); S <- stats::cov(R)
  d <- sqrt(diag(S)); d[d == 0] <- 1
  Tj <- sqrt(n) * m / d
  C <- stats::cov2cor(S + diag(1e-10, ncol(S)))
  L <- chol(C + diag(1e-8, ncol(C)))
  Zs <- matrix(rnorm(nsim * ncol(C)), nsim, ncol(C)) %*% L
  c(sum_signed = 2 * stats::pnorm(-abs(sum(Tj)) / sqrt(sum(C))),
    sum_squares = (1 + sum(rowSums(Zs^2) >= sum(Tj^2))) / (nsim + 1),
    max_abs = (1 + sum(apply(abs(Zs), 1, max) >= max(abs(Tj)))) / (nsim + 1),
    wald = stats::pchisq(n * drop(crossprod(m, solve(S, m))), df = p * q, lower.tail = FALSE))
}

# ---- cross-fitted regression of M on Z: learner = "lin" | "xgb" | "linxgb" ----
crossfit_resid <- function(M, Z, folds, learner, nrounds = 300, esr = 10, val_fraction = 0.2) {
  n <- nrow(M); q <- ncol(M); res <- matrix(NA_real_, n, q)
  prm <- list(objective = "reg:squarederror", tree_method = "hist", multi_strategy = "multi_output_tree",
              max_depth = 3, eta = 0.1, nthread = 1)
  xgb_pred <- function(Ftr_idx, target, te) {      # xgboost on Z, target = matrix for rows Ftr_idx
    mu <- colMeans(target); sdv <- apply(target, 2, stats::sd); sdv[sdv == 0] <- 1
    Ts <- sweep(sweep(target, 2, mu), 2, sdv, "/")
    nv <- max(2L, ceiling(val_fraction * length(Ftr_idx)))
    vi <- sample(seq_along(Ftr_idx), nv); ti <- setdiff(seq_along(Ftr_idx), vi)
    dtr <- xgboost::xgb.DMatrix(Z[Ftr_idx[ti], , drop = FALSE], label = Ts[ti, , drop = FALSE])
    dva <- xgboost::xgb.DMatrix(Z[Ftr_idx[vi], , drop = FALSE], label = Ts[vi, , drop = FALSE])
    fit <- xgboost::xgb.train(params = prm, data = dtr, nrounds = nrounds, evals = list(val = dva),
                              early_stopping_rounds = esr, verbose = 0)
    pr <- matrix(stats::predict(fit, Z[te, , drop = FALSE]), length(te), ncol(target))
    sweep(sweep(pr, 2, sdv, "*"), 2, mu, "+")
  }
  for (k in sort(unique(folds))) {
    te <- which(folds == k); tr <- which(folds != k)
    if (learner == "xgb") {
      pred <- xgb_pred(tr, M[tr, , drop = FALSE], te)
    } else {
      D <- cbind(1, Z)
      beta <- stats::lm.fit(D[tr, , drop = FALSE], M[tr, , drop = FALSE])$coefficients
      lin_tr <- D[tr, , drop = FALSE] %*% beta; lin_te <- D[te, , drop = FALSE] %*% beta
      pred <- if (learner == "lin") lin_te else lin_te + xgb_pred(tr, M[tr, , drop = FALSE] - lin_tr, te)
    }
    res[te, ] <- M[te, , drop = FALSE] - pred
  }
  res
}

cells <- list(list(sim = "sim_1", f = Y ~ X | Z1 + Z2, H0 = TRUE),
              list(sim = "sim_1", f = Y ~ X | Z1,      H0 = FALSE),
              list(sim = "sim_3", f = Y ~ X | Z1 + Z2, H0 = TRUE),
              list(sim = "sim_3", f = Y ~ X | Z1,      H0 = TRUE),
              list(sim = "sim_3", f = Y ~ X | Z2,      H0 = FALSE),
              list(sim = "sim_5", f = Y ~ X | Z1 + Z2, H0 = TRUE),
              list(sim = "sim_5", f = Y ~ X | Z1,      H0 = FALSE))
learners <- c("lin", "xgb", "linxgb")

one_rep <- function(i, cells, learners, gcm_multi, crossfit_resid) {
  ns <- asNamespace("CausalBlocks"); dat <- list()
  do.call(rbind, lapply(seq_along(cells), function(j) {
    ce <- cells[[j]]
    if (is.null(dat[[ce$sim]])) dat[[ce$sim]] <<- get(ce$sim, ns)(500)
    cc <- ns$.canonical(ce$f, dat[[ce$sim]]); Z <- do.call(cbind, cc$V)
    folds <- sample(rep_len(1:5, nrow(Z)))
    do.call(rbind, lapply(learners, function(l) {
      ex <- crossfit_resid(cc$X, Z, folds, l); ey <- crossfit_resid(cc$Y, Z, folds, l)
      data.frame(cell = j, learner = l, t(gcm_multi(ex, ey)))
    }))
  }))
}
reps <- 40
cl <- makeCluster(max(1, detectCores() - 1)); clusterExport(cl, "lib")
clusterEvalQ(cl, { .libPaths(c(lib, .libPaths())); suppressMessages(library(CausalBlocks)) })
clusterSetRNGStream(cl, 321)
t0 <- Sys.time()
res <- do.call(rbind, parLapply(cl, seq_len(reps), one_rep, cells = cells, learners = learners,
                                gcm_multi = gcm_multi, crossfit_resid = crossfit_resid))
stopCluster(cl)
saveRDS(res, file.path(tempdir(), "gcm_linxgb_raw.rds"))
stats_cols <- c("sum_signed", "sum_squares", "max_abs", "wald")
for (l in learners) {
  cat("\n== regression:", l, "| rejection at 0.05,", reps, "reps, n = 500 ==\n")
  out <- do.call(rbind, lapply(seq_along(cells), function(j) {
    r <- res[res$cell == j & res$learner == l, stats_cols]
    data.frame(scenario = cells[[j]]$sim, conditioning = deparse(cells[[j]]$f[[3]]), H0_true = cells[[j]]$H0,
               t(round(colMeans(r < 0.05), 2)))
  }))
  print(out, row.names = FALSE)
}
cat("elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
