# # This is a working script for testing out snippets and tests
# library(devtools)
# devtools::load_all("C:/Users/chtho2119/Documents/GitHub/BlockCI")
# dat <- sim_1()
# 
# pls_ci_test(Y ~ X | Z1 + Z2, data = dat, stat = "rv", resid = "cv", perm = "conditional", B = 99, group_size = 10, segments = 5, seed = 123)
# pls_ci_test(Y ~ X | Z1, data = dat, stat = "rv", resid = "cv", perm = "conditional", B = 99, group_size = 10, segments = 5, seed = 123)
# 
# rosa_ci_test(Y ~ X | Z1 + Z2, data = dat, stat = "rv", resid = "cv", perm = "conditional", B = 99, group_size = 10, segments = 5, seed = 123)
# rosa_ci_test(Y ~ X | Z1, data = dat, stat = "rv", resid = "cv", perm = "conditional", B = 99, group_size = 10, segments = 5, seed = 123)
# 
# SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = dat)
# SOPLS.CCI.test(Y ~ X | Z1, data = dat)
# SOPLS.CCI.fast(Y ~ X | Z1 + Z2, data = dat)
# SOPLS.CCI.fast(Y ~ X | Z1, data = dat)

formula = Y ~ X | Z1 + Z2
data <- sim_1()
nfolds = 5 
nrounds = 1000 
early_stopping_rounds = 10 
val_fraction = 0.2 
stat = "rv" # "cc1" "dcor"
perm = "simple"
B = 999
group_size = 10
  
stat <- match.arg(stat)
perm <- match.arg(perm)
cc <- .canonical(formula, data)            # parses formula, drops incomplete rows

Z <- do.call(cbind, cc$V)


if (!requireNamespace("xgboost", quietly = TRUE))
  stop("The package 'xgboost' (>= 2.0) is required for multixg_ci_test().")
n <- nrow(X)
if (nfolds < 2 || nfolds > n) stop("nfolds must be between 2 and the number of rows")
prm <- utils::modifyList(
  list(objective = "reg:squarederror", tree_method = "hist", multi_strategy = "multi_output_tree",
       max_depth = 3, eta = 0.1, nthread = 1),
  params)
folds <- sample(rep_len(seq_len(nfolds), n))      # shared by the X and the Y regression

crossfit <- function(M) {
  q <- ncol(M)
  res <- matrix(NA_real_, n, q)
  rounds <- numeric(nfolds)
  for (k in seq_len(nfolds)) {
    te <- which(folds == k)
    tr <- which(folds != k)
    mu <- colMeans(M[tr, , drop = FALSE])
    sdv <- apply(M[tr, , drop = FALSE], 2, stats::sd); sdv[sdv == 0] <- 1
    Ms <- sweep(sweep(M, 2, mu), 2, sdv, "/")      # targets standardised with training stats
    nv <- max(2L, ceiling(val_fraction * length(tr)))
    va <- sample(tr, nv)                           # early-stopping rows (inside the training part)
    tr2 <- setdiff(tr, va)
    dtr <- xgboost::xgb.DMatrix(Z[tr2, , drop = FALSE], label = Ms[tr2, , drop = FALSE])
    dva <- xgboost::xgb.DMatrix(Z[va, , drop = FALSE],  label = Ms[va, , drop = FALSE])
    fit <- xgboost::xgb.train(params = prm, data = dtr, nrounds = nrounds,
                              evals = list(val = dva),
                              early_stopping_rounds = early_stopping_rounds, verbose = 0)
    best <- xgboost::xgb.attr(fit, "best_iteration")
    rounds[k] <- if (is.null(best)) xgboost::xgb.get.num.boosted.rounds(fit) else as.integer(best) + 1L
    pr <- matrix(stats::predict(fit, Z[te, , drop = FALSE]), length(te), q)
    res[te, ] <- M[te, , drop = FALSE] - sweep(sweep(pr, 2, sdv, "*"), 2, mu, "+")
  }
  list(res = res, rounds = mean(rounds))


rs <- .multixg_residuals(cc$X, cc$Y, Z, nfolds, nrounds, early_stopping_rounds                           val_fraction, params)


res <- .resid_perm_test(rs$res_x, rs$res_y, Z, stat, perm, B, group_size)
  
  