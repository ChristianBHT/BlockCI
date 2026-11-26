# Multi-output k nearest neighbour regression ----------------------------------------
# Shared by the "KNN" learner of the CCI tests (learners.R) and of the residual test
# (resid_learners.R). The prediction for a row is the mean response (all columns) of
# its k nearest training rows. Distances are Euclidean, with every feature
# standardised by its training mean and sd (scale = TRUE). Neighbours are found with
# FNN::get.knnx (kd-tree). One k is used for all response columns together.

.knn_default_grid <- c(3, 5, 7, 10, 15, 20, 30, 50, 75, 100)

# options from method_args: k ("auto" or a number of neighbours) and k_grid (the
# candidates for k = "auto")
.knn_args <- function(args) {
  args <- args %||% list()
  unknown <- setdiff(names(args), c("k", "k_grid"))
  if (length(args) && (is.null(names(args)) || length(unknown)))
    stop("for KNN, method_args may only contain `k` and `k_grid`")
  k <- args[["k"]] %||% "auto"       # [[ ]]: args$k would also match k_grid
  grid <- args[["k_grid"]] %||% .knn_default_grid
  whole <- function(x) is.numeric(x) && length(x) > 0 && !anyNA(x) && all(x >= 1 & x == round(x))
  if (!identical(k, "auto") && !(length(k) == 1 && whole(k)))
    stop("for KNN, `k` must be \"auto\" or a whole number >= 1")
  if (!whole(grid)) stop("for KNN, `k_grid` must hold whole numbers >= 1")
  list(k = k, grid = sort(unique(grid)))
}

# training and new feature matrices on a common scale
.knn_features <- function(Fm, train, test, scale) {
  Fm <- as.matrix(Fm)
  if (scale) {   # training means and sds; columns constant on the training rows are dropped
    Fs <- .std_cols(Fm, train)
    if (ncol(Fs) > 0) Fm <- Fs
  }
  list(tr = Fm[train, , drop = FALSE], te = Fm[test, , drop = FALSE])
}

# indices (nrow(Fte) x k) of the k nearest rows of Ftr for every row of Fte
.knn_neighbours <- function(Ftr, Fte, k) {
  if (!requireNamespace("FNN", quietly = TRUE))
    stop("The package 'FNN' is required for method = \"KNN\".")
  FNN::get.knnx(Ftr, Fte, k = k, algorithm = "kd_tree")$nn.index
}

# For each k in `ks` (increasing) the predictions from the first k neighbours in `nb`:
# a list of nrow(nb) x ncol(Ytr) matrices. The sums are built up neighbour by neighbour.
.knn_means <- function(Ytr, nb, ks) {
  total <- matrix(0, nrow(nb), ncol(Ytr))
  out <- vector("list", length(ks))
  done <- 0
  for (a in seq_along(ks)) {
    for (j in (done + 1):ks[a]) total <- total + Ytr[nb[, j], , drop = FALSE]
    done <- ks[a]
    out[[a]] <- total / done
  }
  out
}

# k from `grid` with the smallest cross-validated squared error (summed over the
# response columns) on the training data Ftr, Ytr.
.knn_choose_k <- function(Ftr, Ytr, grid, segments = 5) {
  n <- nrow(Ftr)
  segments <- max(2, min(segments, n))
  fold <- sample(rep_len(seq_len(segments), n))
  smallest_fit <- n - max(tabulate(fold))               # fewest rows to search among
  ks <- grid[grid <= smallest_fit]
  if (!length(ks)) ks <- smallest_fit
  if (length(ks) == 1) return(as.integer(ks))
  sse <- numeric(length(ks))
  for (s in seq_len(segments)) {
    out <- fold == s
    nb <- .knn_neighbours(Ftr[!out, , drop = FALSE], Ftr[out, , drop = FALSE], max(ks))
    pred <- .knn_means(Ytr[!out, , drop = FALSE], nb, ks)
    sse <- sse + vapply(pred, function(P) sum((Ytr[out, , drop = FALSE] - P)^2), 0)
  }
  as.integer(ks[which.min(sse)])
}

# Fit on the rows `train` of (Fm, Y) and predict the rows `test`. k = "auto" is chosen
# by cross-validation within the training rows. Returns list(pred, k).
.knn_fit_predict <- function(Fm, Y, train, test, k = "auto", grid = .knn_default_grid,
                             scale = TRUE, segments = 5) {
  Fs <- .knn_features(Fm, train, test, scale)
  Ytr <- as.matrix(Y)[train, , drop = FALSE]
  k <- if (identical(k, "auto")) .knn_choose_k(Fs$tr, Ytr, grid, segments)
       else as.integer(min(k, length(train)))
  nb <- .knn_neighbours(Fs$tr, Fs$te, k)
  list(pred = .knn_means(Ytr, nb, k)[[1]], k = k)
}
