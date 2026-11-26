# Greedy local (k nearest neighbour) permutation -----------------------------------
# Runge, J. (2018). Conditional independence testing based on a nearest-neighbor
# estimator of conditional mutual information. AISTATS (arXiv:1709.01447),
# Algorithm 1. Used by permute_X_cond(), resampling = "conditional" in the CCI tests
# and perm = "conditional" in the residual test.

.check_k_perm <- function(k_perm, n) {
  if (!is.numeric(k_perm) || length(k_perm) != 1 ||
      !isTRUE(k_perm >= 2 && k_perm <= n && k_perm == round(k_perm)))
    stop("`k_perm` must be a whole number from 2 to the number of rows (", n, ")")
  as.integer(k_perm)
}

# n x k_perm matrix: row i holds the indices of the k_perm nearest neighbours of row i
# of V, row i itself included. As in the paper, every column of V is replaced by its
# ranks (ties broken at random) and distances use the maximum norm. If several rows
# are tied at the distance of the k_perm-th neighbour, the remaining places are
# filled at random among them.
# Exact search: the ranks of the first column are 1, ..., n, and a row whose first
# rank differs by more than w is farther away than w in the maximum norm. So only a
# window of rows around row i has to be searched; it is widened until it is certain
# to contain all k_perm neighbours. This saves much time when V has few columns; with
# many columns the window covers nearly all rows (time grows as n^2 * ncol(V)).
.knn_index <- function(V, k_perm = 5) {
  V <- as.matrix(V)
  n <- nrow(V); d <- ncol(V)
  k <- .check_k_perm(k_perm, n)
  R <- matrix(vapply(seq_len(d), function(j) rank(V[, j], ties.method = "random"),
                     integer(n)), n, d)
  first <- R[, 1]
  by_rank <- order(first)                               # by_rank[r]: the row with rank r
  w0 <- min(n, max(k, ceiling(n * (k / n)^(1 / d))))    # initial half-width of the window
  nb <- matrix(0L, n, k)
  for (i in seq_len(n)) {
    r <- first[i]; w <- w0
    repeat {
      lo <- max(1, r - w); hi <- min(n, r + w)
      cand <- by_rank[lo:hi]
      dst <- abs(first[cand] - r)
      for (j in seq_len(d)[-1]) dst <- pmax.int(dst, abs(R[cand, j] - R[i, j]))
      dst[cand == i] <- -1L                             # row i is always its own neighbour
      full <- lo == 1 && hi == n
      if (length(cand) >= k) {
        dk <- sort.int(dst, partial = k)[k]
        if (full || dk <= w) break                      # no row outside can be as close
        w <- dk
      } else w <- 2 * w
    }
    inside <- which(dst < dk)
    tied <- which(dst == dk)
    nb[i, ] <- cand[c(inside, tied[sample.int(length(tied), k - length(inside))])]
  }
  nb
}

# One draw of Algorithm 1 from the neighbour matrix `nb`: every neighbour list is
# shuffled, the rows are visited in random order, and row i takes the first neighbour
# in its list that no earlier row has taken; if all of them are taken, it takes the
# last one in the list. Returns n row indices. Because of that last rule an index can
# occur more than once, so this is not always a strict permutation.
.local_perm <- function(nb) {
  n <- nrow(nb); k <- ncol(nb)
  shuffled <- matrix(nb[order(rep.int(seq_len(n), k), stats::runif(n * k))], k, n)  # column i: row i
  used <- logical(n)
  idx <- integer(n)
  for (i in sample.int(n)) {
    s <- shuffled[, i]
    free <- s[!used[s]]
    j <- if (length(free)) free[1L] else s[k]
    idx[i] <- j
    used[j] <- TRUE
  }
  idx
}

#' Conditional (local nearest-neighbour) permutation of X given V
#'
#' The greedy local permutation scheme of Runge (2018, Algorithm 1). For every row
#' the `k_perm` nearest neighbours in V are found (the row itself included; each
#' column of V rank-transformed, maximum norm). The rows are then visited in random
#' order, and each receives the X row of a randomly chosen neighbour that has not
#' been used yet; if all its neighbours are used, one of them is used again. Every
#' row of X* thus comes from a row with a similar V, which keeps the (possibly
#' nonlinear) dependence between X and V approximately intact while breaking the
#' dependence between X and Y given V.
#'
#' Because a neighbour can be used again, some rows of X may occur more than once in
#' the result and others not at all. A row can also keep its own X.
#'
#' Neighbours are local only when V has few columns relative to the number of rows;
#' with many columns even the nearest rows can be far apart.
#'
#' @param X numeric matrix to resample
#' @param V numeric matrix of conditioning variables, or `NULL` (then a global
#'   permutation is used)
#' @param k_perm size of the neighbourhood, the row itself included (default 5;
#'   Runge recommends about 5 to 10). A smaller value stays closer to the dependence
#'   between X and V but changes X less, which costs power.
#' @return matrix of the same dimension as `X`, made of rows of `X`
#' @references Runge, J. (2018). Conditional independence testing based on a
#'   nearest-neighbor estimator of conditional mutual information. *Proceedings of
#'   AISTATS*. \url{https://arxiv.org/abs/1709.01447}
#' @export
permute_X_cond <- function(X, V = NULL, k_perm = 5) {
  if (is.null(V)) return(permute_X(X))
  X[.local_perm(.knn_index(V, k_perm)), , drop = FALSE]
}
