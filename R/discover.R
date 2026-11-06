# Causal discovery between blocks ---------------------------------------------

# ---- Setting A: same samples, different variables -----------------------------

#' Summarise each block by its first principal component score
#' @param blocks named list of matrices with the same number of rows
#' @return matrix with one column (score) per block
#' @export
block_scores <- function(blocks) {
  S <- sapply(blocks, function(X) stats::prcomp(scale(X), rank. = 1)$x[, 1])
  colnames(S) <- names(blocks); S
}

# Fisher-z partial correlation test; returns p-value
#' @noRd
pcor_test <- function(C, n, i, j, S) {
  idx <- c(i, j, S); P <- tryCatch(solve(C[idx, idx, drop = FALSE]), error = function(e) NULL)
  if (is.null(P)) return(0)
  r <- -P[1, 2] / sqrt(P[1, 1] * P[2, 2]); r <- max(min(r, 0.9999), -0.9999)
  2 * stats::pnorm(-abs(sqrt(n - length(S) - 3) * atanh(r)))
}

# Hyvarinen-Smith pairwise likelihood ratio (assumes super-Gaussian noise, e.g. Laplace).
# Returns TRUE if x -> y is the more plausible direction.
#' @noRd
orient_pair <- function(x, y) {
  x <- as.numeric(scale(x)); y <- as.numeric(scale(y)); rho <- stats::cor(x, y)
  rho * mean(x * tanh(y) - tanh(x) * y) > 0
}

#' PC-style causal discovery between blocks measured on the same samples
#'
#' Summarises each block by its first principal component, runs a PC-style
#' skeleton search (Fisher-z partial correlation tests), orients v-structures,
#' and orients remaining edges with a pairwise non-Gaussian likelihood ratio.
#'
#' @param blocks named list of matrices (same rows, different variables)
#' @param alpha significance level of the conditional independence tests
#' @return directed adjacency matrix between blocks (`[i, j] = 1`: i -> j)
#' @export
discover_multiview <- function(blocks, alpha = 0.01) {
  S <- block_scores(blocks); n <- nrow(S); K <- ncol(S); C <- stats::cor(S)
  adj <- matrix(1, K, K); diag(adj) <- 0; sep <- list()
  for (size in 0:(K - 2)) for (i in 1:(K - 1)) for (j in (i + 1):K) {
    if (!adj[i, j]) next
    nb <- setdiff(which(adj[i, ] == 1), j)
    if (length(nb) < size) next
    subs <- if (size == 0) list(integer(0)) else utils::combn(nb, size, simplify = FALSE)
    for (s in subs) if (pcor_test(C, n, i, j, s) > alpha) {
      adj[i, j] <- adj[j, i] <- 0; sep[[paste(i, j)]] <- s; break
    }
  }
  D <- adj  # undirected edge = both directions set
  for (k in 1:K) {
    nb <- which(adj[k, ] == 1)
    for (i in nb) for (j in nb) if (i < j && !adj[i, j] && !(k %in% sep[[paste(i, j)]])) {
      D[k, i] <- 0; D[k, j] <- 0   # v-structure i -> k <- j
    }
  }
  for (i in 1:(K - 1)) for (j in (i + 1):K) if (D[i, j] == 1 && D[j, i] == 1) {
    if (orient_pair(S[, i], S[, j])) D[j, i] <- 0 else D[i, j] <- 0
  }
  dimnames(D) <- list(colnames(S), colnames(S)); D
}

# ---- Setting B: same variables, different samples (environments) ---------------

#' Invariant Causal Prediction (simplified) across environments
#'
#' Blocks are environments sharing columns. Predictor sets whose residuals look
#' invariant across blocks (mean and variance) are accepted; the estimated
#' parents of a variable are the intersection of its accepted sets.
#'
#' @param blocks list of matrices with identical columns (one per environment)
#' @param alpha significance level
#' @return adjacency matrix `W` with `W[i, j] = 1` for i -> j
#' @export
discover_environments <- function(blocks, alpha = 0.05) {
  p <- ncol(blocks[[1]]); env <- factor(rep(seq_along(blocks), sapply(blocks, nrow)))
  X <- do.call(rbind, blocks); W <- matrix(0, p, p, dimnames = list(colnames(X), colnames(X)))
  for (j in seq_len(p)) {
    others <- setdiff(seq_len(p), j); accepted <- list()
    for (size in 0:length(others)) {
      subs <- if (size == 0) list(integer(0)) else utils::combn(others, size, simplify = FALSE)
      for (s in subs) {
        r <- if (length(s)) stats::resid(stats::lm(X[, j] ~ X[, s, drop = FALSE])) else X[, j] - mean(X[, j])
        p_mean <- stats::anova(stats::lm(r ~ env))[["Pr(>F)"]][1]
        p_var  <- stats::anova(stats::lm(abs(r - stats::ave(r, env, FUN = stats::median)) ~ env))[["Pr(>F)"]][1]
        if (min(p_mean, p_var) * 2 > alpha) accepted[[length(accepted) + 1]] <- s
      }
    }
    if (length(accepted)) { pa <- Reduce(intersect, accepted); W[pa, j] <- 1 }
  }
  W
}
