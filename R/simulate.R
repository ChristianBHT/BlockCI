# Simulators with known ground truth ------------------------------------------

# Laplace noise
#' @noRd
rlap <- function(n) { u <- stats::runif(n) - 0.5; -sign(u) * log(1 - 2 * abs(u)) }

#' Simulate blocks measured on the same samples (multi-view)
#'
#' Setting A: same samples, different variables per block. Each block k is
#' driven by a latent factor z_k; a block-level DAG acts on the factors.
#'
#' @param n number of samples
#' @param p integer vector, number of variables in each block
#' @param A K x K block adjacency; `A[i, j] = 1` means block i -> block j
#'   (must be topologically ordered)
#' @param noise_sd standard deviation of the measurement noise
#' @return list with `blocks` (list of matrices), `truth` (adjacency), `setting`
#' @export
sim_multiview <- function(n = 500, p = c(6, 6, 6),
                          A = rbind(c(0,1,1), c(0,0,1), c(0,0,0)), noise_sd = 0.5) {
  K <- length(p); z <- matrix(0, n, K)
  for (k in seq_len(K)) {
    pa <- which(A[, k] == 1)
    z[, k] <- rlap(n)
    for (j in pa) z[, k] <- z[, k] + stats::runif(1, 0.7, 1.2) * sample(c(-1, 1), 1) * z[, j]
  }
  blocks <- lapply(seq_len(K), function(k) {
    load <- matrix(stats::runif(p[k], 0.6, 1.2) * sample(c(-1, 1), p[k], TRUE), 1)
    z[, k, drop = FALSE] %*% load + matrix(stats::rnorm(n * p[k], sd = noise_sd), n, p[k])
  })
  names(blocks) <- paste0("B", seq_len(K))
  dimnames(A) <- list(names(blocks), names(blocks))
  list(blocks = blocks, truth = A, setting = "multiview")
}

#' Simulate blocks as environments (same variables, different samples)
#'
#' Setting B: linear SEM over p variables; environments 2..E apply
#' interventions (noise shifts) to the `intervened` variables.
#'
#' @param n_env number of environments (blocks)
#' @param n samples per environment
#' @param p number of variables
#' @param edge_prob probability of an edge between two variables
#' @param intervened indices of intervened variables (random if `NULL`)
#' @return list with `blocks`, `truth` (`truth[i, j] = 1` means i -> j),
#'   `intervened`, `setting`
#' @export
sim_environments <- function(n_env = 4, n = 300, p = 5, edge_prob = 0.4, intervened = NULL) {
  W <- matrix(0, p, p)
  for (i in seq_len(p - 1)) for (j in (i + 1):p)
    if (stats::runif(1) < edge_prob) W[i, j] <- stats::runif(1, 0.6, 1.2) * sample(c(-1, 1), 1)
  if (is.null(intervened)) intervened <- sample(p, max(1, p %/% 2))
  shift <- stats::runif(p, 1.5, 2.5); off <- stats::runif(p, 1, 2) * sample(c(-1, 1), p, TRUE)
  blocks <- lapply(seq_len(n_env), function(e) {
    X <- matrix(0, n, p)
    for (j in seq_len(p)) {
      eps <- stats::rnorm(n)
      if (e > 1 && j %in% intervened) eps <- eps * shift[j] + off[j] * (e - 1)
      X[, j] <- eps + if (j > 1) X[, seq_len(j - 1), drop = FALSE] %*% W[seq_len(j - 1), j] else 0
    }
    colnames(X) <- paste0("V", seq_len(p)); X
  })
  names(blocks) <- paste0("E", seq_len(n_env))
  dimnames(W) <- list(colnames(blocks[[1]]), colnames(blocks[[1]]))
  list(blocks = blocks, truth = (W != 0) * 1, intervened = intervened, setting = "environments")
}
