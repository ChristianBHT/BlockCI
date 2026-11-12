#' Simulate blocks Y, X, Z1, Z2 for large-N and high-dimensional studies
#'
#' Memory-friendly generator with a known null. The conditioning blocks `Z1` and
#' `Z2` (together `p_z` variables) are a low-rank signal plus noise,
#' `Z = F A + E`. Both `X` and `Y` are generated from `rank` linear directions of the
#' *observed* `Z`, `S = Z U`, plus their own noise, so with `effect = 0` the null
#' `Y _||_ X | (Z1, Z2)` holds exactly. With `effect > 0` each column of `Y` also
#' receives `effect` times one column of `X` (a direct X -> Y effect), so the null is
#' false.
#'
#' @param n number of rows.
#' @param p_z total number of variables in `Z1` and `Z2` (split evenly; use an even
#'   number).
#' @param q_x,q_y number of variables in `X` and `Y`.
#' @param rank number of directions of Z that drive X and Y (the signal rank).
#' @param noise noise distribution (all rescaled to unit variance):
#'   `"gaussian"`, `"t3"` (heavy tails, Student t with 3 df), `"laplace"`,
#'   `"skewed"` (centred exponential) or `"hetero"` (Gaussian whose standard
#'   deviation depends on the first signal direction of Z; the null still holds).
#' @param effect strength of the direct X -> Y effect (0 gives the null).
#' @param signal which directions of Z drive X and Y: `"factor"` (default; the
#'   directions that carry the low-rank factors of Z, which dominate the variance of
#'   Z and are easy to find even when `p_z > n`) or `"dense"` (random directions in
#'   the whole of Z, mostly noise, which are hard to estimate when `p_z` is large).
#' @return named list of matrices `Y`, `X`, `Z1`, `Z2` with `n` rows.
#' @examples
#' d <- sim_scale(n = 1000, p_z = 20, noise = "t3")
#' sapply(d, dim)
#' @export
sim_scale <- function(n = 1000, p_z = 20, q_x = 4, q_y = 4, rank = 3,
                      noise = c("gaussian", "t3", "laplace", "skewed", "hetero"),
                      effect = 0, signal = c("factor", "dense")) {
  noise <- match.arg(noise)
  signal <- match.arg(signal)
  if (p_z < 2 || p_z %% 2 != 0) stop("p_z must be an even number >= 2")
  # unit-variance noise vector of length m
  draw <- function(m) switch(noise,
    gaussian = stats::rnorm(m),
    hetero   = stats::rnorm(m),
    t3       = stats::rt(m, df = 3) / sqrt(3),
    laplace  = { u <- stats::runif(m) - 0.5; -sign(u) * log(1 - 2 * abs(u)) / sqrt(2) },
    skewed   = stats::rexp(m) - 1)

  Fm <- matrix(stats::rnorm(n * rank), n, rank)
  A  <- matrix(stats::rnorm(rank * p_z), rank, p_z)
  Z  <- Fm %*% A + matrix(draw(n * p_z), n, p_z)               # low-rank signal + noise
  U  <- if (signal == "factor") t(A) else matrix(stats::rnorm(p_z * rank), p_z, rank)
  U  <- sweep(U, 2, sqrt(colSums(U^2)), "/")
  S  <- Z %*% U                                                # directions of the observed Z
  S  <- sweep(S, 2, apply(S, 2, stats::sd), "/")
  Bx <- matrix(stats::rnorm(rank * q_x), rank, q_x) / sqrt(rank)
  By <- matrix(stats::rnorm(rank * q_y), rank, q_y) / sqrt(rank)

  # heteroskedastic noise: scale depends on Z only, so X _||_ Y | Z still holds
  sc <- if (noise == "hetero") sqrt(0.2 + 0.8 * S[, 1]^2) else 1
  X <- S %*% Bx + matrix(draw(n * q_x), n, q_x) * sc
  Y <- S %*% By + matrix(draw(n * q_y), n, q_y) * sc
  if (effect != 0) Y <- Y + effect * X[, (seq_len(q_y) - 1) %% q_x + 1, drop = FALSE]

  h <- p_z / 2
  list(Y = Y, X = X, Z1 = Z[, seq_len(h), drop = FALSE], Z2 = Z[, h + seq_len(h), drop = FALSE])
}
