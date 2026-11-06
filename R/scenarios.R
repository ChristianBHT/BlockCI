# Null scenarios: four blocks Y, X, Z1, Z2 with Y _||_ X | (Z1, Z2) -----------
# X and Y are generated from the OBSERVED Z blocks plus independent noise, so the
# null holds exactly (no latent confounding left outside Z).

.mix <- function(F, p, noise_sd, noise = stats::rnorm) {
  L <- matrix(stats::runif(ncol(F) * p, 0.5, 1.2) * sample(c(-1, 1), ncol(F) * p, TRUE), ncol(F), p)
  S <- scale(F %*% L)                                   # unit-variance signal columns
  S + matrix(noise(nrow(F) * p), nrow(F), p) * noise_sd
}

lap <- function(m) { u <- stats::runif(m) - 0.5; -sign(u) * log(1 - 2 * abs(u)) }

sim_simple <- function(n = 500, p = 4) {
  Z1 <- matrix(rnorm(n * p), n, p)
  Z2 <- matrix(rnorm(n * p), n, p)
  Z  <- cbind(Z1, Z2)
  X  <- Z %*% matrix(rnorm(2 * p * p), 2 * p, p) + matrix(rnorm(n * p), n, p)
  Y  <- Z %*% matrix(rnorm(2 * p * p), 2 * p, p) + matrix(rnorm(n * p), n, p)
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
}


sim_1 <- function(n = 500, p = c(Y = 4, X = 4, Z1 = 5, Z2 = 5)) {
  f1 <- matrix(stats::rnorm(n * 2), n) 
  f2 <- matrix(stats::rnorm(n * 2), n)
  Z1 <- .mix(f1, p[["Z1"]], 0.5) 
  Z2 <- .mix(f2, p[["Z2"]], 0.5)
  Z <- cbind(Z1, Z2)
  X <- .mix(Z, p[["X"]], 0.7) 
  Y <- .mix(Z, p[["Y"]], 0.7)
  return(list(Y = Y, X = X, Z1 = Z1, Z2 = Z2))
  } 
  
sim_2 <- function(n = 500, p = c(Y = 4, X = 4, Z1 = 5, Z2 = 5)) {
    Z1 <- .mix(f1, p[["Z1"]], 0.5)
    Z2 <- .mix(cbind(Z1, f2), p[["Z2"]], 0.5)           # Z1 -> Z2
    X <- .mix(Z1, p[["X"]], 0.7) 
    Y <- .mix(Z2, p[["Y"]], 0.7)
    return(list(Y = Y, X = X, Z1 = Z1, Z2 = Z2))
  } 

sim_3 <- function(){
    Z1 <- .mix(f1, p[["Z1"]], 0.3) 
    Z2 <- .mix(f2, p[["Z2"]], 0.3)
    Z <- cbind(Z1, Z2)
    X <- .mix(Z, p[["X"]], 0.1)
    Y <- .mix(Z, p[["Y"]], 0.7)
    return(list(Y = Y, X = X, Z1 = Z1, Z2 = Z2))
} 


sim_4 <- function(n = 500, p = c(Y = 4, X = 4, Z1 = 5, Z2 = 5)){
  f1 <- matrix(stats::rnorm(n * 2), n) 
  f2 <- matrix(stats::rnorm(n * 2), n)
  lap <- function(m) { 
      u <- stats::runif(m) - 0.5 - sign(u) * log(1 - 2 * abs(u)) 
      }
  Z1 <- .mix(f1, p[["Z1"]], 0.5, lap) 
  Z2 <- .mix(f2, p[["Z2"]], 0.5, lap)
  Z <- scale(cbind(Z1, Z2))
  W <- matrix(stats::rnorm(ncol(Z) * 3), ncol(Z), 3) 
  W <- sweep(W, 2, sqrt(colSums(W^2)), "/")
  u <- Z %*% W
  g <- cbind(tanh(u[, 1]) + u[, 2]^2 / 2, u[, 1] * u[, 3], tanh(u[, 3]) - u[, 2]^2 / 2)
  het <- sqrt(0.3 + 0.5 * abs(u[, 1]))                # noise scale depends on Z
  X <- .mix(g, p[["X"]], 0.6, lap) * het
  Y <- .mix(g[, 3:1], p[["Y"]], 0.6, lap) * het
  return(list(Y = Y, X = X, Z1 = Z1, Z2 = Z2))
  }

  
  
  
  
  
  
  
  
  