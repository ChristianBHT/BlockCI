# Null scenarios: four blocks Y, X, Z1, Z2 with Y _||_ X | (Z1, Z2) -----------
# X and Y are generated from the OBSERVED Z blocks plus independent noise, so the
# null holds exactly (no latent confounding left outside Z).

.mix <- function(F, p, noise_sd, noise = stats::rnorm) {
  L <- matrix(stats::runif(ncol(F) * p, 0.5, 1.2) * sample(c(-1, 1), ncol(F) * p, TRUE), ncol(F), p)
  S <- scale(F %*% L)                                   # unit-variance signal columns
  S + matrix(noise(nrow(F) * p), nrow(F), p) * noise_sd
}

# Laplace noise
.lap <- function(m) { u <- stats::runif(m) - 0.5; -sign(u) * log(1 - 2 * abs(u)) }

sim_1 <- function(n = 500, p = 6) {
  Z1 <- matrix(stats::rnorm(n * p), n, p)
  Z2 <- matrix(stats::rnorm(n * p), n, p)
  Z  <- cbind(Z1, Z2)
  
  X  <- Z %*% matrix(stats::rnorm(2 * p * p), 2 * p, p) + matrix(stats::rnorm(n * p), n, p)
  Y  <- Z %*% matrix(stats::rnorm(2 * p * p), 2 * p, p) + matrix(stats::rnorm(n * p), n, p)
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
}

sim_2 <- function(n = 500, p = c(Y = 4, X = 4, Z1 = 5, Z2 = 5)) {
  f1 <- matrix(stats::rnorm(n * 2), n)
  f2 <- matrix(stats::rnorm(n * 2), n)
  Z1 <- .mix(f1, p[["Z1"]], 0.5)
  Z2 <- .mix(f2, p[["Z2"]], 0.5)
  Z <- cbind(Z1, Z2)
  X <- .mix(Z, p[["X"]], 0.7)
  Y <- .mix(Z, p[["Y"]], 0.7)
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
}

sim_3 <- function(n = 500, p = c(Y = 4, X = 4, Z1 = 5, Z2 = 5)) {
  f1 <- matrix(stats::rnorm(n * 2), n)
  f2 <- matrix(stats::rnorm(n * 2), n)
  Z1 <- .mix(f1, p[["Z1"]], 0.5)
  Z2 <- .mix(cbind(Z1, f2), p[["Z2"]], 0.5)           # Z1 -> Z2
  X <- .mix(Z1, p[["X"]], 0.7)
  Y <- .mix(Z2, p[["Y"]], 0.7)
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
}

sim_4 <- function(n = 500, p = c(Y = 4, X = 4, Z1 = 5, Z2 = 5)) {
  f1 <- matrix(stats::rnorm(n * 2), n)
  f2 <- matrix(stats::rnorm(n * 2), n)
  Z1 <- .mix(f1, p[["Z1"]], 0.3)
  Z2 <- .mix(f2, p[["Z2"]], 0.3)
  Z <- cbind(Z1, Z2)
  X <- .mix(Z, p[["X"]], 0.1)                         
  Y <- .mix(Z, p[["Y"]], 0.7)
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
}


sim_5 <- function(n = 500, p = c(Y = 4, X = 4, Z1 = 5, Z2 = 5)) {
  f1 <- matrix(stats::rnorm(n * 2), n)
  f2 <- matrix(stats::rnorm(n * 2), n)
  Z1 <- .mix(f1, p[["Z1"]], 0.5, .lap)
  Z2 <- .mix(f2, p[["Z2"]], 0.5, .lap)
  Z <- scale(cbind(Z1, Z2))
  W <- matrix(stats::rnorm(ncol(Z) * 3), ncol(Z), 3)
  W <- sweep(W, 2, sqrt(colSums(W^2)), "/")
  u <- Z %*% W
  g <- cbind(tanh(u[, 1]) + u[, 2]^2 / 2, u[, 1] * u[, 3], tanh(u[, 3]) - u[, 2]^2 / 2)
  het <- sqrt(0.3 + 0.5 * abs(u[, 1]))                # noise scale depends on Z
  X <- .mix(g, p[["X"]], 0.6, .lap) * het
  Y <- .mix(g[, 3:1], p[["Y"]], 0.6, .lap) * het
  list(Y = Y, X = X, Z1 = Z1, Z2 = Z2)
}


sim_6 <- function(n = 500, 
                  z1s = 4, 
                  z2s = 4, 
                  Xs = 4, 
                  Ys = 4,
                  sd = 1
                  ){
  noise1 <- rnorm(n) 
  Z1 <- sapply(1:z1s, function(i) noise1 + rnorm(n, sd = sd))
  Z2  <- sapply(1:z2s, function(i) noise1 + rnorm(n, sd = sd))
  
  Z <- cbind(Z1, Z2)
  f2 <- rowMeans(Z)
  X = sapply(1:Xs, function(i) f2 + rnorm(n, Z, sd))
  Y = sapply(1:Ys, function(i) rnorm(n) + rnorm(n, Z, sd))
  
  return(list('Z1' = Z1, 'Z2' = Z2, 'X' = X, 'Y' = Y))
}


sim_7 <- function(n, d = 0, sd = 0.1){
  Z_1 <- rnorm(n, 0, sd)
  Z_2 <- rnorm(n, 0, sd)
  
  z_1 <- rnorm(n, Z_1, sd) 
  z_2 <- rnorm(n, Z_1, sd)
  z_3 <- rnorm(n, Z_1 + Z_2, sd) 
  z_4 <- rnorm(n, Z_1 + Z_2, sd)
  z_5 <- rnorm(n, Z_2, sd) 
  z_6 <- rnorm(n, Z_2, sd)
  
  X <- rnorm(n, z_1 + z_2 + z_3 + z_4, sd)
  Y <- rnorm(n, z_3 + z_4 + z_5 + z_6, sd)
  
  x_1 <- rnorm(n, X, sd)
  x_2 <- rnorm(n, X, sd)  
  x_3 <- rnorm(n, X, sd)
  x_4 <- rnorm(n, X, sd)
  y_1 <- rnorm(n, Y + d*X, sd)
  y_2 <- rnorm(n, Y, sd)
  y_3 <- rnorm(n, Y, sd)
  y_4 <- rnorm(n, Y, sd)
  
  return(data.frame(z_1, z_2, z_3, z_4, z_5, z_6, x_1, x_2, x_3, x_4, y_1, y_2, y_3, y_4))
}



