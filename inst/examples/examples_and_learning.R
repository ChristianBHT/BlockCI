# Example of RV coefficient calculation
# matrix X
X <- matrix(c(-3, -1, 1, 3, 1, -1, -1, 1), nrow = 4, ncol = 2)
Y <- matrix(c(-2, -1, 2, 1, 1, 1, -1, -1), nrow = 4, ncol = 2)
# X transposed
X_t <- t(X)
# MAtrix multiplication
kryss_cov <- X_t %*% Y
# Kvadrer elementene i matrisen 
kryss_cov_squared <- kryss_cov^2
# summer
sum_kryss_cov_squared <- sum(kryss_cov_squared)

# X var
X_var <- diag(X_t %*% X)^2
kryssryss_var_X_sq <- sum(X_var)


# Y var
Y_var <- t(Y) %*% Y
kryss_var_Y_sq <- Y_var^2
sum_kryss_var_Y_sq <- sum(kryss_var_Y_sq)

RV <- sum_kryss_cov_squared / sqrt(sum_kryss_var_X_sq * sum_kryss_var_Y_sq)

# Canonical correlation
# Covariance matrix R_Y and R_X
stats::cancor(X, Y)

