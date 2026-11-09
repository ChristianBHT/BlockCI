#' Test of conditional independence X ⫫ Y | Z via PLS residuals and permutation
#'
#' @param formula `Y ~ X | Z` (the conditioning set is required). Several blocks
#'                may be given after the bar, e.g. `Y ~ X | Z1 + Z2`; they are
#'                combined column-wise.
#' @param data    Named list of matrices (or a data.frame) holding the variables
#'                in `formula`. Rows with missing values are dropped.
#' @param ncomp  Number of PLS components used to remove Z, or "auto" (cross-validated).
#' @param max_ncomp Upper limit for ncomp when ncomp = "auto".
#' @param stat    "rv", "cc1" or "dcor".
#' @param resid   "cv" (cross-validated residuals, recommended) or "fitted".
#' @param perm    "simple" (permute all rows) or "conditional"
#'                (permute within groups of observations with similar Z).
#' @param B       Number of permutations.
#' @param group_size Desired group size for conditional permutation.
#' @param segments Number of CV segments.
#' @param seed    Optional seed for reproducibility.
#' @return An object of class "htest".
#' @export
pls_ci_test <- function(formula, data,
                        ncomp = "auto", max_ncomp = 10,
                        stat = c("rv", "cc1", "dcor"),
                        resid = c("cv", "fitted"),
                        perm = c("conditional", "simple"),
                        B = 999, group_size = 10,
                        segments = 10, seed = NULL) {
  
  stat  <- match.arg(stat)
  resid <- match.arg(resid)
  perm  <- match.arg(perm)
  if (!is.null(seed)) set.seed(seed)
  
  ## ---- Input ----
  cc <- .canonical(formula, data)          # parses the formula, drops incomplete rows
  if (is.null(cc$V)) stop("formula must have a conditioning set: Y ~ X | Z")
  X <- cc$X; Y <- cc$Y
  Z <- do.call(cbind, cc$V)                # several Z blocks are combined column-wise
  n <- nrow(X)
  
  if (stat == "dcor" && !requireNamespace("energy", quietly = TRUE))
    stop("The package 'energy' is required for stat = 'dcor'.")
  
  ## ---- Statistikker ----
  rv <- function(A, Bm) {
    Sab <- crossprod(A, Bm)
    sum(Sab^2) / sqrt(sum(crossprod(A)^2) * sum(crossprod(Bm)^2))
  }
  cc1  <- function(A, Bm) stats::cancor(A, Bm)$cor[1]
  dcor <- function(A, Bm) energy::dcor(A, Bm)
  stat_fun <- switch(stat, rv = rv, cc1 = cc1, dcor = dcor)
  
  ## ---- Fjern Z fra en blokk med PLS ----
  remove_Z <- function(M) {
    dat <- data.frame(M = I(M), Z = I(Z))
    k_max <- min(if (identical(ncomp, "auto")) max_ncomp else ncomp,
                 ncol(Z), n - 2)
    fit <- pls::plsr(M ~ Z, data = dat, ncomp = k_max, scale = TRUE,
                     validation = "CV", segments = segments)
    
    # velg antall komponenter
    if (identical(ncomp, "auto")) {
      cv  <- pls::MSEP(fit, estimate = "CV", intercept = FALSE)$val  # [1, respons, komp]
      nc  <- which.min(apply(cv[1, , , drop = FALSE], 3, sum))
    } else {
      nc <- k_max
    }
    
    # prediksjoner -> n x q matrise
    pred <- if (resid == "cv") fit$validation$pred[, , nc]
    else               predict(fit, newdata = dat, ncomp = nc)
    pred <- matrix(pred, nrow = n)
    
    list(res = unclass(M) - pred, ncomp = nc)
  }
  
  rx <- remove_Z(X)
  ry <- remove_Z(Y)
  X_res <- scale(rx$res, scale = FALSE)
  Y_res <- scale(ry$res, scale = FALSE)
  
  ## ---- Gruppering for betinget permutasjon ----
  if (perm == "conditional") {
    k  <- min(ncol(Z), 3)
    pc <- stats::prcomp(Z, scale. = TRUE)$x[, seq_len(k), drop = FALSE]
    n_groups <- max(2, n %/% group_size)
    grp <- stats::kmeans(pc, centers = n_groups, nstart = 5)$cluster
  }
  permute_rows <- function() {
    if (perm == "simple") return(sample.int(n))
    idx <- seq_len(n)
    for (g in unique(grp)) {
      i <- which(grp == g)
      idx[i] <- i[sample.int(length(i))]
    }
    idx
  }
  
  ## ---- Permutasjonstest ----
  obs   <- stat_fun(X_res, Y_res)
  perms <- replicate(B, stat_fun(X_res, Y_res[permute_rows(), , drop = FALSE]))
  p     <- (1 + sum(perms >= obs)) / (B + 1)
  
  structure(list(
    statistic   = c(stat = obs),
    p.value     = p,
    method      = sprintf("PLS residual test for b X ⫫ Y | Z (%s, %s permutation, %s residuals)",
                          stat, perm, resid),
    data.name   = paste(deparse(formula), collapse = ""),
    null.value  = NULL,
    alternative = "X og Y er avhengige gitt Z",
    ncomp       = c(X = rx$ncomp, Y = ry$ncomp),
    perm_stats  = perms
  ), class = "htest")
}
