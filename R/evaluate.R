#' Compare an estimated directed graph with the truth
#' @param est,truth adjacency matrices (`[i, j] != 0` means i -> j)
#' @return named vector: precision, recall, skeleton_errors, shd
#' @export
evaluate_graph <- function(est, truth) {
  est <- (est != 0) * 1; truth <- (truth != 0) * 1
  tp <- sum(est == 1 & truth == 1); fp <- sum(est == 1 & truth == 0); fn <- sum(est == 0 & truth == 1)
  skel <- function(M) ((M + t(M)) > 0) * 1
  c(precision = if (tp + fp) tp / (tp + fp) else NA, recall = if (tp + fn) tp / (tp + fn) else NA,
    skeleton_errors = sum(skel(est) != skel(truth)) / 2, shd = sum(est != truth))
}
