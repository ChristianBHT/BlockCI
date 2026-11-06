# Repeat Setting A over many seeds to see average accuracy
library(CausalBlocks)
r <- t(sapply(1:30, function(s) {
  set.seed(s); a <- sim_multiview(n = 1000)
  evaluate_graph(discover_multiview(a$blocks, alpha = 0.05), a$truth)
}))
print(round(colMeans(r, na.rm = TRUE), 2))
