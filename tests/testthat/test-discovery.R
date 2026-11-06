test_that("simulators return blocks and truth", {
  set.seed(1)
  a <- sim_multiview(n = 100)
  expect_length(a$blocks, 3); expect_equal(dim(a$truth), c(3, 3))
  b <- sim_environments(n_env = 3, n = 50, p = 4)
  expect_length(b$blocks, 3); expect_equal(ncol(b$blocks[[1]]), 4)
})

test_that("discovery methods return adjacency matrices; evaluation works", {
  set.seed(1)
  a <- sim_multiview(n = 500)
  est <- discover_multiview(a$blocks)
  expect_equal(dim(est), c(3, 3))
  ev <- evaluate_graph(est, a$truth)
  expect_named(ev, c("precision", "recall", "skeleton_errors", "shd"))
  b <- sim_environments(n_env = 3, n = 150, p = 4)
  expect_equal(dim(discover_environments(b$blocks)), c(4, 4))
})
