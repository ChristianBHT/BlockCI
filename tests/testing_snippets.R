# This is a working script for testing out snippets and tests
library(devtools)
devtools::load_all("C:/Users/chtho2119/Documents/GitHub/BlockCI")
dat <- sim_1()

pls_ci_test(Y ~ X | Z1 + Z2, data = dat, stat = "rv", resid = "cv", perm = "conditional", B = 99, group_size = 10, segments = 5, seed = 123)
pls_ci_test(Y ~ X | Z1, data = dat, stat = "rv", resid = "cv", perm = "conditional", B = 99, group_size = 10, segments = 5, seed = 123)
SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = dat)
SOPLS.CCI.test(Y ~ X | Z1, data = dat)
SOPLS.CCI.fast(Y ~ X | Z1 + Z2, data = dat)
SOPLS.CCI.fast(Y ~ X | Z1, data = dat)
