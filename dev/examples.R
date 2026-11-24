# # This is a working script for testing out snippets and tests
# library(devtools)
# devtools::load_all("C:/Users/chtho2119/Documents/GitHub/BlockCI")
dat <- sim_1()

Block.resid.test(Y ~ X | Z1 + Z2, data = dat)
Block.CCI.test(Y ~ X | Z1 + Z2, data = dat)

Block.resid.test(Y ~ X | Z1 + Z2, data = dat, stat = "wgcm")
Block.CCI.test(Y ~ X | Z1 + Z2, data = dat, method = "ROSA")

Block.resid.test(Y ~ X |  Z2, data = dat)
Block.CCI.test(Y ~ X |  Z2, data = dat)

Block.resid.test(Y ~ X |  Z2, data = dat, stat = "wgcm")
Block.CCI.test(Y ~ X |  Z2, data = dat, method = "ROSA")

  