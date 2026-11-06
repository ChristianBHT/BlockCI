# Simulators with known ground truth ------------------------------------------
undebug(sopls_test_gen)
data1 <- sim_simple(200)
test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data1, nperm = 10)

res <- list()
for (i in c(1:100)){
  data1 <- sim_simple(200)
  test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data1, nperm = 100)
  res[[i]] <- test[2]
}

  test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data1)
test
data2 <- sim_2()
test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data2)
test

