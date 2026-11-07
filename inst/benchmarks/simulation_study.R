# Simulators with known ground truth ------------------------------------------
undebug(sopls_test_gen)

test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data1, nperm = 100)
test
data1 <- sim_simple(250)
test <- SOPLS.CCI.test(Y ~ X | Z1, data = data1, nperm = 200, pvalue = "parametric")
plot(test)

tic <- Sys.time()
test1_typeI <- list()
test1_power <- list()
for (i in c(1:100)){
  data <- sim_simple(1000)
  test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data, nperm = 75, pvalue = "parametric")
  test1_typeI[i] <- test[2]
  test <- SOPLS.CCI.test(Y ~ X | Z1, data = data, nperm = 75, pvalue = "parametric")
  test1_power[i] <- test[2]
}
toc <- Sys.time()


test2_typeI <- list()
test2_power <- list()
for (i in c(1:100)){
  data <- sim_1(1000)
  test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data, nperm = 75, pvalue = "parametric")
  test2_typeI[i] <- test[2]
  test <- SOPLS.CCI.test(Y ~ X | Z1, data = data, nperm = 75, pvalue = "parametric")
  test2_power[i] <- test[2]
}

test3_typeI <- list()
test3_power <- list()
for (i in c(1:100)){
  data <- sim_2(1000)
  test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data, nperm = 75, pvalue = "parametric")
  test3_typeI[i] <- test[2]
  test <- SOPLS.CCI.test(Y ~ X | Z1, data = data, nperm = 75, pvalue = "parametric")
  test3_power[i] <- test[2]
}

test4_typeI <- list()
test4_power <- list()
for (i in c(1:100)){
  data <- sim_3(1000)
  test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data, nperm = 75, pvalue = "parametric")
  test4_typeI[i] <- test[2]
  test <- SOPLS.CCI.test(Y ~ X | Z1, data = data, nperm = 75, pvalue = "parametric")
  test4_power[i] <- test[2]
}


test5_typeI <- list()
test5_power <- list()
for (i in c(1:100)){
  data <- sim_3(1000)
  test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data, nperm = 75, pvalue = "parametric")
  test5_typeI[i] <- test[2]
  test <- SOPLS.CCI.test(Y ~ X | Z1, data = data, nperm = 75, pvalue = "parametric")
  test5_power[i] <- test[2]
}

print(paste0("Time for 1000 simple test", toc - tic,""))

