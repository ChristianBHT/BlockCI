# Simulators with known ground truth ------------------------------------------
data1 <- sim_1()
test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data1)
test
data2 <- sim_2()
test <- SOPLS.CCI.test(Y ~ X | Z1 + Z2, data = data2)
test

