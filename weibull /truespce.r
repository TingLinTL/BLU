#single X ~ bernoulli(0.4)
#NO U
t0 <- 2
pX1 <- 0.4
k_true <- 2
eta_a0 <- 0.7; eta_a1 <- 0.2; eta_x <- -0.1

mu0_x0 <- eta_a0 + eta_x*0
mu0_x1 <- eta_a0 + eta_x*1
mu1_x0 <- eta_a1 + eta_x*0
mu1_x1 <- eta_a1 + eta_x*1

b0_x0 <- exp(-k_true * mu0_x0); b0_x1 <- exp(-k_true * mu0_x1)
b1_x0 <- exp(-k_true * mu1_x0); b1_x1 <- exp(-k_true * mu1_x1)

S0_x0 <- exp(-b0_x0 * t0^k_true); S0_x1 <- exp(-b0_x1 * t0^k_true)
S1_x0 <- exp(-b1_x0 * t0^k_true); S1_x1 <- exp(-b1_x1 * t0^k_true)

true_spce <- (1 - pX1) * (S1_x0 - S0_x0) + pX1 * (S1_x1 - S0_x1)
true_spce


#-0.2874432