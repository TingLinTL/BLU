## ---- Parameters----
t0      <- 2
k_true  <- 2
eta_a0  <- 0.7
eta_a1  <- 0.2
eta_x   <- -0.1
eta_u   <- -0.8
pX1     <- 0.4
pU1     <- 0.5   

# Survival under arm a given (x,u)
S_weib <- function(t, k, eta_a, x, u, eta_x, eta_u){
  mu <- eta_a + eta_x * x + eta_u * u
  b  <- exp(-k * mu)
  exp(-b * t^k)
}

## Enumerate (x,u) belongs to {0,1}×{0,1} and average
pX0 <- 1 - pX1
pU0 <- 1 - pU1

w_x0u0 <- pX0 * pU0
w_x0u1 <- pX0 * pU1
w_x1u0 <- pX1 * pU0
w_x1u1 <- pX1 * pU1

S1_x0u0 <- S_weib(t0, k_true, eta_a1, 0, 0, eta_x, eta_u)
S0_x0u0 <- S_weib(t0, k_true, eta_a0, 0, 0, eta_x, eta_u)

S1_x0u1 <- S_weib(t0, k_true, eta_a1, 0, 1, eta_x, eta_u)
S0_x0u1 <- S_weib(t0, k_true, eta_a0, 0, 1, eta_x, eta_u)

S1_x1u0 <- S_weib(t0, k_true, eta_a1, 1, 0, eta_x, eta_u)
S0_x1u0 <- S_weib(t0, k_true, eta_a0, 1, 0, eta_x, eta_u)

S1_x1u1 <- S_weib(t0, k_true, eta_a1, 1, 1, eta_x, eta_u)
S0_x1u1 <- S_weib(t0, k_true, eta_a0, 1, 1, eta_x, eta_u)

true_spce_U <- 
  w_x0u0 * (S1_x0u0 - S0_x0u0) +
  w_x0u1 * (S1_x0u1 - S0_x0u1) +
  w_x1u0 * (S1_x1u0 - S0_x1u0) +
  w_x1u1 * (S1_x1u1 - S0_x1u1)

true_spce_U #-0.1464996
