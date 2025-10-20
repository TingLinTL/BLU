# install.packages("statmod")  # if needed
library(statmod)

true_spce_weibull_aft <- function(t0,
                                  k,
                                  eta_intercept_a0, eta_intercept_a1,
                                  eta_x1, eta_x2, eta_u,
                                  pX1 = 0.4,
                                  pU  = 0.5,
                                  gh_n = 50) {
  # Gauss–Hermite nodes/weights for standard Normal X2 ~ N(0,1)
  gh <- statmod::gauss.quad.prob(gh_n, dist = "normal")
  z  <- gh$nodes     # X2 support points
  w  <- gh$weights   # weights summing to 1
  t0k <- t0^k
  
  # Survival helper: S = exp(-exp(-k*mu)*t0^k)
  S_from_mu <- function(mu) exp( - exp(-k * mu) * t0k )
  
  # For a given x1 in {0,1}, compute E_{X2}[ E_{U|X} S_a ]
  eval_for_x1 <- function(x1) {
    # baseline (U=0) linear predictors across GH nodes z
    mu0_u0 <- eta_intercept_a0 + eta_x1 * x1 + eta_x2 * z
    mu1_u0 <- eta_intercept_a1 + eta_x1 * x1 + eta_x2 * z
    # add U effect
    mu0_u1 <- mu0_u0 + eta_u
    mu1_u1 <- mu1_u0 + eta_u
    
    # E_{U|X} S_a = (1-pU) S(mu_a, U=0) + pU S(mu_a, U=1)
    S0_mix <- (1 - pU) * S_from_mu(mu0_u0) + pU * S_from_mu(mu0_u1)
    S1_mix <- (1 - pU) * S_from_mu(mu1_u0) + pU * S_from_mu(mu1_u1)
    
    # Integrate over X2 ~ N(0,1) via GH
    c(S0 = sum(w * S0_mix), S1 = sum(w * S1_mix))
  }
  
  # Sum over X1 in {0,1}
  out_x0 <- eval_for_x1(0)
  out_x1 <- eval_for_x1(1)
  
  EX_S0 <- (1 - pX1) * out_x0["S0"] + pX1 * out_x1["S0"]
  EX_S1 <- (1 - pX1) * out_x0["S1"] + pX1 * out_x1["S1"]
  
  SPCE  <- as.numeric(EX_S1 - EX_S0)
  return(SPCE)
}


t0 <- 2
k  <- 2
eta_intercept_a0 <- 0.7
eta_intercept_a1 <- 0.2
eta_x1 <- -0.1
eta_x2 <-  0.4
eta_u  <- -0.8

true_spce <- true_spce_weibull_aft(
  t0, k,
  eta_intercept_a0, eta_intercept_a1,
  eta_x1, eta_x2, eta_u,
  pX1 = 0.4, pU = 0.5,
  gh_n = 50   # 30–80 is plenty; 50 is very accurate and fast
)
true_spce #[1] -0.1384194


####################
# Weibull–AFT survival S(t|mu) = exp( - exp(-k*mu) * t^k )
S_from_mu <- function(mu, k, t0k) exp(-exp(-k * mu) * t0k)

# E_U S(t0 | a, x1, x2) with Bernoulli U, P(U=1)=pU (constant)
EU_S_given_X <- function(a, x1, x2, k, t0k,
                         eta0_a0, eta0_a1, eta_x1, eta_x2, eta_u, pU) {
  mu0_u0 <- eta0_a0 + eta_x1 * x1 + eta_x2 * x2
  mu1_u0 <- eta0_a1 + eta_x1 * x1 + eta_x2 * x2
  if (a == 0) {
    (1 - pU) * S_from_mu(mu0_u0, k, t0k) +
      pU     * S_from_mu(mu0_u0 + eta_u, k, t0k)
  } else {
    (1 - pU) * S_from_mu(mu1_u0, k, t0k) +
      pU     * S_from_mu(mu1_u0 + eta_u, k, t0k)
  }
}

# True SPCE via integrate() over X2 ~ N(0,1) and sum over X1
true_spce_integral <- function(t0, k,
                               eta_intercept_a0, eta_intercept_a1,
                               eta_x1, eta_x2, eta_u,
                               pX1 = 0.4, pU = 0.5,
                               rel.tol = 1e-8, abs.tol = 0) {
  t0k <- t0^k
  # integrand for X2: E_U S(t0 | a, x1, x2) * phi(x2)
  int_a_x1 <- function(a, x1) {
    f <- function(x2) {
      EU_S_given_X(a, x1, x2, k, t0k,
                   eta_intercept_a0, eta_intercept_a1,
                   eta_x1, eta_x2, eta_u, pU) * dnorm(x2)
    }
    integrate(f, lower = -Inf, upper = Inf,
              rel.tol = rel.tol, abs.tol = abs.tol, subdivisions = 200L)$value
  }
  
  EX2_S0_x0 <- int_a_x1(0, 0); EX2_S1_x0 <- int_a_x1(1, 0)
  EX2_S0_x1 <- int_a_x1(0, 1); EX2_S1_x1 <- int_a_x1(1, 1)
  
  EX_S0 <- (1 - pX1) * EX2_S0_x0 + pX1 * EX2_S0_x1
  EX_S1 <- (1 - pX1) * EX2_S1_x0 + pX1 * EX2_S1_x1
  EX_S1 - EX_S0
}
t0 <- 2; k <- 2
eta_intercept_a0 <- 0.7; eta_intercept_a1 <- 0.2
eta_x1 <- -0.1; eta_x2 <- 0.4; eta_u <- -0.8
true_spce <- true_spce_integral(t0, k,
                                eta_intercept_a0, eta_intercept_a1,
                                eta_x1, eta_x2, eta_u,
                                pX1 = 0.4, pU = 0.5)
true_spce


