gen_data_sim <- function(n, tau, alpha_0, alpha_x1, alpha_x2, alpha_u,
                         eta_intercept_a0, eta_intercept_a1, eta_x1, eta_x2, eta_u,
                         k, pU) {
  # ---- Covariate ----
  X1 <- rbinom(n, size = 1, prob = 0.4)
  X2 <- rnorm(n, mean = 0, sd = 1)
  
  #unmeasured confounder
  U <- rbinom(n, size = 1, prob = pU)
  
  # ---- Treatment assignment ----
  linpred_A <- alpha_0 + alpha_x1 * X1 + alpha_x2 * X2 + alpha_u * U
  pi_A      <- plogis(linpred_A)
  A         <- rbinom(n, size = 1, prob = pi_A)
  
  # ---- Event times (Weibull) ----
  #generate via inverse CDF
  #T~weibull(k, b), F(t)=1-exp(- b *t^k), where b = exp(-k* lp)
  #T=(-log(U)/b)^(1/k)
  mu1 <- eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
  mu0 <- eta_intercept_a0 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
  b1 <- exp(-k * mu1)
  b0 <- exp(-k * mu0)
  
  U0 <- runif(n)
  U1 <- runif(n)
  D_a0 <- (-log(U0) / b0)^(1 / k)
  D_a1 <- (-log(U1) / b1)^(1 / k)
  
  # ---- Noninformative censoring ----
  C_a0 <- runif(n, 0.1, tau)
  C_a1 <- runif(n, 0.1, tau)
  
  # ---- Observed time & event indicator ----
  Time <- ifelse(A == 1, pmin(tau, C_a1, D_a1),
              pmin(tau, C_a0, D_a0))
  d <- ifelse(A == 1,
              as.numeric(D_a1 <= pmin(tau, C_a1)),
              as.numeric(D_a0 <= pmin(tau, C_a0)))
  
  # ---- Output dataset ----
  data_sim <- data.frame(Time = Time, d = d, A = A, X1 = X1, X2=X2)
}