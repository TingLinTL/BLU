library(future.apply)
library(rjags)
library(coda)
library(LaplacesDemon)

plan(multisession, workers = parallel::detectCores() - 1)  

#JAGS model 
wb_model_U <- "
model {
  
  # ----- Likelihood, non-informative -----
  for (i in 1:N) {
    # Right-censoring via interval indicator
    is_cens[i] ~ dinterval(T[i], C[i])
    
    # Outcome model: Weibull distribution, dweib(shape = k, rate = b[i])
    T[i] ~ dweib(k, b[i])
    
    # AFT linear predictor with latent U 
    # mu[i] = eta0 + eta_x * X[i] + eta_a * A[i] + eta_u * U[i]
    # log b[i] = -k * mu[i]
    # rate parameter b[i] = exp(-k * mu[i])
    
    mu[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] + eta_u*U[i]
    log_b[i] <- -k * mu[i]
    b[i] <- exp(log_b[i])
    
    # ---- Latent U model ----
    #model U as binary
    U[i] ~ dbern(pU[i])
    logit(pU[i]) <- gamma0 + gamma1*X1[i] + gamma2*X2[i] 
    
    # ---- Treatment assignment model A|X,U ----
    A[i] ~ dbern(pA[i])
    logit(pA[i]) <- alpha0 + alpha1*X1[i] + alpha2*X2[i] + alpha3*U[i] 
  }
  
  # ----- Priors -----
  #outcome parameters
  eta0 ~ dnorm(0, 0.001)
  eta_x1 ~ dnorm(0, 0.001)
  eta_x2 ~ dnorm(0, 0.001)
  eta_a ~ dnorm(0, 0.001)
  eta_u ~ dnorm(0, 0.001)      # sensitivity effect of U on outcome
  
  #U parameters
  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)
  gamma2 ~ dnorm(0, 0.01)
  
  #Treatment parameters
  alpha0 ~ dnorm(0, 0.001)
  alpha1 ~ dnorm(0, 0.001)
  alpha2 ~ dnorm(0, 0.001)
  alpha3 ~ dnorm(0, 0.001)    # sensitivity effect of U on treatment
  
  logk ~ dnorm(0, 0.001)
  k <- exp(logk)

}
"

gen_data_sim <- function(n, tau, alpha_0, alpha_x1, alpha_x2, alpha_u,
                         eta_intercept_a0, eta_intercept_a1, eta_x1, eta_x2, eta_u,
                         k, pU) {
  X1 <- rbinom(n, size = 1, prob = 0.4)
  X2 <- rnorm(n, mean = 0, sd = 1)
  U <- rbinom(n, size = 1, prob = pU)
  linpred_A <- alpha_0 + alpha_x1 * X1 + alpha_x2 * X2 + alpha_u * U
  A <- rbinom(n, 1, plogis(linpred_A))
  mu1 <- eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
  mu0 <- eta_intercept_a0 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
  b1  <- exp(-k * mu1); b0 <- exp(-k * mu0)
  D_a0 <- (-log(runif(n)) / b0)^(1 / k)
  D_a1 <- (-log(runif(n)) / b1)^(1 / k)
  C_a0 <- runif(n, 0.1, tau); C_a1 <- runif(n, 0.1, tau)
  M <- ifelse(A == 1, pmin(tau, C_a1, D_a1), pmin(tau, C_a0, D_a0))
  d <- ifelse(A == 1, as.numeric(D_a1 <= pmin(tau, C_a1)),
              as.numeric(D_a0 <= pmin(tau, C_a0)))
  data.frame(M, d, A, X1, X2)
}

compute_spce_bb_with_U <- function(samp, X1, X2, t0, rx, dirichlet_alpha = 1) {
  post   <- as.matrix(samp)
  eta0    <- post[, "eta0"]
  eta_x1   <- post[, "eta_x1"]
  eta_x2   <- post[, "eta_x2"]
  eta_a   <- post[, "eta_a"]
  eta_u  <- post[,"eta_u"]
  k_draws <- post[, "k"]
  gamma0  <- post[, "gamma0"]
  gamma1  <- post[, "gamma1"]
  gamma2  <- post[, "gamma2"]
  M <- nrow(post)
  L <- length(X1)
  S0_marg <- numeric(M); S1_marg <- numeric(M)
  for (m in 1:M) {
    w   <- as.numeric(LaplacesDemon::rdirichlet(1, rep(dirichlet_alpha, L)))
    idx <- sample.int(L, size = rx, replace = TRUE, prob = w)  # sample rows
    X1_star  <- X1[idx]
    X2_star  <- X2[idx]
    pU     <- plogis(gamma0[m] + gamma1[m] * X1_star + gamma2[m] * X2_star)
    U_star <- rbinom(rx, 1, pU)
    mu0 <- eta0[m] + eta_x1[m] * X1_star + eta_x2[m] * X2_star + eta_u[m] * U_star
    mu1 <- mu0 + eta_a[m]
    b0  <- exp(-k_draws[m] * mu0); b1 <- exp(-k_draws[m] * mu1)
    t0k <- t0 ^ k_draws[m]
    S0_marg[m] <- mean(exp(-b0 * t0k))
    S1_marg[m] <- mean(exp(-b1 * t0k))
  }
  list(S0_marg = S0_marg, S1_marg = S1_marg)
}

simulate_once <- function(s, n, tau, t0, true_spce_aft, rx) {
  set.seed(1000 + s)  # per-worker reproducibility
  dat <- gen_data_sim(n = n, tau = tau,
                      alpha_0 = 0.1, alpha_x1 = 0.3, alpha_x2 = 0.7, alpha_u = -0.8,
                      eta_intercept_a0 = 0.7, eta_intercept_a1 = 0.2, eta_x1 = -0.1, eta_x2 = 0.4, eta_u = -0.8,
                      k = 2, pU = 0.5)
  
  T_vec <- dat$M; T_vec[dat$d == 0] <- NA_real_
  eps <- 1e-6; C_vec <- ifelse(dat$d == 1L, dat$M + eps, dat$M)
  is_cens <- 1L - dat$d
  jags_data <- list(N = nrow(dat), T = T_vec, C = C_vec,
                    is_cens = is_cens, A = dat$A, X1 = dat$X1, X2 = dat$X2)
  
  jm <- jags.model(textConnection(wb_model_U), data = jags_data,
                   n.chains = 1)
  params <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","k","gamma0","gamma1","gamma2")
  samp_coda <- coda.samples(jm, params, n.iter = 20000)
  samp = data.frame(samp_coda[[1]][10001:20000, ])
  
  surv <- compute_spce_bb_with_U(samp, dat$X1, dat$X2, t0, rx, dirichlet_alpha = 1)
  spce_draws <- surv$S1_marg - surv$S0_marg
  spce_mean  <- mean(spce_draws)
  spce_ci    <- quantile(spce_draws, c(0.025, 0.975))
  
  list(
    pos_mean = spce_mean,
    sd_draws = sd(spce_draws),
    bias     = spce_mean - true_spce_aft,
    lo       = spce_ci[1],
    hi       = spce_ci[2]
  )
}

# ---- Parallel runs ----
S  <- 18
n  <- 1000
tau <- 5.5
t0  <- 2
true_spce_aft <- -0.1384194
rx <- 10000 #try small first

res_list <- future_lapply(
  X = seq_len(S),
  FUN = simulate_once,
  n = n, tau = tau, t0 = t0, true_spce_aft = true_spce_aft, rx = rx,
  future.seed = TRUE    # reproducible parallel RNG
)

# Collect and summarize
res <- do.call(rbind, lapply(res_list, function(x) as.data.frame(x)))
ESE <- sd(res$pos_mean)
ASD <- mean(res$sd_draws)
CP  <- 100 * mean(true_spce_aft >= res$lo & true_spce_aft <= res$hi)
summary_df <- data.frame(
  avg_pos_mean = mean(res$pos_mean),
  ESE = ESE, ASD = ASD,
  avg_bias = mean(res$bias),
  CP = CP
)
print(summary_df)
