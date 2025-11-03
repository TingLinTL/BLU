t0 <- Sys.time()
library(doParallel)
library(doRNG)            # reproducible foreach RNG streams
library(rjags)
library(MCMCprecision)

## --- data generator function (AFT Weibull) ---
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
  data_sim_noU <- data.frame(Time = Time, d = d, A = A, X1 = X1, X2=X2)
  data_sim_withU <- data.frame(Time = Time, d = d, A = A, X1 = X1, X2=X2, U=U)
  
  list(data_sim_noU = data_sim_noU, data_sim_withU = data_sim_withU)
}

#--------------- funtion to compute spce with U-------------
#compute spce for latU approach
compute_spce_bb_latU <- function(
    samp,        # posterior draws; must have: eta0, eta_x, eta_a, k, gamma0, gamma1
    X1, 
    X2,          # observed covariate vector
    t0,          # time at which to evaluate survival
    rx,          # size of BB resample X* B
    dirichlet_alpha = 1)  # BB prior mass
{
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
  
  S0_marg  <- numeric(M)
  S1_marg  <- numeric(M)
  
  for (m in 1:M) {
    #BB over X -> X*
    w <- as.numeric(MCMCprecision::rdirichlet(1, rep(dirichlet_alpha, L)))
    idx <- sample.int(L, size = rx, replace = TRUE, prob = w)  # sample rows
    X1_star  <- X1[idx]
    X2_star  <- X2[idx]
    
    #U* | X1*, X2*
    pU     <- plogis(gamma0[m] + gamma1[m] * X1_star + gamma2[m] * X2_star)
    U_star <- rbinom(rx, size = 1, prob = pU)
    
    #AFT linear predictors
    mu0 <- eta0[m] + eta_x1[m] * X1_star + eta_x2[m] * X2_star + eta_u[m] * U_star
    mu1 <- mu0 + eta_a[m]
    
    #Weibull survival at t0 and average over X*
    b0  <- exp(-k_draws[m] * mu0)
    b1  <- exp(-k_draws[m] * mu1)
    t0k <- t0 ^ k_draws[m]
    
    S0_marg[m] <- mean(exp(-b0 * t0k))
    S1_marg[m] <- mean(exp(-b1 * t0k))
  }
  return(list(S0_marg = S0_marg, S1_marg = S1_marg))
}

#compute spce for observed U approach
compute_spce_bb_obsU <- function(
    samp,        # posterior draws; must have: eta0, eta_x, eta_a, k, gamma0, gamma1
    X1, 
    X2,          # observed covariate vector
    U,           # including unobserved covariate vector
    t0,          # time at which to evaluate survival
    rx,          # size of BB resample X* B
    dirichlet_alpha = 1)  # BB prior mass
{
  post   <- as.matrix(samp)
  eta0    <- post[, "eta0"]
  eta_x1   <- post[, "eta_x1"]
  eta_x2   <- post[, "eta_x2"]
  eta_a   <- post[, "eta_a"]
  eta_u  <- post[,"eta_u"]
  k_draws <- post[, "k"]
  M <- nrow(post)
  L <- length(X1)
  
  S0_marg  <- numeric(M)
  S1_marg  <- numeric(M)
  
  for (m in 1:M) {
    #BB over X -> X*
    w <- as.numeric(MCMCprecision::rdirichlet(1, rep(dirichlet_alpha, L)))
    idx <- sample.int(L, size = rx, replace = TRUE, prob = w)  # sample rows
    X1_star  <- X1[idx]
    X2_star  <- X2[idx]
    U_star   <- U[idx]
    
    #AFT linear predictors
    mu0 <- eta0[m] + eta_x1[m] * X1_star + eta_x2[m] * X2_star + eta_u[m] * U_star
    mu1 <- mu0 + eta_a[m]
    
    #Weibull survival at t0 and average over X*
    b0  <- exp(-k_draws[m] * mu0)
    b1  <- exp(-k_draws[m] * mu1)
    t0k <- t0 ^ k_draws[m]
    
    S0_marg[m] <- mean(exp(-b0 * t0k))
    S1_marg[m] <- mean(exp(-b1 * t0k))
  }
  return(list(S0_marg = S0_marg, S1_marg = S1_marg))
}

#compute spce for naive approach
compute_spce_bb_naive <- function(
    samp,        # posterior draws; must have: eta0, eta_x, eta_a, k, gamma0, gamma1
    X1, 
    X2,          # observed covariate vector
    t0,          # time at which to evaluate survival
    rx,          # size of BB resample X* B
    dirichlet_alpha = 1)  # BB prior mass
{
  post   <- as.matrix(samp)
  eta0    <- post[, "eta0"]
  eta_x1   <- post[, "eta_x1"]
  eta_x2   <- post[, "eta_x2"]
  eta_a   <- post[, "eta_a"]
  k_draws <- post[, "k"]
  M <- nrow(post)
  L <- length(X1)
  
  S0_marg  <- numeric(M)
  S1_marg  <- numeric(M)
  
  for (m in 1:M) {
    #BB over X -> X*
    w <- as.numeric(MCMCprecision::rdirichlet(1, rep(dirichlet_alpha, L)))
    idx <- sample.int(L, size = rx, replace = TRUE, prob = w)  # sample rows
    X1_star  <- X1[idx]
    X2_star  <- X2[idx]
    
    #AFT linear predictors
    mu0 <- eta0[m] + eta_x1[m] * X1_star + eta_x2[m] * X2_star 
    mu1 <- mu0 + eta_a[m]
    
    #Weibull survival at t0 and average over X*
    b0  <- exp(-k_draws[m] * mu0)
    b1  <- exp(-k_draws[m] * mu1)
    t0k <- t0 ^ k_draws[m]
    
    S0_marg[m] <- mean(exp(-b0 * t0k))
    S1_marg[m] <- mean(exp(-b1 * t0k))
  }
  return(list(S0_marg = S0_marg, S1_marg = S1_marg))
}
## --- single-replicate runner ---
run_one_rep <- function(r,
                        nburn = 3000, M = 5000,
                        N = 1000,  B = 10000, t0 = 2) {
  
  set.seed(1000 + r)  
  
  # --- JAGS model  ---
  # #bias parameter: gamma0-Unif(-5,5) and rest - Unif(-2,2) 
  #remaning parameters: intercept N(0, 5^2), remaning N(0,2^2)
  jags_model <- "
  model {
    for (i in 1:N) {
      is_cens[i] ~ dinterval(T[i], C[i])
      T[i] ~ dweib(k, b[i])

      mu[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] + eta_u*U[i]
      log_b[i] <- -k * mu[i]
      b[i] <- exp(log_b[i])

      U[i] ~ dbern(pU[i])
      logit(pU[i]) <- gamma0 + gamma1*X1[i] + gamma2*X2[i]

      A[i] ~ dbern(pA[i])
      logit(pA[i]) <- alpha0 + alpha1*X1[i] + alpha2*X2[i] + alpha3*U[i] 
    }

    eta0 ~ dnorm(0, 0.1)
    eta_x1 ~ dnorm(0, 0.25)
    eta_x2 ~ dnorm(0, 0.25)
    eta_a ~ dnorm(0, 0.25)
    eta_u ~ dunif(-2, 2)      # sensitivity effect of U on outcome
  
    #U parameters
    gamma0 ~ dunif(-5, 5)
    gamma1 ~ dunif(-2, 2)
    gamma2 ~ dunif(-2, 2)

    #Treatment parameters
    alpha0 ~ dnorm(0, 0.1)
    alpha1 ~ dnorm(0, 0.25)
    alpha2 ~ dnorm(0, 0.25)
    alpha3 ~ dunif(-2, 2)   # sensitivity effect of U on treatment
    
    #logk ~ dnorm(log(2), 4) #increasing hazard
    #k <- exp(logk) 
    k ~ dgamma(4, 2)
  }"
  
  jags_model_obsU <- "
  model {
    for (i in 1:N) {
      is_cens[i] ~ dinterval(T[i], C[i])
      T[i] ~ dweib(k, b[i])

      mu[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] + eta_u*U[i]
      log_b[i] <- -k * mu[i]
      b[i] <- exp(log_b[i])

    }

    eta0 ~ dnorm(0, 0.1)
    eta_x1 ~ dnorm(0, 0.25)
    eta_x2 ~ dnorm(0, 0.25)
    eta_a ~ dnorm(0, 0.25)
    eta_u ~ dunif(-2, 2)      # sensitivity effect of U on outcome
    
    #logk ~ dnorm(log(2), 4) #increasing hazard
    #k <- exp(logk) 
    k ~ dgamma(4, 2)
  }"
  
  jags_model_naive <- "
  model {
    for (i in 1:N) {
      is_cens[i] ~ dinterval(T[i], C[i])
      T[i] ~ dweib(k, b[i])

      mu[i] <- eta0 + eta_x1*X1[i] + eta_x2*X2[i] + eta_a*A[i] 
      log_b[i] <- -k * mu[i]
      b[i] <- exp(log_b[i])
    }

    eta0 ~ dnorm(0, 0.1)
    eta_x1 ~ dnorm(0, 0.25)
    eta_x2 ~ dnorm(0, 0.25)
    eta_a ~ dnorm(0, 0.25)
    
    #logk ~ dnorm(log(2), 4) #increasing hazard
    #k <- exp(logk) 
    k ~ dgamma(4, 2)
  }"
  
  # --- simulate one dataset  ---
  sim <- gen_data_sim(
    n = N, tau = 5.5,
    alpha_0 = 0.1, alpha_x1 = 0.3, alpha_x2 = 0.7, alpha_u = -0.8,
    eta_intercept_a0 = 0.7, eta_intercept_a1 = 0.2, eta_x1 = -0.1, eta_x2 = 0.4, eta_u = -0.8,
    k = 2, pU = 0.5
  )
  data_sim <- sim$data_sim_noU
  data_sim_U <- sim$data_sim_withU
  
  # --- prepare JAGS data  ---
  T_vec <- data_sim$Time
  T_vec[data_sim$d == 0] <- NA_real_
  eps <- 1e-6
  C_vec <- ifelse(data_sim$d == 1L, data_sim$Time + eps, data_sim$Time)
  is_cens <- 1L - data_sim$d
  X1 <- data_sim$X1
  X2 <- data_sim$X2
  U <- data_sim_U$U
  
  jags_data <- list(
    N = N, T = T_vec, C = C_vec, is_cens = is_cens,
    A = data_sim$A, X1 = X1, X2 = X2
  ) #same for both naive and latU
  
  jags_data_obsU <- list(
    N = N, T = T_vec, C = C_vec, is_cens = is_cens,
    A = data_sim$A, X1 = X1, X2 = X2, U=U
  ) 
  
  # --- initial values ----
  # Base init template
  base_inits <- list(
    eta0 = 0, eta_x1 = 0, eta_x2 = 0, eta_a = 0, eta_u = 0,
    k = 1,
    gamma0 = 0, gamma1 = 0, gamma2 = 0,
    .RNG.name = "base::Wichmann-Hill", .RNG.seed = r
  )
  
  # Subsets per model
  inits_latU  <- base_inits  # full set, includes gamma's
  inits_obsU  <- base_inits[names(base_inits) %in% c("eta0","eta_x1","eta_x2","eta_a","eta_u","k",".RNG.name",".RNG.seed")]
  inits_naive <- base_inits[names(base_inits) %in% c("eta0","eta_x1","eta_x2","eta_a","k",".RNG.name",".RNG.seed")]
  
  
  # --- monitor parameters  ---
  params_latU <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","k","gamma0","gamma1","gamma2")
  params_obsU <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","k")
  params_naive <- c("eta0","eta_x1","eta_x2","eta_a","k")
  
  # --- fit JAGS once per replicate  ---
  #sim1 - latentU approach
  sim1 <- jags.model(textConnection(jags_model), data = jags_data, inits = inits_latU, n.chains = 1)
  post1 <- coda.samples(sim1, params_latU, n.iter = nburn + M)
  post_1 <- data.frame(post1[[1]][(nburn + 1):(nburn + M), ])
  rm(post1)
  
  #sim2 - obsU approach
  sim2 <- jags.model(textConnection(jags_model_obsU), data = jags_data_obsU, inits = inits_obsU, n.chains = 1)
  post2 <- coda.samples(sim2, params_obsU, n.iter = nburn + M)
  post_2 <- data.frame(post2[[1]][(nburn + 1):(nburn + M), ])
  rm(post2)
  
  #sim3 - naive approach
  sim3 <- jags.model(textConnection(jags_model_naive), data = jags_data, inits = inits_naive, n.chains = 1)
  post3 <- coda.samples(sim3, params_naive, n.iter = nburn + M)
  post_3 <- data.frame(post3[[1]][(nburn + 1):(nburn + M), ])
  rm(post3)
  
  
  
  #------------compute spce by predefined function---------
  
  surv_latU <- compute_spce_bb_latU(samp=post_1, X1=data_sim$X1, X2=data_sim$X2, t0=t0, rx=B, 
                                       dirichlet_alpha = 1) 
  surv_obsU <- compute_spce_bb_obsU(samp=post_2, X1=data_sim$X1, X2=data_sim$X2, U=data_sim_U$U, t0=t0, rx=B, 
                                       dirichlet_alpha = 1) 
  surv_naive <- compute_spce_bb_naive(samp=post_3, X1=data_sim$X1, X2=data_sim$X2, t0=t0, rx=B, 
                                    dirichlet_alpha = 1) 
  
  SPCE_draws_latU <- surv_latU$S1_marg - surv_latU$S0_marg
  SPCE_sd_latU    <- sd(SPCE_draws_latU)
  SPCE_mean_latU  <- mean(SPCE_draws_latU)
  SPCE_CI_latU    <- quantile(SPCE_draws_latU, c(0.025, 0.975), names = FALSE)
  
  SPCE_draws_obsU <- surv_obsU$S1_marg - surv_obsU$S0_marg
  SPCE_sd_obsU    <- sd(SPCE_draws_obsU)
  SPCE_mean_obsU  <- mean(SPCE_draws_obsU)
  SPCE_CI_obsU   <- quantile(SPCE_draws_obsU, c(0.025, 0.975), names = FALSE)
  
  SPCE_draws_naive <- surv_naive$S1_marg - surv_naive$S0_marg
  SPCE_sd_naive   <- sd(SPCE_draws_naive)
  SPCE_mean_naive  <- mean(SPCE_draws_naive)
  SPCE_CI_naive   <- quantile(SPCE_draws_naive, c(0.025, 0.975), names = FALSE)
  
  list(
    spce_draws_latU = SPCE_draws_latU,
    spce_mean_latU  = SPCE_mean_latU,
    spce_sd_latU    = SPCE_sd_latU,
    spce_ci_latU    = SPCE_CI_latU,
    spce_draws_obsU  = SPCE_draws_obsU ,
    spce_mean_obsU   = SPCE_mean_obsU ,
    spce_sd_obsU    = SPCE_sd_obsU ,
    spce_ci_obsU    = SPCE_CI_obsU,
    spce_draws_naive  = SPCE_draws_naive ,
    spce_mean_naive  = SPCE_mean_naive,
    spce_sd_naive    = SPCE_sd_naive ,
    spce_ci_naive    = SPCE_CI_naive
  )
}

# -------- parallel across replicates r (no inner parallel) --------
R <- 9

n.cores <- max(1, parallel::detectCores() - 1)
cl <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

set.seed(123)  # master seed for %dorng%

results <- foreach(r = 1:R,
                   .packages = c("rjags", "MCMCprecision"),
                   .export   = c("gen_data_sim", "compute_spce_bb_latU", "compute_spce_bb_obsU","run_one_rep")) %dorng% {
                     run_one_rep(r)
                   }

parallel::stopCluster(cl)

t1 <- Sys.time()
runtime_min <- as.numeric(difftime(t1, t0, units = "mins"))
runtime_min
# combine to a data.frame with columns: spce_mean, spce_ci_lo, spce_ci_hi
library(dplyr)
library(tidyr)
true_spce_aft <- -0.1384194

comb_long <- do.call(bind_rows, lapply(results, function(x) {
  if (is.null(x)) return(NULL)
  tibble(
    method = c("latU","obsU","naive"),
    mean   = c(x$spce_mean_latU,  x$spce_mean_obsU,  x$spce_mean_naive),
    sd     = c(x$spce_sd_latU,    x$spce_sd_obsU,    x$spce_sd_naive),
    ci_lo  = c(x$spce_ci_latU[1], x$spce_ci_obsU[1], x$spce_ci_naive[1]),
    ci_hi  = c(x$spce_ci_latU[2], x$spce_ci_obsU[2], x$spce_ci_naive[2])
  )
}))

write.csv(comb_long, "results_long.csv", row.names = FALSE)

#-------------#-------------#-------------#-------------#-------------#-------------#-------------

# Compute coverage, bias, and summary stats all together
summary_comb <- comb_long %>%
  mutate(
    method = factor(method, levels = c("naive", "latU", "obsU")), 
    bias   = mean - true_spce_aft,
    cover  = (true_spce_aft >= ci_lo & true_spce_aft <= ci_hi),
    width  = ci_hi - ci_lo,
    z_abs  = abs(bias / sd)
  ) %>%
  group_by(method) %>%
  summarise(
    coverage   = mean(cover),
    mean_est   = mean(mean),
    mean_bias  = mean(bias),
    mean_sd    = mean(sd),
    mc_sd      = sd(mean),
    mean_width = mean(width),
    mean_abs_z = mean(z_abs),
    .groups = "drop"
  ) %>%
  arrange(method)

summary_comb
#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------

#draft
t0 <- Sys.time()
res1 <- run_one_rep(1)
t1 <- Sys.time()
runtime_min <- as.numeric(difftime(t1, t0, units = "mins"))
runtime_min

