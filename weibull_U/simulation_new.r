library(survival)
library(rjags)
library(LaplacesDemon)
library(coda)


#-------------- function to generate the data for each simulation run --------------
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
  M <- ifelse(A == 1, pmin(tau, C_a1, D_a1),
              pmin(tau, C_a0, D_a0))
  d <- ifelse(A == 1,
              as.numeric(D_a1 <= pmin(tau, C_a1)),
              as.numeric(D_a0 <= pmin(tau, C_a0)))
  
  # ---- Output dataset ----
  data_sim <- data.frame(M = M, d = d, A = A, X1 = X1, X2=X2)
}


#--------------- funtion to compute spce with U-------------
compute_spce_bb_with_U <- function(
    samp,        # posterior draws; must have: eta0, eta_x, eta_a, k, gamma0, gamma1
    X1, 
    X2,          # observed covariate vector
    t0,          # time at which to evaluate survival
    rx,          # size of BB resample X*
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
    w   <- as.numeric(LaplacesDemon::rdirichlet(1, rep(dirichlet_alpha, L)))
    idx <- sample.int(L, size = rx, replace = TRUE, prob = w)  # sample rows
    X1_star  <- X1[idx]
    X2_star  <- X2[idx]
    # w1_m    <- as.numeric(LaplacesDemon::rdirichlet(1, rep(dirichlet_alpha, L)))
    # w2_m    <- as.numeric(LaplacesDemon::rdirichlet(1, rep(dirichlet_alpha, L)))
    # X1_star <- sample(X1, size = rx, replace = TRUE, prob = w1_m)
    # X2_star <- sample(X2, size = rx, replace = TRUE, prob = w2_m)
    
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

#--------------- funtion to compute spce no U-------------
compute_spce_bb_no_U <- function(
    samp,        # posterior draws; must have: eta0, eta_x, eta_a, k, gamma0, gamma1
    X1, 
    X2,          # observed covariate vector
    t0,          # time at which to evaluate survival
    rx,          # size of BB resample X*
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
    w   <- as.numeric(LaplacesDemon::rdirichlet(1, rep(dirichlet_alpha, L)))
    idx <- sample.int(L, size = rx, replace = TRUE, prob = w)  # sample rows
    X1_star  <- X1[idx]
    X2_star  <- X2[idx]
    # w1_m    <- as.numeric(LaplacesDemon::rdirichlet(1, rep(dirichlet_alpha, L)))
    # w2_m    <- as.numeric(LaplacesDemon::rdirichlet(1, rep(dirichlet_alpha, L)))
    # X1_star <- sample(X1, size = rx, replace = TRUE, prob = w1_m)
    # X2_star <- sample(X2, size = rx, replace = TRUE, prob = w2_m)
    
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

#----------------- Simulation -----------------------
#Simulation runs s=1,2....,S
n <- 1000 #sample size 
tau <- 5.5 #maximum follow-up time
t0 <- 2   #predict at specific t0
true_spce_aft <- -0.1384194
S <- 100 #number of simulation runs
rx <- 5000 #number of resample of covariate X

pos_mean_withU <- sd_withU <- bias_withU <- lo_withU <- hi_withU <- numeric(S)
pos_mean_noU <- sd_noU <- bias_noU <- lo_noU <- hi_noU <- numeric(S)

for (s in 1:S) {
  #----------------generate data----------------
  data_sim <- gen_data_sim(n = 1000, tau = 5.5,
                           alpha_0 = 0.1, alpha_x1 = 0.3, alpha_x2 = 0.7, alpha_u = -0.8,
                           eta_intercept_a0 = 0.7, eta_intercept_a1 = 0.2, eta_x1 = -0.1, eta_x2 = 0.4, eta_u = -0.8,
                           k = 2, pU = 0.5)

  #---------- prepare for jags-----------------
  N <- nrow(data_sim)           #number of observations
  T_vec <- data_sim$M           # 'M' is the observed time 
  T_vec[data_sim$d == 0] <- NA_real_    # For subjects with d == 0 (censored), 
  # set event time to NA since the true event time is unknown.
  # Censoring time equals observed time for everyone
  eps <- 1e-6 
  C_vec <- ifelse(data_sim$d == 1L, data_sim$M + eps, data_sim$M)
  
  # 1 = censored, 0 = event
  is_cens <- 1L - data_sim$d
  
  jags_data <- list(
    N = N,
    T = T_vec,
    C = C_vec,
    is_cens = is_cens,
    A = data_sim$A,
    X1 = data_sim$X1,
    X2 = data_sim$X2
  )

  #-------------------JAGS---------------------
  #considering U
  model_withU <- jags.model("wbmodel_U.txt", data=jags_data, n.chains=1)
  params_U <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","k","gamma0","gamma1","gamma2")
  samp_U <- coda.samples(model_withU, variable.names=params_U, n.iter=20000)
  samp_U = data.frame(samp_U[[1]][10001:20000, ])
  
  #not considering U
  model_noU <- jags.model("wbmodel_noU.txt", data=jags_data, n.chains=1)
  params_noU <- c("eta0","eta_x1","eta_x2","eta_a","k")
  samp_noU <- coda.samples(model_noU, variable.names=params_noU, n.iter=20000)
  samp_noU = data.frame(samp_noU[[1]][10001:20000, ])
  
  #considering U indep
  model_withU <- jags.model("wbmodel_U_indep.txt", data=jags_data, n.chains=1)
  params_U <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","k","pU")
  samp_U <- coda.samples(model_withU, variable.names=params_U, n.iter=20000)
  samp_U = data.frame(samp_U[[1]][10001:20000, ])

  #------------compute spce by predefined function---------
  
  surv_withU <- compute_spce_bb_with_U(samp=samp_U, X1=data_sim$X1, X2=data_sim$X2, t0=t0, rx=rx, 
                                        dirichlet_alpha = 1) 
  surv_noU <- compute_spce_bb_no_U(samp=samp_noU, X1=data_sim$X1, X2=data_sim$X2, t0=t0, rx=rx, 
                                       dirichlet_alpha = 1) 
  
  SPCE_draws_withU <- surv_withU$S1_marg-surv_withU$S0_marg
  SPCE_mean_withU  <- mean(surv_withU$S1_marg)-mean(surv_withU$S0_marg)
  SPCE_CI_withU    <- quantile(SPCE_draws_withU, c(0.025, 0.975))
  
  SPCE_draws_noU <- surv_noU$S1_marg-surv_noU$S0_marg
  SPCE_mean_noU  <- mean(surv_noU$S1_marg)-mean(surv_noU$S0_marg)
  SPCE_CI_noU    <- quantile(SPCE_draws_noU, c(0.025, 0.975))
  
  # --- Posterior summaries of marginal SPCE ---
  pos_mean_withU[s] <- SPCE_mean_withU
  sd_withU[s] <- sd(SPCE_draws_withU)
  bias_withU[s] <- SPCE_mean_withU - true_spce_aft
  lo_withU[s] <- SPCE_CI_withU[1]
  hi_withU[s] <- SPCE_CI_withU[2] #95% equal-tail credible interval from draws 
  
  pos_mean_noU[s] <- SPCE_mean_noU
  sd_noU[s] <- sd(SPCE_draws_noU)
  bias_noU[s] <- SPCE_mean_noU - true_spce_aft
  lo_noU[s] <- SPCE_CI_noU[1]
  hi_noU[s] <- SPCE_CI_noU[2] #95% equal-tail credible interval from draws 
}

pos_mean_withU;sd_withU;bias_withU;lo_withU;hi_withU
pos_mean_noU;sd_noU;bias_noU;lo_noU;hi_noU

#----summary---
#withU
avg_pos_mean_withU <- mean(pos_mean_withU)#average posterior mean across S simulation runs
ESE_withU <- sd(pos_mean_withU)#empirical standard error
ASD_withU <- mean(sd_withU)#mean standard deviation across S simulation runs
avg_bias_withU <- mean(bias_withU) #average of bias
CP_withU <- 100 * mean(true_spce_aft >= lo_withU & true_spce_aft <= hi_withU)
within_CI_withU <- ifelse(true_spce_aft >= lo_withU & true_spce_aft <= hi_withU, "Yes", "No")

#noU
avg_pos_mean_noU <- mean(pos_mean_noU)#average posterior mean across S simulation runs
ESE_noU <- sd(pos_mean_noU)#empirical standard error
ASD_noU <- mean(sd_noU)#mean standard deviation across S simulation runs
avg_bias_noU <- mean(bias_noU) #average of bias
CP_noU <- 100 * mean(true_spce_aft >= lo_noU & true_spce_aft <= hi_noU)
within_CI_noU <- ifelse(true_spce_aft >= lo_noU & true_spce_aft <= hi_noU, "Yes", "No")


#withU summary
summary_df_withU <- data.frame(
  avg_pos_mean_withU = avg_pos_mean_withU, 
  ESE_withU = ESE_withU,
  ASD_withU = ASD_withU,
  avg_bias_withU = avg_bias_withU,
  CP_withU = CP_withU
)

#noU summary
summary_df_noU <- data.frame(
  avg_pos_mean_noU = avg_pos_mean_noU, 
  ESE_noU = ESE_noU,
  ASD_noU = ASD_noU,
  avg_bias_noU = avg_bias_noU,
  CP_noU = CP_noU
)

print(summary_df_withU)
print(summary_df_noU)
#true -0.1384194



# Collect results
results <- data.frame(
  sim_withU = 1:length(pos_mean_withU),
  pos_mean_withU = pos_mean_withU,
  sd_withU = sd_withU,
  bias_withU = bias_withU,
  lower_ci_withU = lo_withU, 
  upper_ci_withU = hi_withU,
  within_CI_withU = within_CI_withU,
  pos_mean_noU = pos_mean_noU,
  sd_noU = sd_noU,
  bias_noU = bias_noU,
  lower_ci_noU = lo_noU,
  upper_ci_noU = hi_noU,
  within_CI_noU = within_CI_noU
)

write.csv(results, "spce_results_both_100_WeaklyInfo.csv", row.names = FALSE)
