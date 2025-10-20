library(survival)
library(rjags)
library(LaplacesDemon)
library(coda)
n <- 1000 #sample size 
tau <- 5.5 #maximum follow-up time

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

#----------------- Simulation -----------------------
#Simulation runs s=1,2....,S
t0 <- 2   #predict at specific t0
true_spce_aft <- -0.1384194
S <- 5 #number of simulation runs
rx <- 100000 #number of resample of covariate X

pos_mean <- sd <- bias <- lo <- hi<- numeric(S)


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
  params <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","k","gamma0","gamma1","gamma2")
  samp <- coda.samples(model_withU, variable.names=params, n.iter=20000)
  samp = data.frame(samp[[1]][10001:20000, ])
  
  #not considering U
  #model_noU <- jags.model("wbmodel_noU.txt", data=jags_data, n.chains=1)
  #params_noU <- c("eta0","eta_x","eta_a","eta_u","k")
  #samp_noU <- coda.samples(m, variable.names=params, n.iter=20000)
  #samp_noU = data.frame(samp[[1]][10001:20000, ])
  #summary(samp_noU)
  #------------compute spce by predefined function---------
  
  surv_withU <- compute_spce_bb_with_U(samp=samp, X1=data_sim$X1, X2=data_sim$X2, t0=t0, rx=rx, 
                                        dirichlet_alpha = 1) 
  
  SPCE_draws <- surv_withU$S1_marg-surv_withU$S0_marg
  SPCE_mean  <- mean(surv_withU$S1_marg)-mean(surv_withU$S0_marg)
  SPCE_CI    <- quantile(SPCE_draws, c(0.025, 0.975))
  
  # --- Posterior summaries of marginal SPCE ---
  pos_mean[s] <- SPCE_mean
  sd[s] <- sd(SPCE_draws)
  bias[s] <- SPCE_mean - true_spce_aft
  lo[s] <- SPCE_CI[1]
  hi[s] <- SPCE_CI[2] #95% equal-tail credible interval from draws 
  
}



#----summary---
avg_pos_mean <- mean(pos_mean)#avergae posterior mean across S simulation runs
ESE <- sd(pos_mean)#empirical standard error
ASD <- mean(sd)#mean standard deviation across S simulation runs
avg_bias <- mean(bias) #average of bias
CP <- 100 * mean(true_spce_aft >= lo & true_spce_aft <= hi)

summary_df <- data.frame(
  avg_pos_mean = avg_pos_mean,
  ESE = ESE,
  ASD = ASD,
  avg_bias = avg_bias,
  CP = CP
)
print(summary_df)


# Collect results
results <- data.frame(
  sim = 1:length(pos_mean),
  pos_mean = pos_mean,
  sd = sd,
  bias = bias,
  lower_ci = lo,
  upper_ci = hi
)

write.csv(results, "spce_results_U.csv", row.names = FALSE)
