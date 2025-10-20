library(survival)
library(rjags)
library(LaplacesDemon)
library(coda)
n <- 1000#sample size 
tau <- 5.5 #maximum follow-up time
true_spce_aft <- -0.1384194
t0 <- 2
rx <- 100000

#----------------generate data----------------
#covariate
X1 <- rbinom(n, size = 1, prob = 0.4)
X2 <- rnorm(n, mean = 0, sd = 1)
#unmeasured confounder
U <- rbinom(n, size = 1, prob = 0.5)

#treatment assignment
alpha_0  <- 0.1; alpha_x1 <- 0.3; alpha_x2 <- 0.7; alpha_u <- -0.8
linpred_A <- alpha_0 + alpha_x1 * X1 + alpha_x2 * X2 + alpha_u * U
pi_A      <- plogis(linpred_A)
A         <- rbinom(n, size = 1, prob = pi_A)

#time
eta_intercept_a0 <- 0.7; eta_intercept_a1 <- 0.2; eta_x1 <- -0.1; eta_x2 <- 0.4; eta_u  <- -0.8
k <- 2
mu1 <- eta_intercept_a1 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
mu0 <- eta_intercept_a0 + eta_x1 * X1 + eta_x2 * X2 + eta_u * U
b1 <- exp(-k * mu1) 
b0 <- exp(-k * mu0)   
U0 <- runif(n)
U1 <- runif(n)
D_a0 <- (-log(U0) / b0)^(1/k)     
D_a1 <- (-log(U1) / b1)^(1/k) 


# ---- Non-informative censoring ----
C_a0 <- runif(n, 0.1, tau)
C_a1 <- runif(n, 0.1, tau)

# ---- Observed time & event indicator ----
M <- ifelse(A == 1, pmin(tau, C_a1, D_a1),
            pmin(tau, C_a0, D_a0))
d<- ifelse(A == 1,
           as.numeric(D_a1 <= pmin(tau, C_a1)),
           as.numeric(D_a0 <= pmin(tau, C_a0)))

## ---- Observed dataset ----
data_sim <- data.frame(
  M = M,
  d = d,
  A = A,
  X1 = X1,
  X2 = X2
)

#mean(D_a1 > t0) - mean(D_a0 > t0)

#---------- prepare for jags-----------------
N <- nrow(data_sim)           #number of observations
# Observed time
T_vec <- data_sim$M
T_vec[data_sim$d == 0] <- NA_real_   # For subjects with d == 0 (censored), 
# set event time to NA since the true event time is unknown.
eps <- 1e-6 
C_vec <- ifelse(data_sim$d == 1L, data_sim$M + eps, data_sim$M)
#C_vec <- ifelse(data_sim$d == 1L, data_sim$M + 1, data_sim$M)
# If the event occurred (d == 1), define censoring time slightly *after* M 
# ensures event time < censoring time in survival analysis in JAGS input.
# If censored (d == 0), the censoring time equals the observed time M.


# 1 = censored, 0 = event
is_cens <- 1L - data_sim$d

jags_data <- list(
  N = nrow(data_sim),
  T = T_vec,
  C = C_vec,
  is_cens = is_cens,
  A = data_sim$A,
  X1 = data_sim$X1,
  X2 = data_sim$X2
)

#-------------------JAGS---------------------
m <- jags.model("wbmodel_U.txt", data=jags_data, n.chains=1)
params <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","k","gamma0","gamma1","gamma2")
samp <- coda.samples(m, variable.names=params, n.iter=20000)
samp = data.frame(samp[[1]][10001:20000, ])
summary(samp)


#m <- jags.model("wbmodel_U.txt", data=jags_data, n.chains=1)
#params <- c("eta0","eta_x","eta_a","eta_u","k","gamma0","gamma1")
#samp <- coda.samples(m, variable.names=params, n.iter=20000)
#samp = data.frame(samp[[1]][10001:20000, ])
#summary(samp)

#-------------- sample new X, X* first ---------------------
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
    w1_m    <- as.numeric(LaplacesDemon::rdirichlet(1, rep(dirichlet_alpha, L)))
    w2_m    <- as.numeric(LaplacesDemon::rdirichlet(1, rep(dirichlet_alpha, L)))
    X1_star <- sample(X1, size = rx, replace = TRUE, prob = w1_m)
    X2_star <- sample(X2, size = rx, replace = TRUE, prob = w2_m)
    
    #U* | X*
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

surv_result <- compute_spce_bb_with_U(samp=samp, X1=data_sim$X1, X2=data_sim$X2, t0=t0, rx=rx, 
                                      dirichlet_alpha = 1) 

SPCE_draws <- surv_result$S1_marg - surv_result$S0_marg 
SPCE_mean  <- mean(SPCE_draws)
SPCE_sd <-sd(SPCE_draws) #sd over M draws
SPCE_CI    <- quantile(SPCE_draws, c(0.025, 0.975))
bias <- SPCE_mean - true_spce_aft 

SPCE_mean;SPCE_sd;SPCE_CI;bias
#--------------------------------------
#true: -0.1384194



#--------------draft-------------
#sampling U from P(U|X) as Umat 
# spce_drawU_bb <- function(m, X, t0,
#                           eta0, eta_x, eta_a, eta_u, k_draws, pU,
#                           LU ) {
#   # m: index of posterior draw
#   N   <- length(X)
#   t0k <- t0^k_draws[m]
#   #pU  <- plogis(gamma0[m] + gamma1[m] * X)       
#   
#   # Monte Carlo over U|X: N x LU
#   Umat <- matrix(rbinom(N*LU, 1, rep(pU[m], each = LU)), nrow = N, ncol = LU)
#   Xmat <- matrix(X, nrow = N, ncol = LU)
#   
#   # Linear predictors (A=0,1)
#   mu0 <- eta0[m] + eta_x[m]*Xmat                 + eta_u[m]*Umat
#   mu1 <- eta0[m] + eta_x[m]*Xmat + eta_a[m]      + eta_u[m]*Umat
#   
#   b0 <- exp(-k_draws[m] * mu0);  b1 <- exp(-k_draws[m] * mu1)
#   S0 <- exp(-b0 * t0k);    S1 <- exp(-b1 * t0k)
#   
#   # Average over LU draws of U for each i
#   SPCE_bar <- rowMeans(S1) - rowMeans(S0)        # length N
#   
#   # One Bayesian-bootstrap draw over X
#   w <- as.numeric(LaplacesDemon::rdirichlet(1, rep(1, N)))
#   sum(w * SPCE_bar)                               # scalar psi for draw m
# }
# 
# 
# #no sampling U, Using P(U|X) directly
# spce_drawU_bb2 <- function(m, X, t0,
#                            eta0, eta_x, eta_a, eta_u, k_draws, pU
# ){
#   
#   N   <- length(X)
#   t0k <- t0^k_draws[m]
#   #pU  <- plogis(gamma0[m] + gamma1[m]*X)
#   
#   mu0_u0 <- eta0[m] + eta_x[m]*X
#   mu0_u1 <- mu0_u0 + eta_u[m]
#   mu1_u0 <- eta0[m] + eta_x[m]*X + eta_a[m]
#   mu1_u1 <- mu1_u0 + eta_u[m]
#   
#   b0_u0 <- exp(-k_draws[m]*mu0_u0); b0_u1 <- exp(-k_draws[m]*mu0_u1)
#   b1_u0 <- exp(-k_draws[m]*mu1_u0); b1_u1 <- exp(-k_draws[m]*mu1_u1)
#   
#   S0_marg <- (1-pU[m])*exp(-b0_u0*t0k) + pU[m]*exp(-b0_u1*t0k)
#   S1_marg <- (1-pU[m])*exp(-b1_u0*t0k) + pU[m]*exp(-b1_u1*t0k)
#   SPCE_marg  <- S1_marg - S0_marg
#   
#   w <- as.numeric(LaplacesDemon::rdirichlet(1, rep(1, length(X))))
#   sum(w *  SPCE_marg)
#   
# }

#set.seed(2025)
# LU <- 500
# psi_list <- lapply(seq_len(Mdraws), function(m) {
#   spce_drawU_bb2(m, X, t0, eta0,eta_x,eta_a,eta_u,k_draws,pU
#                  #, LU
#   )
# })
# psi_all <- unlist(psi_list)

# mean_SPCE <- mean(psi_all)
# sd_SPCE   <- sd(psi_all)
# ci_SPCE   <- quantile(psi_all, c(0.025, 0.975))
# bias <- mean_SPCE - true_spce_aft
# mean_SPCE; sd_SPCE; ci_SPCE; bias