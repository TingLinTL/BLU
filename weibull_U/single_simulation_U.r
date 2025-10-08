library(survival)
library(rjags)
library(LaplacesDemon)
library(coda)
n <- 1000 #sample size 
tau <- 5.5 #maximum follow-up time
true_spce_aft <- -0.1464996
t0 <- 2

#----------------generate data----------------
#covariate
X <- rbinom(n, size = 1, prob = 0.4)

#unmeasured confounder
U <- rbinom(n, size = 1, prob = 0.5)

#treatment assignment
alpha_0  <- 0.1; alpha_x <- 0.3; alpha_u <- -0.8
linpred_A <- alpha_0 + alpha_x * X + alpha_u * U
pi_A      <- plogis(linpred_A)
A         <- rbinom(n, size = 1, prob = pi_A)

#time
eta_intercept_a0 <- 0.7; eta_intercept_a1 <- 0.2; eta_x <- -0.1; eta_u  <- -0.8
k <- 2
mu1 <- eta_intercept_a1 + eta_x * X + eta_u * U
mu0 <- eta_intercept_a0 + eta_x * X + eta_u * U
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
  X = X
)


#---------- prepare for jags-----------------
N <- nrow(data_sim)           #number of observations
T_vec <- data_sim$M           # 'M' is the observed time 
T_vec[data_sim$d == 0] <- NA_real_    # For subjects with d == 0 (censored), 
# set event time to NA since the true event time is unknown.
C_vec <- ifelse(data_sim$d == 1L, data_sim$M + 1, data_sim$M)
# If the event occurred (d == 1), define censoring time slightly *after* M 
# ensures event time < censoring time in survival analysis in JAGS input.
# If censored (d == 0), the censoring time equals the observed time M.

jags_data <- list(
  N = N,
  T = T_vec,
  C = C_vec,
  is_cens = 1L - data_sim$d,     # 1=censored, 0=event
  A = data_sim$A,
  X = data_sim$X
)

#-------------------JAGS---------------------
m <- jags.model("wbmodel_U_indep.txt", data=jags_data, n.chains=1)
params <- c("eta0","eta_x","eta_a","eta_u","k","pU")
samp <- coda.samples(m, variable.names=params, n.iter=20000)
samp = data.frame(samp[[1]][10001:20000, ])
summary(samp)


#m <- jags.model("wbmodel_U.txt", data=jags_data, n.chains=1)
#params <- c("eta0","eta_x","eta_a","eta_u","k","gamma0","gamma1")
#samp <- coda.samples(m, variable.names=params, n.iter=20000)
#samp = data.frame(samp[[1]][10001:20000, ])
#summary(samp)
#-------------------------
#Extract posterior matrix
post <- as.matrix(samp)  

# pull parameter vectors
eta0  <- post[,"eta0"]
eta_x <- post[,"eta_x"]
eta_a <- post[,"eta_a"]
eta_u  <- post[,"eta_u"]
k_draws  <- post[,"k"]         
#gamma0 <- post[,"gamma0"]
#gamma1 <- post[,"gamma1"]
pU <- post[,"pU"]
Mdraws <- nrow(post) 


#-----------------
#sampling U from P(U|X) as Umat 
spce_drawU_bb <- function(m, X, t0,
                          eta0, eta_x, eta_a, eta_u, k_draws, pU,
                          LU ) {
  # m: index of posterior draw
  N   <- length(X)
  t0k <- t0^k_draws[m]
  #pU  <- plogis(gamma0[m] + gamma1[m] * X)       
  
  # Monte Carlo over U|X: N x LU
  Umat <- matrix(rbinom(N*LU, 1, rep(pU[m], each = LU)), nrow = N, ncol = LU)
  Xmat <- matrix(X, nrow = N, ncol = LU)
  
  # Linear predictors (A=0,1)
  mu0 <- eta0[m] + eta_x[m]*Xmat                 + eta_u[m]*Umat
  mu1 <- eta0[m] + eta_x[m]*Xmat + eta_a[m]      + eta_u[m]*Umat
  
  b0 <- exp(-k_draws[m] * mu0);  b1 <- exp(-k_draws[m] * mu1)
  S0 <- exp(-b0 * t0k);    S1 <- exp(-b1 * t0k)
  
  # Average over LU draws of U for each i
  SPCE_bar <- rowMeans(S1) - rowMeans(S0)        # length N
  
  # One Bayesian-bootstrap draw over X
  w <- as.numeric(LaplacesDemon::rdirichlet(1, rep(1, N)))
  sum(w * SPCE_bar)                               # scalar psi for draw m
}


#no sampling U, Using P(U|X) directly
spce_drawU_bb2 <- function(m, X, t0,
                           eta0, eta_x, eta_a, eta_u, k_draws, pU
                            ){
  
  N   <- length(X)
  t0k <- t0^k_draws[m]
  #pU  <- plogis(gamma0[m] + gamma1[m]*X)
  
  mu0_u0 <- eta0[m] + eta_x[m]*X
  mu0_u1 <- mu0_u0 + eta_u[m]
  mu1_u0 <- eta0[m] + eta_x[m]*X + eta_a[m]
  mu1_u1 <- mu1_u0 + eta_u[m]
  
  b0_u0 <- exp(-k_draws[m]*mu0_u0); b0_u1 <- exp(-k_draws[m]*mu0_u1)
  b1_u0 <- exp(-k_draws[m]*mu1_u0); b1_u1 <- exp(-k_draws[m]*mu1_u1)
  
  S0_marg <- (1-pU[m])*exp(-b0_u0*t0k) + pU[m]*exp(-b0_u1*t0k)
  S1_marg <- (1-pU[m])*exp(-b1_u0*t0k) + pU[m]*exp(-b1_u1*t0k)
  SPCE_marg  <- S1_marg - S0_marg
  
  w <- as.numeric(LaplacesDemon::rdirichlet(1, rep(1, length(X))))
  sum(w *  SPCE_marg)
  
}
#--------------------------------------

#set.seed(2025)
LU <- 500
psi_list <- lapply(seq_len(Mdraws), function(m) {
  spce_drawU_bb2(m, X, t0, eta0,eta_x,eta_a,eta_u,k_draws,pU
                #, LU
               )
})
psi_all <- unlist(psi_list)
#psi_all
# Summaries
mean_SPCE <- mean(psi_all)
sd_SPCE   <- sd(psi_all)
ci_SPCE   <- quantile(psi_all, c(0.025, 0.975))
bias <- mean_SPCE - true_spce_aft
mean_SPCE; sd_SPCE; ci_SPCE; bias
#true: -0.1464996

#single run
#bb1---------
#mean_SPCE; sd_SPCE; ci_SPCE; bias
# -0.1320271
#0.03266953
#2.5%       97.5% 
#  -0.20128426 -0.07334746 
#0.01447249

#bb2--------
#mean_SPCE; sd_SPCE; ci_SPCE; bias
# -0.1320206
# 0.0326568
#2.5%       97.5% 
#  -0.20119914 -0.07325163 
# 0.01447905
#

