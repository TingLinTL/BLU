library(survival)
library(LaplacesDemon)
n <- 1000 #sample size 
tau <- 5.5 #maximum follow-up time
true_spce_aft <-  -0.2874432
t0 <- 2 #predict time
rx <- 1000000 #number of resample of covariate X
#covariate
X <- rbinom(n, size = 1, prob = 0.4)

#treatment assignment
alpha_0  <- 0.1; alpha_x <- 0.3
linpred_A <- alpha_0 + alpha_x * X
pi_A      <- plogis(linpred_A)
A         <- rbinom(n, size = 1, prob = pi_A)

#time
eta_intercept_a0 <- 0.7; eta_intercept_a1 <- 0.2; eta_x <- -0.1
k <- 2
mu1 <- eta_intercept_a1 + eta_x * X
mu0 <- eta_intercept_a0 + eta_x * X
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

# --------------- JAGS model --------------- 
library(rjags)
library(coda)

N <- nrow(data_sim)
T_vec <- data_sim$M
T_vec[data_sim$d == 0] <- NA_real_
C_vec <- ifelse(data_sim$d == 1L, data_sim$M + 1, data_sim$M)

jags_data <- list(
    N = N,
    T = T_vec,
    C = C_vec,
    is_cens = 1L - data_sim$d,     # 1=censored, 0=event
    A = data_sim$A,
    X = data_sim$X
  )

#c("beta0","beta1","beta2","k","S0","S1","SPCE")

#-------------------------------------------
#m <- jags.model("wbmodel.txt", data=jags_data, n.chains=3, n.adapt=2000)
#update(m, 4000)
#samp <- coda.samples(m, c("beta0","beta1","beta2","k","SPCE"), n.iter=12000, thin=2)
#summary(samp[, c("beta0","beta1","beta2","k")])

#-------------------------------------------
m <- jags.model("wbmodel.txt", data=jags_data, n.chains=1)
params <-c("eta0","eta_x","eta_a","k")
samp <- coda.samples(m, variable.names=params, n.iter=20000)
samp = data.frame(samp[[1]][10001:20000, ])
summary(samp)

#---------------------------resample X new without matrix ---------------------------------------
#Extract posterior matrix
post <- as.matrix(samp)  

# pull parameter vectors
eta0  <- post[,"eta0"]
eta_x <- post[,"eta_x"]
eta_a <- post[,"eta_a"]
k_draws  <- post[,"k"]         
M<- nrow(post) 

#Bayesian bootstrap over X
w <- as.numeric(LaplacesDemon::rdirichlet(1, rep(1, N)))   
#resample new X from reweighted X
X_star <- sample(X, size = rx, replace = TRUE, prob = w)

# --- precompute t0^k for all posterior draws ---
t0k <- t0 ^ k_draws      
S0_marg <- numeric(M)          # running sums over X*
S1_marg<- numeric(M)

           
for (m in 1:M) {
  # linear predictors over the rx X* values
  mu0 <- eta0[m] + eta_x[m] * X_star
  mu1 <- eta0[m] + eta_a[m] + eta_x[m] * X_star
  
  # Weibull scale terms
  b0 <- exp(-k_draws[m] * mu0)
  b1 <- exp(-k_draws[m] * mu1)
  
  # Survival at t0 for each X* then average over X* 
  S0_marg[m] <- mean( exp(- b0 * t0k[m]) )
  S1_marg[m] <- mean( exp(- b1 * t0k[m]) )
}

SPCE_draws <- S1_marg - S0_marg 
SPCE_mean  <- mean(SPCE_draws)
SPCE_sd <-sd(SPCE_draws) #sd over M draws
SPCE_CI    <- quantile(SPCE_draws, c(0.025, 0.975))
bias <- SPCE_mean - true_spce_aft 

SPCE_mean;SPCE_sd;SPCE_CI;bias
#---------------------------resample X new with matrix ---------------------------------------
#Extract posterior matrix
post <- as.matrix(samp)  

# pull parameter vectors
eta0  <- post[,"eta0"]
eta_x <- post[,"eta_x"]
eta_a <- post[,"eta_a"]
k_draws  <- post[,"k"]         
M<- nrow(post) 

# Convert to matrices
eta0_mat <- matrix(eta0, ncol = 1)
eta_x_mat <- matrix(eta_x, ncol = 1)
eta_a_mat <- matrix(eta_a, ncol = 1)
k_draws_mat <- matrix(k_draws, ncol = 1)

#Bayesian bootstrap over X

w <- as.numeric(LaplacesDemon::rdirichlet(1, rep(1, N)))   
#resample new X from reweighted X
X_star <- sample(X, size = rx, replace = TRUE, prob = w)
X_star_mat <- matrix(X_star, nrow = 1)
L <- length(X_star) 


# Compute the linear predictors
mu0_mat <- eta0_mat %*% matrix(1, nrow = 1, ncol = L) + eta_x_mat %*% X_star_mat
mu1_mat <- eta0_mat %*% matrix(1, nrow = 1, ncol = L) + eta_x_mat %*% X_star_mat + eta_a_mat %*% matrix(1, nrow = 1, ncol = L)  

k_grid <- k_draws_mat %*% matrix(1, nrow = 1, ncol = L)
# Convert to rates
log_b0_mat <- -(k_grid * mu0_mat)
log_b1_mat <- -(k_grid * mu1_mat)
b0_mat <- exp(log_b0_mat)
b1_mat <- exp(log_b1_mat)
# Weibull survival with shape k and rate b S(t) = exp( - b * t^k )
S0_mat <- exp(- b0_mat * (t0 ^ k_grid))
S1_mat <- exp(- b1_mat * (t0 ^ k_grid))

#over L for X
S0_marg <- rowMeans(S0_mat)   # length M
S1_marg <- rowMeans(S1_mat)   

# --- Posterior summaries of marginal SPCE ---
#Across all posterior draws m = 1,...,M
SPCE_draws <- S1_marg - S0_marg   # length M
SPCE_mean  <- mean(SPCE_draws)    # posterior mean
SPCE_sd <-sd(SPCE_draws) #sd over M draws
bias <- SPCE_mean - true_spce_aft ; bias
SPCE_CI <- quantile(SPCE_draws, c(0.025, 0.975))  # 95% CI

SPCE_mean;SPCE_sd;bias;SPCE_CI

#easy trick
#---------------------------BB no resample X---------------------------------------
mat <- as.matrix(samp)  

####not complete, need to be updated later#####

#Bayesian bootstrap over X
w <- rdirichlet(M, rep(1, N))    
SPCE_marg <- rowSums(w * SPCE_mat)             #  M * 1 column

# --- Posterior summaries of marginal SPCE ---
pos_mean <- mean(SPCE_marg);pos_mean #posterior mean of maginal SPCE
sd(SPCE_marg) #sd over 10000 draws
bias <- pos_mean - true_spce_aft ; bias
quantile(SPCE_marg, c(.025, .975))

