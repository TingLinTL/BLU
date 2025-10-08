library(survival)
library(LaplacesDemon)

#tau - maximum follow-up time
#alpha_0 & alpha_x, parameters for treatment model
#outcome parameter: eta_intercept_a0 for A=0, eta_intercept_a1 for A=1, eta_x for X
#k is weibull shape

#function to generate the data for each simulation run
gen_data_sim <- function(n, tau, alpha_0, alpha_x,
                              eta_intercept_a0, eta_intercept_a1, eta_x,
                              k, p_X) {
  # ---- Covariate ----
  X <- rbinom(n, size = 1, prob = p_X) 
  
  # ---- Treatment assignment ----
  linpred_A <- alpha_0 + alpha_x * X
  pi_A      <- plogis(linpred_A)
  A         <- rbinom(n, size = 1, prob = pi_A)
  
  # ---- Event times (Weibull) ----
  #generate via inverse CDF
  #T~weibull(k, b), F(t)=1-exp(- b *t^k), where b = exp(-k* lp)
  #T=(-log(U)/b)^(1/k)
  mu1 <- eta_intercept_a1 + eta_x * X
  mu0 <- eta_intercept_a0 + eta_x * X
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
  data_sim <- data.frame(M = M, d = d, A = A, X = X)
}

#Simulation runs s=1,2....,S
t0 <- 2   #predict at specific t0
true_spce <- -0.2874432 #true spce 
S <- 10

pos_mean <- sd <- bias <- lo <- hi<- numeric(S)

for (s in 1:S) {
  data_sim <- gen_data_sim(
    n = 1000, tau = 5.5,
    alpha_0 = 0.1, alpha_x = 0.3,
    eta_intercept_a0 = 0.7, eta_intercept_a1 = 0.2, eta_x = -0.1,
    k = 2, p_X = 0.4)
  
  N <- nrow(data_sim)           #number of observations
  T_vec <- data_sim$M           # 'M' is the observed time 
  T_vec[data_sim$d == 0] <- NA_real_    # For subjects with d == 0 (censored), 
                                         # set event time to NA since the true event time is unknown.
  C_vec <- ifelse(data_sim$d == 1L, data_sim$M + 1, data_sim$M)
  # If the event occurred (d == 1), define censoring time slightly *after* M (e.g., M + 1)
  # ensures event time < censoring time in survival analysis in JAGS input.
  # If censored (d == 0), the censoring time equals the observed time M.
  
  jags_data <- list(
    N = N,
    T = T_vec,
    C = C_vec,
    is_cens = 1L - data_sim$d,     # 1=censored, 0=event
    A = data_sim$A,
    X = data_sim$X,
    t0 = t0      
  )
  
  #-------------------JAGS---------------------
  m <- jags.model("wbmodel.txt", data=jags_data, n.chains=1)
  params <- c("S0","S1","SPCE")
  #params <- c("SPCE")
  samp <- coda.samples(m, variable.names=params, n.iter=20000)
  samp = data.frame(samp[[1]][10001:20000, ])
  
  #------------prepare matrix for next step BB---------
  
  mat <- as.matrix(samp)  
  spce_cols <- grep("^SPCE(\\[|\\.)", colnames(mat), value = TRUE)
  idx       <- as.integer(sub(".*?(\\d+).*", "\\1", spce_cols))  
  ord       <- order(idx)
  SPCE_mat  <- mat[ , spce_cols[ord], drop = FALSE]  
  # optional, check to see if dimention is correct, dim(SPCE_mat)
  Mdraws <- nrow(SPCE_mat) #number of total draws 
  
  
  #------------Bayesian bootstrap over X------------
  w <- rdirichlet(Mdraws, rep(1, N))    #Direchelt weights for each subject and each draw, matrix M rows * N columns
  SPCE_marg <- rowSums(w * SPCE_mat)             #  M rows * 1 column
  
  
  # --- Posterior summaries of marginal SPCE ---
  pos_mean[s] <- mean(SPCE_marg)      #posterior mean of maginal SPCE, store each s
  sd[s] <- sd(SPCE_marg) #sd over draws for each s
  bias[s] <- pos_mean[s] - true_spce #bias for each s
  
  # --- 95% equal-tail credible interval from draws ---
  q <- quantile(SPCE_marg, c(0.025, 0.975), names = FALSE)
  lo[s] <- q[1]; hi[s] <- q[2]
}

avg_pos_mean <- mean(pos_mean); avg_pos_mean#avergae posterior mean across S simulation runs
ESE <- sd(pos_mean);ESE #empirical standard error
ASD <- mean(sd);ASD #mean standard deviation across S simulation runs
avg_bias <- mean(bias);avg_bias #average of bias
#coverage probability
#CP <- 100.00 * mean((true_spce > (pos_mean + qnorm(0.025) * sd)) & (true_spce < (pos_mean + qnorm(0.975) * sd)))
#(true_spce > (pos_mean + qnorm(0.025) * sd)) & (true_spce < (pos_mean + qnorm(0.975) * sd))

CP <- 100 * mean(true_spce >= lo & true_spce <= hi)

CP


# Collect results
results <- data.frame(
  sim = 1:length(pos_mean),
  pos_mean = pos_mean,
  sd = sd,
  bias = bias,
  lower_ci = lo,
  upper_ci = hi
)

results

write.csv(results, "spce_results.csv", row.names = FALSE)
