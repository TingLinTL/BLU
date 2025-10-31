run_one_rep <- function(r,
                        nburn = 10000, M = 10000,
                        N = 1000, B = 10000, t0 = 2) {
  
  # --- libs & sources ---
  library(doParallel)
  library(MCMCprecision)
  library(rjags)
  source("data_gen_func.R")
  
  set.seed(1000 + r)  # replicate-level RNG (data gen)
  
  # --- JAGS model  ---
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

    eta0 ~ dnorm(0, 0.25)
    eta_x1 ~ dnorm(0, 0.25)
    eta_x2 ~ dnorm(0, 0.25)
    eta_a ~ dnorm(0, 0.25)
    eta_u ~ dnorm(0, 0.25)

    gamma0 ~ dnorm(0, 0.25)
    gamma1 ~ dnorm(0, 0.25)
    gamma2 ~ dnorm(0, 0.25)

    alpha0 ~ dnorm(0, 0.25)
    alpha1 ~ dnorm(0, 0.25)
    alpha2 ~ dnorm(0, 0.25)
    alpha3 ~ dnorm(0, 0.25)

    logk ~ dnorm(log(2), 4)
    k <- exp(logk)
  }"
  
  # --- Generate one replicate of data ---
  data_sim <- gen_data_sim(
    n = N, tau = 5.5,
    alpha_0 = 0.1, alpha_x1 = 0.3, alpha_x2 = 0.7, alpha_u = -0.8,
    eta_intercept_a0 = 0.7, eta_intercept_a1 = 0.2, eta_x1 = -0.1, eta_x2 = 0.4, eta_u = -0.8,
    k = 2, pU = 0.5
  )
  
  T_vec <- data_sim$Time
  T_vec[data_sim$d == 0] <- NA_real_
  eps <- 1e-6
  C_vec <- ifelse(data_sim$d == 1L, data_sim$Time + eps, data_sim$Time)
  is_cens <- 1L - data_sim$d
  X1 <- data_sim$X1
  X2 <- data_sim$X2
  
  jags_data <- list(
    N = N, T = T_vec, C = C_vec, is_cens = is_cens,
    A = data_sim$A, X1 = X1, X2 = X2
  )
  
  inits <- list(
    eta0 = 0, eta_x1 = 0, eta_x2 = 0, eta_a = 0, eta_u = 0,
    logk = log(1.5),
    gamma0 = 0, gamma1 = 0, gamma2 = 0,
    .RNG.name = "base::Wichmann-Hill", .RNG.seed = r
  )
  
  params <- c("eta0","eta_x1","eta_x2","eta_a","eta_u","k","gamma0","gamma1","gamma2")
  
  # --- Fit JAGS once per replicate ---
  sim <- jags.model(textConnection(jags_model), data = jags_data, inits = inits, n.chains = 1)
  post <- coda.samples(sim, params, n.iter = nburn + M)
  post1 <- data.frame(post[[1]][(nburn + 1):(nburn + M), ])
  rm(post)
  
  # --- Define inner work for posterior draw m  ---
  func <- function(m) {
    set.seed(m)
    pi <- as.numeric(MCMCprecision::rdirichlet(n = 1, a = rep(1, N)))
    idx <- sample.int(N, size = B, replace = TRUE, prob = pi)
    X1_star <- X1[idx]; X2_star <- X2[idx]
    pU      <- plogis(post1$gamma0[m] + post1$gamma1[m]*X1_star + post1$gamma2[m]*X2_star)
    U_star  <- rbinom(B, 1, pU)
    
    mu0 <- post1$eta0[m] + post1$eta_x1[m]*X1_star + post1$eta_x2[m]*X2_star + post1$eta_u[m]*U_star
    mu1 <- mu0 + post1$eta_a[m]
    k_m <- post1$k[m]
    b0  <- exp(-k_m * mu0)
    b1  <- exp(-k_m * mu1)
    t0k <- t0^k_m
    
    S0_marg <- mean(exp(-b0 * t0k))
    S1_marg <- mean(exp(-b1 * t0k))
    S1_marg - S0_marg
  }
  
  # --- Parallel over m (inside one replicate) ---
  n.cores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(n.cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  
  Results <- foreach(
    m = 1:M, .combine = 'c',
    .packages = c("MCMCprecision"),
    .export   = c("func","post1","X1","X2","N","B","t0")
  ) %dopar% {
    func(m)
  }
  
  parallel::stopCluster(cl)
  Results
}
