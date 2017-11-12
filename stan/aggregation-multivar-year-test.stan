data {
  int<lower=0> N;                // N total
  int<lower=0> M;                // N species
  int<lower=0> K;                // N BCR
  int<lower=0> T;                // N time periods
  
  vector[N] est_mean;              // estimated effects
  vector<lower=0>[N] est_sigma;    // s.e. of estimated effects
  
  int<lower=0> M_index[N];       // Index for species
  int<lower=0> K_index[N];       // Index for BCR
  int<lower=0> T_index[N];       // Index for time period
  
  real period[N];                // binary predictor for whether pre- or post-neonic period
}

parameters {
  vector[K] eta_K;
  vector[M] eta_M;
  matrix[K,T] eta_T;
  
  vector[K] eta_beta;
  
  real alpha;
  
  real mu_beta;
  real<lower=0> sigma_beta;
  
  real<lower=0> sigma_theta_K;
  real<lower=0> sigma_theta_M;
  vector<lower=0>[K] sigma_theta_T;
}

transformed parameters {
  vector[K] theta_K;
  vector[M] theta_M;
  matrix[K,T] theta_T;
  
  vector[K] beta;
  vector[N] theta_hat;
  
  beta = mu_beta + sigma_beta * eta_beta;
  theta_K = eta_K * sigma_theta_K;
  theta_M = eta_M * sigma_theta_M;
  
  for(k in 1:K) {
    theta_T[k,] = eta_T[k,] * sigma_theta_T[k];
  }
  
  for(i in 1:N) {
    theta_hat[i] = alpha + beta[K_index[i]]*period[i] + theta_K[K_index[i]] + theta_M[M_index[i]] + theta_T[K_index[i],T_index[i]];
  }
  
}

model {
  eta_K ~ normal(0, 1);
  eta_M ~ normal(0, 1);
  for(k in 1:K) eta_T[k,] ~ normal(0, 1);
  
  alpha ~ normal(0, 30);
  
  sigma_theta_K ~ cauchy(0, 30);
  sigma_theta_M ~ cauchy(0, 30);
  sigma_theta_T ~ cauchy(0, 30);
  
  // hierarchical change model
  eta_beta ~ normal(0, 1);
  mu_beta ~ normal(0, 20);
  sigma_beta ~ cauchy(0, 20);
  
  est_mean ~ normal(theta_hat, est_sigma);
}

generated quantities {
  vector[N] log_lik;
  matrix[K,T] theta_hat_bcr_t;
  
  for(i in 1:N) {
    log_lik[i] = normal_lpdf(est_mean[i] | theta_hat[i], est_sigma[i]);
  }
    
  for(k in 1:K) {
    for(t in 1:(T-1)) {
      theta_hat_bcr_t[k,t] = alpha + theta_K[k] + theta_T[k,t];
    }
    theta_hat_bcr_t[k,T] = alpha + beta[k] + theta_K[k] + theta_T[k,T];
  }
}
