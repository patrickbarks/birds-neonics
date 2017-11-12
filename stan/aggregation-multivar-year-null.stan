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
}

parameters {
  vector[M] eta_M;
  vector[K] eta_K;
  matrix[K,T] eta_T;
  
  real alpha;
  
  real<lower=0> sigma_theta_M;
  real<lower=0> sigma_theta_K;
  vector<lower=0>[K] sigma_theta_T;
}

transformed parameters {
  vector[M] theta_M;
  vector[K] theta_K;
  matrix[K,T] theta_T;
  
  vector[N] theta_hat;
  
  theta_M = eta_M * sigma_theta_M;
  theta_K = eta_K * sigma_theta_K;
  
  for(k in 1:K) {
    theta_T[k,] = eta_T[k,] * sigma_theta_T[k];
  }
  
  for(i in 1:N) {
    theta_hat[i] = alpha + theta_K[K_index[i]] + theta_M[M_index[i]] + theta_T[K_index[i],T_index[i]];
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
  
  est_mean ~ normal(theta_hat, est_sigma);
}

generated quantities {
  vector[N] log_lik;
  matrix[K,T] theta_hat_bcr_t;
  
  for(i in 1:N) {
    log_lik[i] = normal_lpdf(est_mean[i] | theta_hat[i], est_sigma[i]);
  }
  
  for(k in 1:K) {
    for(t in 1:T) {
      theta_hat_bcr_t[k,t] = alpha + theta_K[k] + theta_T[k,t];
    }
  }
}
