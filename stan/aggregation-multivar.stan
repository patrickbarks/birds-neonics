data {
  int<lower=0> N;                // N total
  int<lower=0> J;                // N groups in factor 1
  int<lower=0> K;                // N group in factor 2
  
  vector[N] est_mean;              // estimated effects
  vector<lower=0>[N] est_sigma;    // s.e. of estimated effects
  
  int<lower=0> J_index[N];       // Index for factor 1
  int<lower=0> K_index[N];       // Index for factor 2
}

parameters {
  vector[J] eta_J;
  vector[K] eta_K;
  
  real alpha;
  
  real<lower=0> sigma_theta_J;
  real<lower=0> sigma_theta_K;
}

transformed parameters {
  vector[J] theta_J;
  vector[K] theta_K;
  
  vector[N] theta_hat;
  
  theta_J = eta_J * sigma_theta_J;
  theta_K = eta_K * sigma_theta_K;
  
  for(i in 1:N) {
    theta_hat[i] = alpha + theta_J[J_index[i]] + theta_K[K_index[i]];
  }
  
}

model {
  eta_J ~ normal(0, 1);
  eta_K ~ normal(0, 1);
  
  alpha ~ normal(0, 30);
  
  sigma_theta_J ~ cauchy(0, 30);
  sigma_theta_K ~ cauchy(0, 30);
  
  est_mean ~ normal(theta_hat, est_sigma);
}

generated quantities {
  vector[N] log_lik;
  
  for(i in 1:N) {
    log_lik[i] = normal_lpdf(est_mean[i] | theta_hat[i], est_sigma[i]);
  }
}
