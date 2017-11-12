data {
  int<lower=0> N;                 // number of species
  real est_mean[N];               // estimated treatment effects
  real<lower=0> est_sigma[N];     // s.e. of effect estimates 
}

parameters {
  vector[N] eta;
  real mu_theta;
  real<lower=0> sigma_theta;
}

transformed parameters {
  vector[N] theta;
  theta = mu_theta + sigma_theta * eta;
}

model {
  eta ~ normal(0, 1);
  mu_theta ~ normal(0, 20);
  sigma_theta ~ cauchy(0, 20); 
  
  est_mean ~ normal(theta, est_sigma);
}
