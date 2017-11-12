
data {
  int<lower=0> ncounts;
  int<lower=0> nquadrats;
  int<lower=0> nobservers;
  int<lower=0> nyears;
  int<lower=0> bird_count[ncounts];
  int<lower=0> quad[ncounts];
  vector[ncounts] year;
  int<lower=0> obser[ncounts];
  int<lower=0> firstyr[ncounts];
}

transformed data {
  vector[ncounts] year_cent;
  year_cent = year - mean(year);      // center year
}

parameters {
  // parameters within non-centered priors
  vector[nobservers] z_obs;
  vector[ncounts] z_overdisp;

  // fixed effect parameters
  vector[nquadrats] int_quadrat;      // quadrat-specific intercepts
  vector[nquadrats] beta_quadrat;     // quadrat-specific slopes (with respect to time)
  real eta;                           // first-year observer effect
  
  // hyperparameters
  real<lower=0> sd_obs;
  real<lower=0> sd_overdisp;
}

transformed parameters {
  vector[ncounts] lambda;
  
  // random effect parameters
  vector[nobservers] obs;
  vector[ncounts] overdisp;

  // hierarchical priors (non-centered parameterization)
  obs = z_obs * sd_obs;
  overdisp = z_overdisp * sd_overdisp;

  // mean of log-poisson distribution
  for(i in 1:ncounts) {
    lambda[i] = int_quadrat[quad[i]] +
                beta_quadrat[quad[i]] * year_cent[i] +
                obs[obser[i]] +
                eta * firstyr[i] +
                overdisp[i];
  }
}

model {
  // unit normal priors for non-centered terms
  z_obs ~ normal(0, 1);
  z_overdisp ~ normal(0, 1);
  
  // priors for fixed effects
  int_quadrat ~ normal(0, 30);
  beta_quadrat ~ normal(0, 30);
  eta ~ normal(0, 30);
  
  // priors for hyperparameters
  sd_obs ~ cauchy(0, 20);
  sd_overdisp ~ cauchy(0, 20);
  
  // likelihood
  bird_count ~ poisson_log(lambda);
}

generated quantities {
  vector[nquadrats] beta_geom;
  beta_geom = 100 * (exp(beta_quadrat) - 1);
}

