
data {
  int<lower=0> ncounts;                   // N observations (i.e. bird counts)
  int<lower=0> ncounties;                 // N counties (i.e. intersection of bcr and unique counties)
  int<lower=0> nbcrs;                     // N bird conservation regions
  int<lower=0> nobservers;                // N observers
  int<lower=0> nyears;                    // N years
  int<lower=0> bird_count[ncounts];       // bird count (response variable)
  int<lower=0> county[ncounts];           // county identifier
  int<lower=0> bcr[ncounts];              // bcr identifier
  int<lower=0> observer[ncounts];         // observer identifier
  int<lower=0> firstyr[ncounts];          // binary indicator for whether overserver's first survey
  vector[ncounts] year;                   // year
  int<lower=0> county_by_bcr[ncounties];  // bcr within which each county is located
  real neonic[ncounties];                 // neonic measure for each county
}

transformed data {
  vector[ncounts] year_cent;
  year_cent = year - mean(year);     // center year
}

parameters {
  // parameters within non-centered priors
  vector[nbcrs] z_int_bcr;
  vector[nbcrs] z_beta_bcr;
  vector[ncounties] z_int_county;
  vector[ncounties] z_beta_county;
  vector[nobservers] z_obs;
  vector[ncounts] z_overdisp;
  
  // fixed effect parameters
  real eta;                      // first-year observer effect
  
  // hyperparameters
  vector[nbcrs] gamma;           // neonicotinoid effect within each bcr
  real mu_int_bcr;               // mean of bcr-specific intercepts
  real mu_beta_bcr;              // mean of bcr-specific slopes (with respect to time)
  real<lower=0> sd_int_bcr;      // sd of bcr-specific intercepts
  real<lower=0> sd_beta_bcr;     // sd of bcr-specific slopes (with respect to time)
  real<lower=0> sd_int_county;   // sd of county-specific intercept deviations
  real<lower=0> sd_beta_county;  // sd of county-specific slope deviations
  real<lower=0> sd_obs;          // sd of observer intercepts
  real<lower=0> sd_overdisp;     // sd of overdispersion term
}

transformed parameters {
  vector[ncounts] lambda;                 // mean of log-poisson distribution
  vector[ncounties] dev_beta_county_hat;  // mean of county-specific slope deviations
  
  // random effect parameters
  vector[nbcrs] int_bcr;                  // bcr-specific intercepts
  vector[nbcrs] beta_bcr;                 // bcr-specific slopes
  vector[ncounties] dev_int_county;       // county-specific intercept deviations
  vector[ncounties] dev_beta_county;      // county-specific slope deviations
  vector[nobservers] obs;                 // observer-specific intercepts
  vector[ncounts] overdisp;               // overdispersion intercepts

  // county-level neonicotinoid model
  for (j in 1:ncounties) {
    dev_beta_county_hat[j] = gamma[county_by_bcr[j]] * neonic[j];
  }
  
  // hierarchical priors (non-centered parameterization)
  int_bcr = mu_int_bcr + z_int_bcr * sd_int_bcr;
  beta_bcr = mu_beta_bcr + z_beta_bcr * sd_beta_bcr;
  dev_int_county = z_int_county * sd_int_county;
  dev_beta_county = dev_beta_county_hat + z_beta_county * sd_beta_county;
  obs = z_obs * sd_obs;
  overdisp = z_overdisp * sd_overdisp;

  // mean of log-poisson distribution
  for(i in 1:ncounts) {
    lambda[i] = int_bcr[bcr[i]] + dev_int_county[county[i]] +
                year_cent[i] * (beta_bcr[bcr[i]] + dev_beta_county[county[i]]) +
                obs[observer[i]] +
                eta * firstyr[i] +
                overdisp[i];
  }
}

model {
  // unit normal priors for non-centered terms
  z_int_bcr ~ normal(0, 1);
  z_beta_bcr ~ normal(0, 1);
  z_int_county ~ normal(0, 1);
  z_beta_county ~ normal(0, 1);
  z_obs ~ normal(0, 1);
  z_overdisp ~ normal(0, 1);
  
  // priors for fixed effects
  eta ~ normal(0, 30);
  
  // priors for hyperparameters
  gamma ~ normal(0, 30);
  mu_int_bcr ~ normal(0, 30);
  mu_beta_bcr ~ normal(0, 30);
  sd_int_bcr ~ cauchy(0, 20);
  sd_beta_bcr ~ cauchy(0, 20);
  sd_int_county ~ cauchy(0, 20);
  sd_beta_county ~ cauchy(0, 20);
  sd_obs ~ cauchy(0, 20);
  sd_overdisp ~ cauchy(0, 20);

  // likelihood
  bird_count ~ poisson_log(lambda);
}

