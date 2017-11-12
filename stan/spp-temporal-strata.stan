
data {
  int<lower=0> ncounts;
  int<lower=0> nstrata;
  int<lower=0> nobservers;
  int<lower=0> nyears;
  int<lower=0> bird_count[ncounts];
  int<lower=0> strat[ncounts];
  vector[ncounts] year;
  vector[nyears] year_uniq;
  int<lower=0> obser[ncounts];
  int<lower=0> firstyr[ncounts];
  real<lower=0> area[nstrata];
  real<lower=0> nonzero[nstrata];
  int<lower=0> n_year_strat;
  int<lower=0> year_by_strat[n_year_strat, 2];
  int<lower=0> year_strat_index[ncounts];
}

parameters {
  vector[n_year_strat] yeareffect;
  vector[nobservers] obs;
  real eta;
  vector[ncounts] overdisp;
  real<lower=0> sd_obs;
  real<lower=0> sd_overdisp;
}

transformed parameters {
  vector[ncounts] lambda;

  for(i in 1:ncounts) {
    lambda[i] = yeareffect[year_strat_index[i]] +
                obs[obser[i]] +
                eta * firstyr[i] +
                overdisp[i];
  }
}

model {
  // priors for fixed effects
  eta ~ normal(0, 30);
  yeareffect ~ normal(0, 30);
  
  // priors for hierarhical parameters (centered parameterization)
  obs ~ normal(0, sd_obs);
  overdisp ~ normal(0, sd_overdisp);
  
  // priors for hyperparameters
  sd_obs ~ cauchy(0, 20);
  sd_overdisp ~ cauchy(0, 20);
  
  // likelihood
  bird_count ~ poisson_log(lambda);
}

generated quantities {
  real totarea;
  matrix[nstrata, nyears] n;
  matrix[nstrata, nyears] N;
  matrix[nstrata, nyears] yeareffect_mat;
  vector[nyears] CompIndex;

  totarea = sum(area);

  for(j in 1:nstrata) {
    for(k in 1:nyears) {
      yeareffect_mat[j,k] = 0.0;
    }
  }

  for (l in 1:n_year_strat) {
    int j;
    int k;
    k = year_by_strat[l,1];
    j = year_by_strat[l,2];

    yeareffect_mat[j,k] = yeareffect[l];
  }

  for(j in 1:nstrata) {
    for(k in 1:nyears) {
      n[j,k] = nonzero[j] * exp(yeareffect_mat[j,k] +
                                0.5 * pow(sd_obs, 2) +
                                0.5 * pow(sd_overdisp, 2));

      N[j,k] = area[j] * n[j,k] / totarea;
    }
  }

  for(k in 1:nyears) {
    CompIndex[k] = sum(N[1:nstrata,k]);
  }
}
