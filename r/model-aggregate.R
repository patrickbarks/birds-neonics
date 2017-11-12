

# load libraries
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(rstan)
library(loo)

# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# setwd
setwd('~/desktop/bbs/')


# read species and bcr data
species_df <- read_csv('data/species-list.csv')
bcr_dat <- read_csv('data/bbs_raw_2015/BCR.csv') %>% mutate(bcr = as.character(BCR))


# read posterior summaries from species-level models
summary_spp_neonic <-          read_csv('analysis/post-summary-spp-spatial-neonic.csv')
summary_spp_neonic_mismatch <- read_csv('analysis/post-summary-spp-spatial-neonic-mismatch.csv')
summary_spp_imidacloprid <-    read_csv('analysis/post-summary-spp-spatial-imidacloprid.csv')
summary_spp_cropland <-        read_csv('analysis/post-summary-spp-spatial-cropland.csv')





################# Spatial models

#### Functionalize
AggregationModel <- function(data) {
  
  # arrange data in list
  stan_dat <- list(
    N = nrow(data),
    J = length(unique(data$bcr)),
    K = length(unique(data$aou)),
    est_mean = data$gamma_mean,
    est_sigma = data$gamma_se,
    J_index = as.numeric(as.factor(data$bcr)),
    K_index = as.numeric(as.factor(data$aou))
  )
  
  # model pars
  stan_pars <- c(
    'alpha', 'theta_J', 'theta_K',
    'sigma_theta_J', 'sigma_theta_K', 'log_lik'
  )
  
  # stan fit
  fit <- stan(
    file = 'stan/aggregation-multivar.stan',
    data = stan_dat, 
    pars = stan_pars,
    warmup = 1500,
    iter = 4000,
    thin = 2,
    chains = 2,
    control = list(adapt_delta = 0.99, stepsize = 0.01)
  )
  
  # extract posterior samples
  alpha_df <- tibble(alpha = rstan::extract(fit, 'alpha')$alpha)
  theta_int_spp <- rstan::extract(fit, 'theta_K')$theta_K + as.numeric(alpha_df$alpha)
  theta_int_bcr <- rstan::extract(fit, 'theta_J')$theta_J + as.numeric(alpha_df$alpha)
  
  # convert to geometric
  alpha_df$alpha_geom <- 100 * (exp(alpha_df$alpha) - 1)
  theta_int_bcr_geom <- 100 * (exp(theta_int_bcr) - 1)
  theta_int_spp_geom <- 100 * (exp(theta_int_spp) - 1)

  # arrange in dataframes by species and by bcr
  df_spp <- tibble(
    aou = sort(unique(data$aou)),
    theta_mean =  apply(theta_int_spp_geom, 2, mean),
    theta_se =    apply(theta_int_spp_geom, 2, sd),
    theta_med =   apply(theta_int_spp_geom, 2, function(x) quantile(x, 0.500)),
    theta_low90 = apply(theta_int_spp_geom, 2, function(x) quantile(x, 0.050)),
    theta_upp90 = apply(theta_int_spp_geom, 2, function(x) quantile(x, 0.950)),
    theta_low99 = apply(theta_int_spp_geom, 2, function(x) quantile(x, 0.005)),
    theta_upp99 = apply(theta_int_spp_geom, 2, function(x) quantile(x, 0.995)),
    theta_raw_med =   apply(theta_int_spp, 2, function(x) quantile(x, 0.500)),
    theta_raw_low90 = apply(theta_int_spp, 2, function(x) quantile(x, 0.050)),
    theta_raw_upp90 = apply(theta_int_spp, 2, function(x) quantile(x, 0.950)),
    theta_raw_low99 = apply(theta_int_spp, 2, function(x) quantile(x, 0.005)),
    theta_raw_upp99 = apply(theta_int_spp, 2, function(x) quantile(x, 0.995)))

  df_bcr <- tibble(
    bcr = as.character(sort(unique(data$bcr))),
    theta_mean =  apply(theta_int_bcr_geom, 2, mean),
    theta_se =    apply(theta_int_bcr_geom, 2, sd),
    theta_med =   apply(theta_int_bcr_geom, 2, function(x) quantile(x, 0.500)),
    theta_low90 = apply(theta_int_bcr_geom, 2, function(x) quantile(x, 0.050)),
    theta_upp90 = apply(theta_int_bcr_geom, 2, function(x) quantile(x, 0.950)),
    theta_low99 = apply(theta_int_bcr_geom, 2, function(x) quantile(x, 0.005)),
    theta_upp99 = apply(theta_int_bcr_geom, 2, function(x) quantile(x, 0.995)),
    theta_raw_med =   apply(theta_int_bcr, 2, function(x) quantile(x, 0.500)),
    theta_raw_low90 = apply(theta_int_bcr, 2, function(x) quantile(x, 0.050)),
    theta_raw_upp90 = apply(theta_int_bcr, 2, function(x) quantile(x, 0.950)),
    theta_raw_low99 = apply(theta_int_bcr, 2, function(x) quantile(x, 0.005)),
    theta_raw_upp99 = apply(theta_int_bcr, 2, function(x) quantile(x, 0.995)))

  return(tibble(alpha_df = lst(alpha_df), df_spp = lst(df_spp), df_bcr = lst(df_bcr)))
}


summary_agg_neonic <- summary_spp_neonic %>%
  group_by(year) %>%
  do(AggregationModel(.)) %>%
  ungroup()

summary_agg_neonic_mismatch <- summary_spp_neonic_mismatch %>%
  group_by(year) %>%
  do(AggregationModel(.)) %>%
  ungroup()

summary_agg_imidacloprid <- summary_spp_imidacloprid %>%
  group_by(year) %>%
  do(AggregationModel(.)) %>%
  ungroup()

summary_agg_cropland <- summary_spp_cropland %>%
  group_by(year) %>%
  do(AggregationModel(.)) %>%
  ungroup()


# write to file
save(summary_agg_neonic,          file = 'analysis/post-summary-agg-spatial-neonic.RData')
save(summary_agg_neonic_mismatch, file = 'analysis/post-summary-agg-spatial-neonic-mismatch.RData')
save(summary_agg_imidacloprid,    file = 'analysis/post-summary-agg-spatial-imidacloprid.RData')
save(summary_agg_cropland,        file = 'analysis/post-summary-agg-spatial-cropland.RData')





################# Temporal model (strata-level)

# posterior summaries from species-level models
summary_spp_trend_7years <- read_csv('analysis/post-summary-spp-temporal-strata-trend-7years.csv')
summary_spp_trend_allyears <- read_csv('analysis/post-summary-spp-temporal-strata-trend-allyears.csv')


### Get aggregate population trend
AggregateModelTrend <- function(mean, sigma, return_agg = F, return_species = F) {
  
  # arrange data in list
  stan_dat <- list(
    N = length(mean),
    est_mean = mean,
    est_sigma = sigma
  )
  
  # stan fit
  fit <- stan(
    file = 'stan/aggregation-univar.stan',
    data = stan_dat,
    warmup = 2000,
    iter = 4000,
    thin = 2,
    chains = 2,
    control = list(adapt_delta = 0.99, stepsize = 0.01)
  )
  
  # extract posterior samples
  df_mu_theta <- tibble(mu_theta = rstan::extract(fit, 'mu_theta')$mu_theta)
  theta <- rstan::extract(fit, 'theta')$theta + df_mu_theta$mu_theta
  
  # arrange in dataframes by species and by bcr
  df_out <- tibble(mu_theta_med =   quantile(df_mu_theta$mu_theta, 0.500),
                   mu_theta_low90 = quantile(df_mu_theta$mu_theta, 0.050),
                   mu_theta_upp90 = quantile(df_mu_theta$mu_theta, 0.950),
                   mu_theta_low99 = quantile(df_mu_theta$mu_theta, 0.005),
                   mu_theta_upp99 = quantile(df_mu_theta$mu_theta, 0.995))
  
  df_spp <- tibble(aou = sort(unique(summary_spp_trend_allyears$aou)),
                   theta_med =   apply(theta, 2, function(x) quantile(x, 0.500)),
                   theta_low90 = apply(theta, 2, function(x) quantile(x, 0.050)),
                   theta_upp90 = apply(theta, 2, function(x) quantile(x, 0.950)),
                   theta_low99 = apply(theta, 2, function(x) quantile(x, 0.005)),
                   theta_upp99 = apply(theta, 2, function(x) quantile(x, 0.995)))
  
  if(return_agg == TRUE)     df_out$df_mu_theta <- lst(df_mu_theta)
  if(return_species == TRUE) df_out$df_spp <- lst(df_spp)
  
  return(df_out)
}

# get posterior summaries
summary_agg_trend_allyears <- summary_spp_trend_allyears %>%
  do(AggregateModelTrend(mean = .$beta_mean, sigma = .$beta_se,
                         return_agg = T, return_species = T)) %>% ungroup()

summary_agg_trend_7years <- summary_spp_trend_7years %>%
  group_by(year_start) %>%
  do(AggregateModelTrend(mean = .$beta_mean, sigma = .$beta_se,
                         return_agg = F, return_species = F)) %>% ungroup()

# write to file
save(summary_agg_trend_allyears, file = 'analysis/post-summary-agg-temporal-strata-trend-allyears.RData')
save(summary_agg_trend_7years,   file = 'analysis/post-summary-agg-temporal-strata-trend-7years.RData')





################# Temporal model (quadrat-level)

# posterior summaries from species-level models
summary_temporal_quadrat <- read_csv('analysis/post-summary-spp-temporal-quadrat.csv')


#### Functionalize
AggregationModel <- function(data) {
  
  # arrange data in list
  stan_dat <- list(
    N = nrow(data),
    J = length(unique(data$quadrat)),
    K = length(unique(data$aou)),
    est_mean = data$beta_mean,
    est_sigma = data$beta_se,
    J_index = as.numeric(as.factor(data$quadrat)),
    K_index = as.numeric(as.factor(data$aou))
  )
  
  # model pars
  stan_pars <- c(
    'alpha', 'theta_J', 'theta_K',
    'sigma_theta_J', 'sigma_theta_K', 'log_lik'
  )
  
  # stan fit
  fit <- stan(
    file = 'stan/aggregation-multivar.stan',
    data = stan_dat, 
    pars = stan_pars,
    warmup = 1500,
    iter = 4000,
    thin = 2,
    chains = 2,
    control = list(adapt_delta = 0.99, stepsize = 0.01)
  )
  
  # extract posteriors
  alpha_df <- tibble(alpha = rstan::extract(fit, 'alpha')$alpha)
  theta_int_spp <- rstan::extract(fit, 'theta_K')$theta_K + as.numeric(alpha_df$alpha)
  theta_int_quad <- rstan::extract(fit, 'theta_J')$theta_J + as.numeric(alpha_df$alpha)
  
  # convert to geometric
  theta_int_quad_geom <- 100 * (exp(theta_int_quad) - 1)
  
  # arrange in dataframes by quadrat
  df_quad <- tibble(
    quadrat = as.character(sort(unique(data$quadrat))),
    theta_med =   apply(theta_int_quad_geom, 2, function(x) quantile(x, 0.500)),
    theta_low90 = apply(theta_int_quad_geom, 2, function(x) quantile(x, 0.050)),
    theta_upp90 = apply(theta_int_quad_geom, 2, function(x) quantile(x, 0.950)),
    theta_low99 = apply(theta_int_quad_geom, 2, function(x) quantile(x, 0.005)),
    theta_upp99 = apply(theta_int_quad_geom, 2, function(x) quantile(x, 0.995)))
  
  return(df_quad)
}

# get posterior summaries
summary_agg_temporal_quadrat <- summary_temporal_quadrat %>%
  group_by(year) %>%
  do(AggregationModel(.)) %>%
  ungroup()

# write to file
save(summary_agg_temporal_quadrat, file = 'analysis/post-summary-agg-temporal-quadrat.RData')





################# Spatio-temporal analyses (randomization)

# load relevant data
load('analysis/post-summary-agg-spatial-neonic.RData')  # summary_agg_neonic

# organize data
df_agg_bcr <- summary_agg_neonic %>%
  dplyr::select(year, df_bcr) %>% unnest() %>%
  left_join(bcr_dat, by = 'bcr')

# empirical z-score difference in mean neonic effect between modern period and pre-neonic period
z_empirical <- df_agg_bcr %>%
  dplyr::select(year, bcr_short, theta_mean) %>%
  group_by(bcr_short) %>%
  mutate(mean_pre1 = mean(theta_mean[year %in% c(1980, 1985, 1990)]),
         sd_pre1 = sd(theta_mean[year %in% c(1980, 1985, 1990)]),
         mean_pre2 = mean(theta_mean[year != 2010]),
         sd_pre2 = sd(theta_mean[year != 2010])) %>%
  ungroup() %>% filter(year == 2010) %>%
  mutate(z1 = (theta_mean - mean_pre1) / sd_pre1,
         z2 = (theta_mean - mean_pre2) / sd_pre2) %>%
  arrange(z1)

# Randomization test
RandomizationTest <- function(pre_neonic_period) {

  df_out <- df_agg_bcr %>%
    dplyr::select(year, bcr_short, theta_mean) %>%
    group_by(bcr_short) %>%
    mutate(year = sample(1:7),                                        # randomize year labels within bcr
           mean_pre = mean(theta_mean[year %in% pre_neonic_period]),  # mean in 'pre-neonic' period
           sd_pre = sd(theta_mean[year %in% pre_neonic_period])) %>%  # sd in 'pre-neonic' period
    ungroup() %>%
    filter(year == 7) %>%
    mutate(z = (theta_mean - mean_pre) / sd_pre)   # z-score difference between pre- and post-neonic periods

  return(max(abs(range(df_out$z))))
}

randomization_out1 <- replicate(1000, RandomizationTest(pre_neonic_period = 1:3))
randomization_out2 <- replicate(1000, RandomizationTest(pre_neonic_period = 1:6))

crit1 <- filter(z_empirical, year == 2010, bcr_short == 'Mississippi Valley')$z1 %>% abs()
crit2 <- filter(z_empirical, year == 2010, bcr_short == 'Mississippi Valley')$z2 %>% abs()

length(which(randomization_out1 > crit1)) / length(randomization_out1)   # pre-neonic period is {1975-1984, 1980-1989, 1985-1994}
length(which(randomization_out2 > crit2)) / length(randomization_out2)   # pre-neonic period is all intervals prior to 2005-2014





################# Spatio-temporal analyses (full bayesian model)


# examining difference in BCR-specific neonic effects between pre- and post-neonic periods
FitSpatioTemporal <- function(pre_neonic_period) {
  # organize data
  summary_spp_neonic_sub <- summary_spp_neonic %>% 
    filter(year %in% c(pre_neonic_period, 2010)) %>%
    mutate(post_neonic = ifelse(year == 2010, 1, 0))
  
  # arrange data for stan
  stan_dat <- list(
    N = nrow(summary_spp_neonic_sub),
    K = length(unique(summary_spp_neonic_sub$bcr)),
    M = length(unique(summary_spp_neonic_sub$aou)),
    T = length(unique(summary_spp_neonic_sub$year)),
    est_mean = summary_spp_neonic_sub$gamma_mean,
    est_sigma = summary_spp_neonic_sub$gamma_se,
    K_index = as.numeric(as.factor(summary_spp_neonic_sub$bcr)),
    M_index = as.numeric(as.factor(summary_spp_neonic_sub$aou)),
    T_index = as.numeric(as.factor(summary_spp_neonic_sub$year)),
    period = summary_spp_neonic_sub$post_neonic
  )
  
  # stan parameters
  pars_null <- c(
    'alpha', 'theta_M', 'theta_K',
    'sigma_theta_M', 'sigma_theta_K', 'log_lik', 'theta_hat_bcr_t'
  )
  
  pars_test <- c(
    'alpha', 'theta_M', 'theta_K',
    'sigma_theta_M', 'sigma_theta_K',
    'beta', 'mu_beta', 'sigma_beta', 'log_lik', 'theta_hat_bcr_t'
  )
  
  # fit null and alternative models
  fit_null <- stan(
    file = 'stan/aggregation-multivar-year-null.stan',
    data = stan_dat,
    pars = pars_null,
    warmup = 1500,
    iter = 4000,
    thin = 2,
    chains = 2,
    control = list(adapt_delta = 0.99, stepsize = 0.01)
  )
  
  fit_test <- stan(
    file = 'stan/aggregation-multivar-year-test.stan',
    data = stan_dat,
    pars = pars_test,
    warmup = 1500,
    iter = 4000,
    thin = 2,
    chains = 2,
    control = list(adapt_delta = 0.99, stepsize = 0.01)
  )
  
  # extract log-likelihoods
  log_lik_null <- extract_log_lik(fit_null)
  log_lik_test <- extract_log_lik(fit_test)
  
  # extract bcr-specific estimates of beta
  beta <- rstan::extract(fit_test, 'beta')$beta
  
  return(tibble(ll_null = list(log_lik_null), ll_test = list(log_lik_test), beta = list(beta)))
}

# fit null and alternative models, based on two different definitions of the pre-neonicotinoid period
spatiotemporal1 <- FitSpatioTemporal(pre_neonic_period = c(1980, 1985, 1990))
spatiotemporal2 <- FitSpatioTemporal(pre_neonic_period = c(1980, 1985, 1990, 1995, 2000, 2005))

# save(spatiotemporal1, file = 'analysis/post-summary-agg-spatialtemporal-1.RData')
# save(spatiotemporal2, file = 'analysis/post-summary-agg-spatialtemporal-2.RData')

# load('analysis/post-summary-agg-spatialtemporal-1.RData')  # spatiotemporal1
# load('analysis/post-summary-agg-spatialtemporal-2.RData')  # spatiotemporal2

# perform leave-one out cross validation
loo_null1 <- loo(spatiotemporal1$ll_null[[1]])
loo_test1 <- loo(spatiotemporal1$ll_test[[1]])

loo_null2 <- loo(spatiotemporal2$ll_null[[1]])
loo_test2 <- loo(spatiotemporal2$ll_test[[1]])

# perform null and alternative models, by definition of the pre-neonic period
compare(loo_null1, loo_test1)
compare(loo_null2, loo_test2)


### plot change in neonic effect between pre- and post-neonic periods, by bcr
# load plotting libraries
library(ggplot2)
library(gridExtra)

# extract betas for each definition of pre-neonic period, and apply geometric transform
beta_1 <- 100 * (exp(spatiotemporal1$beta[[1]]) - 1)  # pre-neonic period definition 1
beta_2 <- 100 * (exp(spatiotemporal2$beta[[1]]) - 1)  # pre-neonic period definition 2

# organize dataframes for plotting
df_plot_1 <- data.frame(bcr = sort(unique(summary_spp_neonic$bcr)),   # pre-neonic period definition 1
                        pb0 =   apply(beta_1, 2, function(x)  length(which(x < 0)) / length(x)),
                        pg0 =   apply(beta_1, 2, function(x)  length(which(x > 0)) / length(x)),
                        med =   apply(beta_1, 2, function(x) quantile(x, 0.500)),
                        low90 = apply(beta_1, 2, function(x) quantile(x, 0.050)),
                        upp90 = apply(beta_1, 2, function(x) quantile(x, 0.950)),
                        low99 = apply(beta_1, 2, function(x) quantile(x, 0.005)),
                        upp99 = apply(beta_1, 2, function(x) quantile(x, 0.995))) %>%
  left_join(dplyr::select(bcr_dat, BCR, bcr_short), by = c('bcr' = 'BCR')) %>% 
  arrange(pb0) %>% 
  mutate(bcr_short = factor(bcr_short, levels = bcr_short[order(med)]))

df_plot_2 <- data.frame(bcr = sort(unique(summary_spp_neonic$bcr)),   # pre-neonic period definition 2
                        pb0 =   apply(beta_2, 2, function(x)  length(which(x < 0)) / length(x)),
                        pg0 =   apply(beta_2, 2, function(x)  length(which(x > 0)) / length(x)),
                        med =   apply(beta_2, 2, function(x) quantile(x, 0.500)),
                        low90 = apply(beta_2, 2, function(x) quantile(x, 0.050)),
                        upp90 = apply(beta_2, 2, function(x) quantile(x, 0.950)),
                        low99 = apply(beta_2, 2, function(x) quantile(x, 0.005)),
                        upp99 = apply(beta_2, 2, function(x) quantile(x, 0.995))) %>%
  left_join(dplyr::select(bcr_dat, BCR, bcr_short), by = c('bcr' = 'BCR')) %>% 
  arrange(pb0) %>% 
  mutate(bcr_short = factor(bcr_short, levels = bcr_short[order(med)]))

# first panel, based on pre-neonic period definition 1
p1 <- ggplot(df_plot_1) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5, alpha = 0.3) +
  geom_point(aes(y = bcr_short, x = med), col = 'black', shape = 16, size = 2, alpha = 1) +
  geom_errorbarh(aes(y = bcr_short, x = med, xmin = low90, xmax = upp90), height = 0, size = 0.9) +
  geom_errorbarh(aes(y = bcr_short, x = med, xmin = low99, xmax = upp99), height = 0, size = 0.25) +
  xlab(expression(paste(Delta, Neonicotinoid~effect, ' (', italic(b[k]), ')'))) +
  ylab(NULL) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        text = element_text(size = 13),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.2, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .6, 0, 0, unit = 'cm'))) +
  theme(axis.text.x = element_text(size = 10.5, color = 'grey30'),
        axis.text.y = element_text(size = 9.5, color = 'grey30'),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.2, 'lines'),
        plot.margin = unit(c(0.5, 0.5, 1.1, 0.5), 'lines'))

# second panel, based on pre-neonic period definition 2
p2 <- ggplot(df_plot_2) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5, alpha = 0.3) +
  geom_point(aes(y = bcr_short, x = med), col = 'black', shape = 16, size = 2, alpha = 1) +
  geom_errorbarh(aes(y = bcr_short, x = med, xmin = low90, xmax = upp90), height = 0, size = 0.9) +
  geom_errorbarh(aes(y = bcr_short, x = med, xmin = low99, xmax = upp99), height = 0, size = 0.25) +
  scale_x_continuous(breaks = seq(-3, 2, 1)) +
  xlab(expression(paste(Delta, Neonicotinoid~effect, ' (', italic(b[k]), ')'))) +
  ylab(NULL) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        text = element_text(size = 13),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.2, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .6, 0, 0, unit = 'cm'))) +
  theme(axis.text.x = element_text(size = 10.5, color = 'grey30'),
        axis.text.y = element_text(size = 9.5, color = 'grey30'),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.2, 'lines'),
        plot.margin = unit(c(0.5, 0.5, 1.1, 0.5), 'lines'))

# arrange full plot
p_full <- arrangeGrob(p1, p2, nrow = 1)

# view full plot
dev.off()
quartz(height = 6.5, width = 9)
grid.arrange(p_full)

# write to file
ggsave('figures/appendix-delta-neonic-bcr.png', p_full, height = 6.5, width = 9, units = 'in', dpi = 300)


