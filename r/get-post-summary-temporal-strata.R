
# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(rstan)

# setwd
setwd('~/birds-neonics/')

# read species of interest
species_df <- read_csv('data/species-list.csv')

# create df of relevant stanfit files
path <- 'stanfit/temporal-strata/'
fits <- list.files(path)

fit_aous <- sapply(fits, function (x) strsplit(x, '-|\\.')[[1]][4], USE.NAMES = F)
fit_df <- tibble(aou = fit_aous, file = fits, path = path)



############## Get species abundance by year over the period 1985-2014
GetAbundance85to14 <- function(path, aou, file) {
  
  load(paste0(path, file))
  CompIndex <- rstan::extract(fit, pars = 'CompIndex', permuted = T)$CompIndex
  
  diverg <- do.call(rbind, args = get_sampler_params(fit, inc_warmup = F))[,5]
  n_diverg <- length(which(diverg == 1))
  df_summ <- as.data.frame(summary(fit)$summary)
  rhat_high <- length(which(df_summ$Rhat > 1.1))
  
  df <- data.frame(year = 1985:2014,
                   comp_med = apply(CompIndex, 2, function(x) quantile(x, 0.500)),
                   comp_low90 = apply(CompIndex, 2, function(x) quantile(x, 0.05)),
                   comp_upp90 = apply(CompIndex, 2, function(x) quantile(x, 0.95)),
                   n_diverg = n_diverg, rhat_high = rhat_high)
  
  return(df)
}

# get posterior summary for species abundance by year over the period 1985-2014
summary_abundance_85_14 <- group_by(fit_df, aou) %>% 
  do(GetAbundance85to14(path = .$path, aou = .$aou, file = .$file)) %>%
  ungroup()

# diagnostics
summary_abundance_85_14 %>% 
  group_by(aou) %>% 
  summarize(n_diverg = unique(n_diverg),
            rhat_high = unique(rhat_high)) %>% 
  as.data.frame()

# write to file
# write.csv(summary_abundance_85_14, 'analysis/post-summary-spp-temporal-strata-abundance.csv', row.names = F)




############## Get population trends over N-year interval[s]
GetTrendsHelper <- function(path, aou, file, year_indices_l) {
  load(paste0(path, file))
  CompIndex <- rstan::extract(fit, pars = 'CompIndex', permuted = T)$CompIndex
  
  out <- lapply(year_indices_l, GetTrends, CompIndex = CompIndex) %>%
    do.call(rbind.data.frame, .)
  
  return(out)
}

GetTrends <- function(year_indices, CompIndex) {
  years <- 1985:2014
  
  coefs <- apply(CompIndex[,year_indices], 1, GetBetaGeometric, len_x = length(year_indices)) %>%
    do.call(rbind.data.frame, .) %>%
    rownames_to_column(var = 'rep')
  
  beta_df <- data.frame(
    year_start = min(years[year_indices]),
    year_end = max(years[year_indices]),
    beta_mean = mean(coefs$beta),
    beta_se = sd(coefs$beta),
    beta_med = as.numeric(quantile(coefs$beta, 0.500)),
    beta_low90 = as.numeric(quantile(coefs$beta, 0.050)),
    beta_upp90 = as.numeric(quantile(coefs$beta, 0.950)),
    beta_low99 = as.numeric(quantile(coefs$beta, 0.005)),
    beta_upp99 = as.numeric(quantile(coefs$beta, 0.995))
  )
  
  return(beta_df)
}

GetBetaGeometric <- function(y_raw, len_x) {
  x <- 1:len_x
  y <- log(y_raw)
  xhat = mean(x)
  yhat = mean(y)
  bhat = sum((x - xhat) * (y - yhat)) / sum((x - xhat)^2)
  ahat = yhat - bhat * xhat
  bhat_geom = 100 * (exp(bhat) - 1)
  return(data.frame(alpha = ahat, beta = bhat_geom))
}

# create lists of relevant years over which to calculate population trends
years_list_full <- list(1:30)                             # 1985 to 2014
years_list_7 <- lapply(1:24, function(x) seq(x, x + 6))   # 7-year intervals

# get posterior summary for population trends over period 1985-2014
summary_trend_allyears <- group_by(fit_df, aou) %>%
  do(GetTrendsHelper(path = .$path, aou = .$aou, file = .$file, year_indices_l = years_list_full)) %>%
  ungroup()

# get posterior summary for population trends in 7-year intervals (moving window)
summary_trend_7years <- group_by(fit_df, aou) %>%
  do(GetTrendsHelper(path = .$path, aou = .$aou, file = .$file, year_indices_l = years_list_7)) %>%
  ungroup()

# write to file
# write.csv(summary_trend_allyears, 'analysis/post-summary-spp-temporal-strata-trend-allyears.csv', row.names = F)
# write.csv(summary_trend_7years, 'analysis/post-summary-spp-temporal-strata-trend-7years.csv', row.names = F)

