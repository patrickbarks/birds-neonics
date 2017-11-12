

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(rstan)


# setwd
setwd('~/desktop/bbs/')


# read relevant data
species_df <- read_csv('data/species-list.csv')
route_quads <- read_csv('data/route-quadrat.csv')


# load bbs data
load('data/bbs-use.RData')


# create list of relevant stan files
path <- 'stanfit/temporal-quadrat/'
fit_file <- list.files(path)


# extract aou and year from filenames
ParseString <- function(x, index) { strsplit(x, '-|\\.')[[1]][index] }

fit_aou <- sapply(fit_file, ParseString, index = 4, USE.NAMES = F)
fit_year <- sapply(fit_file, ParseString, index = 5, USE.NAMES = F) %>% 
  substr(1, 4) %>% as.numeric() + 5


# arrange file info into dataframe
fit_df <- data.frame(aou = fit_aou,
                     year = fit_year,
                     file = fit_file,
                     path = path,
                     stringsAsFactors = F)


# function to extract posterior summaries of interest from rstan objects
GetPars <- function(path, aou, year, file) {
  aou_focal <- aou
  years_focal <- seq(year - 5, year + 4)
  
  # filter bbs data for species of interest
  bbs_filter <- bbs_use %>% 
    filter(aou == aou_focal, runtype == 1, year %in% years_focal) %>% 
    left_join(route_quads, by = 'route_uniq') %>% 
    group_by(quad_poly_id) %>%        # remove quadrats lacking nonzero counts
    mutate(allzero = ifelse(all(count == 0), T, F)) %>%
    filter(allzero == F) %>% 
    mutate(n = n()) %>%               # remove quadrats with < 3 counts (incl. zeros)
    filter(n >= 3) %>%
    ungroup()
  
  load(paste0(path, file))
  
  diverg <- do.call(rbind, args = get_sampler_params(fit, inc_warmup = F))[,5]
  n_diverg <- length(which(diverg == 1))
  df_summ <- as.data.frame(summary(fit)$summary)
  rhat_high <- length(which(df_summ$Rhat > 1.1))
  n_eff_low <- length(which(df_summ$n_eff / length(diverg) < 0.1))
  mcse_high <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
  
  n_quadrats <- length(unique(bbs_filter$quad_poly_id))
  
  beta <- rstan::extract(fit, pars = 'beta_quadrat', permuted = T)$beta_quadrat
  beta_mean <- apply(beta, 2, mean)
  beta_se <- apply(beta, 2, sd)
  
  df <- data.frame(
    aou, year, quadrat = sort(unique(bbs_filter$quad_poly_id)),
    n_diverg, rhat_high, mcse_high, n_eff_low, n_quadrats,
    beta_mean, beta_se, stringsAsFactors = F
  )
  
  return(df)
}


# get posterior summaries
summary_temporal_quadrat <- fit_df %>% 
  group_by(aou, year) %>% 
  do(GetPars(path = .$path, year = .$year, aou = .$aou, file = .$file)) %>%
  ungroup()


# diagnostics summary
summary_temporal_quadrat %>%
  group_by(aou, year) %>%
  summarize(n_diverg = unique(n_diverg),
            rhat_high = unique(rhat_high),
            mcse_high = unique(mcse_high),
            n_quadrats = unique(n_quadrats)) %>%
  as.data.frame()


# write.csv(summary_temporal_quadrat, 'analysis/post-summary-spp-temporal-quadrat.csv', row.names = F)



# ### diagnostics
# library(shinystan)
# load('stanfit/temporal-quadrat/fit-temporal-quadrat-07660-2005to2014.RData')
# launch_shinystan(fit)


