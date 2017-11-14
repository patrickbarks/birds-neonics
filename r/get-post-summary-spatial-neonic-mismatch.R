
# load libraries
library(dplyr)
library(readr)
library(docopt)
library(tidyr)
library(rstan)

# setwd
setwd('~/birds-neonics/')

# species of interest
species_df <- read_csv('data/species-list.csv')


# load bbs
load('data/bbs-use.RData')

# read neonic data
neonics_counties <- read_csv('data/neonic-county-mean-2005-2012.csv') %>% 
  mutate(neonic_cubert = neonic_mean^(1/3)) %>% dplyr::select(-neonic_mean)

# merge county and neonic data
route_counties <- read_csv('data/route-county-start.csv') %>% 
  left_join(neonics_counties, by = 'fips') %>% 
  group_by(route_uniq) %>% 
  summarize(fips_merge = paste(fips, collapse = '-'),
            neonic_cubert = mean(neonic_cubert))

neonics_counties_merge <- route_counties %>% 
  group_by(fips_merge) %>% 
  summarize(neonic_cubert = unique(neonic_cubert))

# create list of relevant stanfit files
path <- 'stanfit/spatial-neonic-mismatch/'
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


# function to extract parameters of interest from rstan objects
GetPars <- function(path, aou, year, file) {
  aou_focal <- aou
  years_focal <- seq(year - 5, year + 4)
  
  # filter bbs data for species of interest
  bbs_filter <- bbs_use %>% 
    filter(aou == aou_focal, runtype == 1, year %in% years_focal) %>% 
    left_join(route_counties, by = 'route_uniq') %>% 
    mutate(fips_bcr = paste(fips_merge, bcr, sep = '-')) %>% 
    group_by(fips_bcr) %>%        # remove merged counties lacking nonzero counts
    mutate(allzero = ifelse(all(count == 0), T, F)) %>%
    filter(allzero == F) %>% 
    mutate(n = n()) %>%           # remove merged counties with < 3 counts (incl. zeros)
    filter(n >= 3) %>%
    ungroup() %>%
    group_by(bcr) %>%             # remove BCR with < 3 merged counties
    mutate(n_fips = length(unique(fips_bcr))) %>%
    filter(n_fips >= 3) %>% 
    ungroup()
  
  bcr_df <- data.frame(bcr = sort(unique(bbs_filter$bcr))) %>% 
    mutate(bcr_no = as.numeric(as.factor(bcr)))
  
  load(paste0(path, file))
  if(!exists('fit')) fit <- fit_redo
  
  diverg <- do.call(rbind, args = get_sampler_params(fit, inc_warmup = F))[,5]
  n_diverg <- length(which(diverg == 1))
  df_summ <- as.data.frame(summary(fit)$summary)
  rhat_high <- length(which(df_summ$Rhat > 1.1))
  n_eff <- df_summ$n_eff / length(diverg)
  n_eff_low <- length(which(n_eff < 0.1))
  mcse_high <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
  
  n_strata <- length(unique(bbs_filter$fips_bcr))
  n_bcr <- length(unique(bbs_filter$bcr))
  
  gamma <- rstan::extract(fit, pars = 'gamma', permuted = T)$gamma
  gamma <- gamma * sd(neonics_counties_merge$neonic_cubert)  # forgot to initially z-standardize neonic use in spatial mismatch models
  gamma_mean <- apply(gamma, 2, mean)
  gamma_se <- apply(gamma, 2, sd)
  
  df <- data.frame(
    aou, year, bcr = bcr_df$bcr,
    n_diverg, rhat_high, mcse_high, n_eff_low, n_strata, n_bcr,
    gamma_mean, gamma_se, stringsAsFactors = F
  )
  
  return(df)
}

# get posterior summaries
summary_neonic_mismatch <- group_by(fit_df, aou) %>% 
  do(GetPars(path = .$path, year = .$year, aou = .$aou, file = .$file)) %>% ungroup()

# summary
summary_neonic_mismatch %>%
  group_by(aou, year) %>%
  summarize(n_diverg = unique(n_diverg),
            rhat_high = unique(rhat_high),
            mcse_high = unique(mcse_high),
            n_bcr = unique(n_bcr),
            n_strata = unique(n_strata)) %>%
  as.data.frame()

# write.csv(summary_neonic_mismatch, 'analysis/post-summary-spp-spatial-neonic-mismatch.csv', row.names = F)

