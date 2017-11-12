

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(rstan)


# setwd
setwd('~/desktop/bbs/')


# species of interest
species_df <- read_csv('data/species-list.csv')


# load bbs
load('data/bbs-use.RData')


# read county data
route_counties <- read_csv('data/route-county-start.csv') %>%
  dplyr::select(route_uniq, fips)


# create list of relevant stanfit files
path <- 'stanfit/spatial-neonic/'
fit_file <- list.files(path)


# extract aou and year from filenames
ParseString <- function(x, index) { strsplit(x, '-|\\.')[[1]][index] }

fit_aou <- sapply(fit_file, ParseString, index = 5, USE.NAMES = F)
fit_year <- sapply(fit_file, ParseString, index = 6, USE.NAMES = F) %>% 
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
    left_join(route_counties, by = 'route_uniq') %>% 
    mutate(fips_bcr = paste(fips, bcr, sep = '_')) %>% 
    group_by(fips_bcr) %>%        # remove counties lacking nonzero counts
    mutate(allzero = ifelse(all(count == 0), T, F)) %>%
    filter(allzero == F) %>% 
    mutate(n = n()) %>%           # remove counties with < 3 counts (incl. zeros)
    filter(n >= 3) %>%
    ungroup() %>%
    group_by(bcr) %>%             # remove BCR with < 3 counties
    mutate(n_fips = length(unique(fips_bcr))) %>%
    filter(n_fips >= 3) %>% 
    ungroup()
  
  load(paste0(path, file))
  
  diverg <- do.call(rbind, args = get_sampler_params(fit, inc_warmup = F))[,5]
  n_diverg <- length(which(diverg == 1))
  df_summ <- as.data.frame(summary(fit)$summary)
  rhat_high <- length(which(df_summ$Rhat > 1.1))
  n_eff_low <- length(which(df_summ$n_eff / length(diverg) < 0.1))
  mcse_high <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
  
  n_counties <- length(unique(bbs_filter$fips_bcr))
  n_bcr <- length(unique(bbs_filter$bcr))
  
  # note gamma was called 'g1' in the code run on AWS (but otherwise same model)
  gamma <- rstan::extract(fit, pars = 'g1', permuted = T)$g1 
  gamma_mean <- apply(gamma, 2, mean)
  gamma_se <- apply(gamma, 2, sd)
  
  df <- data.frame(
    aou, year, bcr = sort(unique(bbs_filter$bcr)),
    n_diverg, rhat_high, mcse_high, n_eff_low, n_counties, n_bcr,
    gamma_mean, gamma_se, stringsAsFactors = F
  )
  
  return(df)
}


# get posterior summaries
summary_neonic <- fit_df %>% 
  group_by(aou, year) %>% 
  do(GetPars(path = .$path, year = .$year, aou = .$aou, file = .$file)) %>%
  ungroup()


# diagnostics summary
summary_neonic %>%
  group_by(aou, year) %>%
  summarize(n_diverg = unique(n_diverg),
            rhat_high = unique(rhat_high),
            n_bcr = unique(n_bcr),
            n_counties = unique(n_counties)) %>%
  as.data.frame()


# write.csv(summary_neonic, 'analysis/post-summary-spp-spatial-neonic.csv', row.names = F)


