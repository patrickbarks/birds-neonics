
# load libraries
library(dplyr)
library(readr)
library(docopt)
library(tidyr)
library(rstan)

# setwd
setwd('~/birds-neonics/')

# read species of interest
species_df <- read_csv('data/species-list.csv')

# load bbs data
load('data/bbs-use.RData')

# read cropland data
landcover <- read_csv('data/landcover-route.csv') %>% 
  dplyr::select(route_uniq, pct_crop) %>% 
  filter(route_uniq %in% bbs_use$route_uniq) %>% 
  mutate(crop_cubert = pct_crop^(1/3)) %>% 
  mutate(crop_scale = (crop_cubert - mean(crop_cubert)) / sd(crop_cubert))

# create list of relevant stan files
path <- 'stanfit/spatial-cropland/'
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
    group_by(route_uniq) %>%      # remove routes lacking nonzero counts
    mutate(allzero = ifelse(all(count == 0), T, F)) %>%
    filter(allzero == F) %>% 
    mutate(n = n()) %>%           # remove routes with < 3 counts (incl. zeros)
    filter(n >= 3) %>%
    ungroup() %>%
    group_by(bcr) %>%             # remove BCR with < 3 routes
    mutate(n_routes = length(unique(route_uniq))) %>%
    filter(n_routes >= 3) %>% 
    ungroup()
  
  bcr_df <- data.frame(bcr = sort(unique(bbs_filter$bcr))) %>% 
    mutate(bcr_no = as.numeric(as.factor(bcr)))
  
  load(paste0(path, file))
  
  diverg <- do.call(rbind, args = get_sampler_params(fit, inc_warmup = F))[,5]
  n_diverg <- length(which(diverg == 1))
  df_summ <- as.data.frame(summary(fit)$summary)
  rhat_high <- length(which(df_summ$Rhat > 1.1))
  n_eff_low <- length(which(df_summ$n_eff / length(diverg) < 0.1))
  mcse_high <- length(which(df_summ$se_mean / df_summ$sd > 0.1))
  
  n_routes <- length(unique(bbs_filter$route_uniq))
  n_bcr <- length(unique(bbs_filter$bcr))
  
  # note gamma was called 'g1' in the code run on AWS (but otherwise same model)
  gamma <- rstan::extract(fit, pars = 'g1', permuted = T)$g1
  gamma_mean <- apply(gamma, 2, mean)
  gamma_se <- apply(gamma, 2, sd)
  
  df <- data.frame(
    aou, year, bcr = bcr_df$bcr,
    n_diverg, rhat_high, mcse_high, n_eff_low, n_routes, n_bcr,
    gamma_mean, gamma_se, stringsAsFactors = F
  )
  
  return(df)
}

# get posterior summaries
summary_cropland <- fit_df %>% 
  group_by(aou, year) %>% 
  do(GetPars(path = .$path, year = .$year, aou = .$aou, file = .$file)) %>%
  ungroup()

# diagnostics summary
summary_cropland %>%
  group_by(aou, year) %>%
  summarize(n_diverg = unique(n_diverg),
            rhat_high = unique(rhat_high),
            n_bcr = unique(n_bcr),
            n_routes = unique(n_routes)) %>%
  as.data.frame()

# write.csv(summary_cropland, 'analysis/post-summary-spp-spatial-cropland.csv', row.names = F)

