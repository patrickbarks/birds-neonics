#!/usr/bin/Rscript

# load libraries
require(rstan)
require(dplyr)
require(readr)
require(docopt)

# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# docopt
opts <- docopt('
  Usage:
  stan-spatial-cropland.R [-y <year> -c <cores> -s <species>]

  Options:
  -y Year range [default: 2005]
  -c Number of cores to use [default: 4]
  -s Which portion of spp list to run [default: full]
')

# set options specific to my macbook...
if(any(grep('Darwin', Sys.info()))) {
  setwd('~/birds-neonics/')
  opts$c <- '2'
}

# load bbs
load('data/bbs-use.RData')

# params
if (opts$y == '2005') years_focal <- 2005:2014
if (opts$y == '2000') years_focal <- 2000:2009
if (opts$y == '1995') years_focal <- 1995:2004
if (opts$y == '1990') years_focal <- 1990:1990
if (opts$y == '1985') years_focal <- 1985:1994
if (opts$y == '1980') years_focal <- 1980:1989
if (opts$y == '1975') years_focal <- 1975:1984

# species of interest
species_df <- read_csv('data/species-list.csv')

# read landcover data
landcover <- read_csv('data/landcover-route.csv') %>% 
  dplyr::select(route_uniq, pct_crop) %>% 
  filter(route_uniq %in% bbs_use$route_uniq) %>% 
  mutate(crop_cubert = pct_crop^(1/3)) %>% 
  mutate(crop_scale = (crop_cubert - mean(crop_cubert)) / sd(crop_cubert))

# which portion of species list to run?
indices <- 1:nrow(species_df)
indices_cuts2 <- cut(indices, 2)

if (opts$s == 'full') indices <- indices  # full list
if (opts$s == 'half1') indices <- indices[indices_cuts2 == levels(indices_cuts2)[1]]  # first half
if (opts$s == 'half2') indices <- indices[indices_cuts2 == levels(indices_cuts2)[2]]  # second half


## for each species...
for (i in indices) {

  # get aou
  # i <- sample(indices, 1)  # if want to do test run on a single species
  # aou_focal <- '06050'
  aou_focal <- species_df$aou[i]
  
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
  
  # route by bcr and cropland
  route_by_bcr <- bbs_filter %>%
    dplyr::select(route_uniq, bcr) %>% unique() %>% 
    left_join(landcover, by = 'route_uniq') %>% 
    group_by(bcr) %>% 
    mutate(crop_mean_cent = crop_scale - mean(crop_scale)) %>% ungroup() %>% 
    mutate(route_no = as.numeric(as.factor(route_uniq)), bcr_no = as.numeric(as.factor(bcr))) %>% 
    arrange(route_no)
  
  # # summary
  # nrow(bbs_filter)
  # length(unique(bbs_filter$bcr))
  # length(unique(bbs_filter$route_uniq))
  # filter(species_df, aou == aou_focal)
  
  # arrange data in list
  stan_dat <- list(
    ncounts = nrow(bbs_filter),
    nroutes = length(unique(bbs_filter$route_uniq)),
    nbcrs = length(unique(bbs_filter$bcr)),
    nobservers = length(unique(bbs_filter$observer)),
    nyears = length(unique(bbs_filter$year)),
    bird_count = bbs_filter$count,
    route = as.numeric(as.factor(bbs_filter$route_uniq)),
    bcr = as.numeric(as.factor(bbs_filter$bcr)),
    year = as.numeric(as.factor(bbs_filter$year)),
    observer = as.numeric(as.factor(bbs_filter$observer)),
    firstyr = bbs_filter$firstyr,
    route_by_bcr = route_by_bcr$bcr_no,
    cropland = route_by_bcr$crop_mean_cent
  )
  
  # model pars and stan file
  stan_pars <- c('gamma', 'eta', 'int_bcr', 'beta_bcr',
                 'dev_int_route', 'dev_beta_route', 'dev_beta_route_hat',
                 'mu_int_bcr', 'mu_beta_bcr', 'sd_int_bcr', 'sd_beta_bcr',
                 'sd_int_route', 'sd_beta_route', 'sd_obs', 'sd_overdisp')
  
  # stan fit
  fit <- stan(
    file = 'stan/spp-spatial-cropland.stan',
    data = stan_dat,
    par = stan_pars,
    warmup = 2000,
    iter = 3000,
    thin = 2,
    chains = as.numeric(opts$c)
  )

  # check for divergent transitions or rhat high
  diverg <- do.call(rbind, args = get_sampler_params(fit, inc_warmup = F))[,5]
  n_diverg <- length(which(diverg == 1))
  df_summ <- as.data.frame(summary(fit)$summary)
  rhat_high <- length(which(df_summ$Rhat > 1.1))
  
  # save and remove stanfit object
  file_out <- paste0('stanfit/spatial-cropland/fit-cropland-route-nested4-', aou_focal, '-', min(years_focal), 'to', max(years_focal), '.RData')
  save(fit, file = file_out)
  rm(fit)
  
  # rerun if necessary
  if(n_diverg > 30 | rhat_high > 0) {
    
    fit <- stan(
      file = 'stan/spp-spatial-cropland.stan',
      data = stan_dat,
      par = stan_pars,
      warmup = 2000,
      iter = 4000,
      thin = 2,
      chains = as.numeric(opts$c),
      control = list(adapt_delta = 0.99, stepsize = 0.01)
    )
    
    file_out_redo <- paste0('stanfit/spatial-cropland/fit-cropland-route-nested4-redo-', aou_focal, '-', min(years_focal), 'to', max(years_focal), '.RData')
    save(fit, file = file_out_redo)
    rm(fit)
  }
}


# ################### Diagnostics
# library(shinystan)
# launch_shinystan(fit)

