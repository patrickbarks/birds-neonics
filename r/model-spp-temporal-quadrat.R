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
    spp-model-temporal-quadrat.R [-c <cores> -s <species>]

  Options:
    -c Number of cores to use [default: 4]
    -s Which portion of spp list to run [default: full]
')


# set options specific to macbook...
if(any(grep('Darwin', Sys.info()))) {
  setwd('~/desktop/bbs/')
  opts$c <- '2'
}


# load relevant data
load('data/bbs-use.RData')
species_df <- read_csv('data/species-list.csv')
route_quads <- read_csv('data/route-quadrat.csv')


# focal year range
years_focal <- 2005:2014


# which portion of species list to run?
indices <- 1:nrow(species_df)
indices_cuts2 <- cut(indices, 2)

if (opts$s == 'full') indices <- indices  # full list
if (opts$s == 'half1') indices <- indices[indices_cuts2 == levels(indices_cuts2)[1]]  # first half
if (opts$s == 'half2') indices <- indices[indices_cuts2 == levels(indices_cuts2)[2]]  # second half


# for each species...
for (i in indices) {
  
  # get aou
  # i <- sample(1:nrow(species_df), 1)  # if want to do test run on a single species
  aou_focal <- species_df$aou[i]
  
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
  
  # # summary
  # nrow(bbs_filter)
  # length(unique(bbs_filter$quad_poly_id))
  # filter(species_df, aou == aou_focal)
  
  # arrange data in list for stan
  stan_dat <- list(
    ncounts = nrow(bbs_filter),
    nquadrats = length(unique(bbs_filter$quad_poly_id)),
    nobservers = length(unique(bbs_filter$obs_route)),
    nyears = length(unique(bbs_filter$year)),
    bird_count = bbs_filter$count,
    quad = as.numeric(as.factor(bbs_filter$quad_poly_id)),
    year = as.numeric(as.factor(bbs_filter$year)),
    obser = as.numeric(as.factor(bbs_filter$obs_route)),
    firstyr = bbs_filter$firstyr
  )
  
  # model pars
  stan_pars <- c('eta', 'int_quadrat', 'beta_quadrat', 'sd_obs', 'sd_overdisp', 'beta_geom')

  # fit stan model
  fit <- stan(
    file = 'stan/spp-temporal-quadrat.stan',
    data = stan_dat,
    par = stan_pars,
    warmup = 2000,
    iter = 3000,
    thin = 2,
    chains = as.numeric(opts$c)
  )
  
  # save and remove stanfit object
  file_out <- paste0('stanfit/temporal-quadrat/fit-temporal-quadrat-', aou_focal, '-', min(years_focal), 'to', max(years_focal), '.RData')
  save(fit, file = file_out)
  rm(fit)
}


