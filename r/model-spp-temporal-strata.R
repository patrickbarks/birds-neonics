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
    spp-model-temporal-strata.R [-c <cores> -s <species>]

  Options:
    -c Number of cores to use [default: 4]
    -s Which portion of spp list to run [default: full]
')

# set options specific to my macbook...
if(any(grep('Darwin', Sys.info()))) {
  setwd('~/birds-neonics/')
  opts$c <- '2'
}


# load bbs data
load('data/bbs-use.RData')

# read species list
species_df <- read_csv('data/species-list.csv')

# read strata info
strat_info <- read_csv('data/bbs_raw_2015/bbs_strata_info.csv') %>% 
  dplyr::select(stratum = strat_no, strat_area)

# which portion of species list to run?
indices <- 1:nrow(species_df)
indices_cuts2 <- cut(indices, 2)

if (opts$s == 'full') indices <- indices  # full list
if (opts$s == 'half1') indices <- indices[indices_cuts2 == levels(indices_cuts2)[1]]  # first half
if (opts$s == 'half2') indices <- indices[indices_cuts2 == levels(indices_cuts2)[2]]  # second half

# focal year range
years_focal <- 1985:2014


## for each species...
for (i in indices) {
  
  # get aou
  # i <- sample(1:nrow(species_df), 1)  # if want to do test run on a single species
  aou_focal <- species_df$aou[i]
  
  # filter bbs data
  bbs_filter1 <- bbs_use %>% 
    filter(aou == aou_focal, runtype == 1, year %in% years_focal) %>%
    group_by(stratum) %>%        # remove strata lacking nonzero counts
    mutate(strat_zero = ifelse(all(count == 0), T, F)) %>%
    filter(strat_zero == F) %>% 
    ungroup()
  
  # find routes with non-zero counts
  route_nonzero <- bbs_filter1 %>% 
    group_by(stratum) %>% 
    summarize(route_nonzero = length(unique(route_uniq[which(count > 0)])),
              route_total = length(unique(route_uniq)),
              prop_nonzero = route_nonzero/route_total) %>% 
    left_join(strat_info, by = 'stratum')
  
  # remove routes lacking non-zero counts
  bbs_filter2 <- bbs_filter1 %>% 
    group_by(route_uniq) %>%
    mutate(route_zero = ifelse(all(count == 0), T, F)) %>%
    filter(route_zero == F) %>% 
    ungroup() %>% 
    mutate(strat_num = as.numeric(as.factor(stratum)), year_num = as.numeric(as.factor(year)))

  # every combination of year and strata with bird counts (zero or otherwise)
  year_by_strat <- bbs_filter2 %>% 
    dplyr::select(year_num, strat_num) %>%
    unique() %>% arrange(strat_num, year_num) %>% 
    mutate(year_strat_index = 1:n())
  
  # join year_by_strat to rest of bbs data
  bbs_filter3 <- bbs_filter2 %>% 
    left_join(year_by_strat, by = c('strat_num', 'year_num'))
  
  # # summary
  # nrow(bbs_filter3)
  # length(unique(bbs_filter3$stratum))
  # filter(species_df, aou == aou_focal)
  
  # arrange data in list
  stan_dat <- list(
    ncounts = nrow(bbs_filter3),
    nstrata = length(unique(bbs_filter3$stratum)),
    nobservers = length(unique(bbs_filter3$obs_route)),
    nyears = length(unique(bbs_filter3$year)),
    bird_count = bbs_filter3$count,
    strat = bbs_filter3$strat_num,
    year = bbs_filter3$year_num,
    year_int = bbs_filter3$year_num,
    year_uniq = sort(unique(bbs_filter3$year_num)),
    obser = as.numeric(as.factor(bbs_filter3$obs_route)),
    firstyr = bbs_filter3$firstyr,
    area = route_nonzero$strat_area,
    nonzero = route_nonzero$prop_nonzero,
    n_year_strat = nrow(year_by_strat),
    year_by_strat = as.matrix(year_by_strat[,1:2]),
    year_strat_index = bbs_filter3$year_strat_index
  )
  
  # model pars
  stan_pars <- c('yeareffect', 'eta', 'sd_obs', 'sd_overdisp', 'n', 'N', 'CompIndex')
  
  # stanfit
  fit <- stan(
    file = 'stan/spp-temporal-strata.stan',
    data = stan_dat, 
    par = stan_pars,
    warmup = 2000,
    iter = 3000,
    thin = 2,
    chains = as.numeric(opts$c)
  )
  
  # save and remove stanfit object
  file_out <- paste0('stanfit/temporal-strata/fit-temporal-strata-', aou_focal, '-', min(years_focal), 'to', max(years_focal), '.RData')
  save(fit, file = file_out)
  rm(fit)
}


# ################### Diagnostics
# library(shinystan)
# launch_shinystan(fit)

