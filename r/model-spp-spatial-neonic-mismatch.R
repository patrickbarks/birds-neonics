#!/usr/bin/Rscript


# load libraries
require(rstan)
require(dplyr)
require(readr)
require(docopt)


# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# docopt options
opts <- docopt('
  Usage:
  stan-spatial-neonic-spatialmismatch.R [-y <year> -c <cores> -s <species>]

  Options:
  -y Year range [default: 2005]
  -c Number of cores to use [default: 4]
  -s Which portion of spp list to run [default: full]
')


# set options specific to macbook...
if(any(grep('Darwin', Sys.info()))) {
  setwd('~/desktop/bbs/')
  opts$c <- '2'
}


# load bbs data
load('data/bbs-use.RData')


# read species of interest
species_df <- read_csv('data/species-list.csv')


# read neonic data
neonics_counties <- read_csv('data/neonic-county-mean-2005-2012.csv') %>% 
  mutate(neonic_cubert = neonic_mean^(1/3)) %>% dplyr::select(-neonic_mean)


# merge county and neonic data
route_counties <- read_csv('data/route-county-full.csv') %>% 
  left_join(neonics_counties, by = 'fips') %>% 
  group_by(route_uniq) %>% 
  summarize(fips_merge = paste(fips, collapse = '-'),
            neonic_cubert = mean(neonic_cubert),
            type = unique(type))

neonics_counties_merge <- route_counties %>% 
  group_by(fips_merge) %>% 
  summarize(neonic_cubert = unique(neonic_cubert))


# params
if (opts$y == '2005') years_focal <- 2005:2014



# which portion of species list to run?
indices <- 1:nrow(species_df)
indices_cuts2 <- cut(indices, 2)

if (opts$s == 'full') indices <- indices  # full list
if (opts$s == 'half1') indices <- indices[indices_cuts2 == levels(indices_cuts2)[1]]  # first half
if (opts$s == 'half2') indices <- indices[indices_cuts2 == levels(indices_cuts2)[2]]  # second half


# for each species...
for (i in indices) {

  # get species aou
  # i <- sample(indices, 1)  # if want to do test run on a single species
  aou_focal <- species_df$aou[i]
  
  # filter bbs data for species of interest
  bbs_filter <- bbs_use %>% 
    filter(aou == aou_focal, runtype == 1, year %in% years_focal) %>% 
    left_join(route_counties, by = 'route_uniq') %>% 
    mutate(fips_bcr = paste(fips_merge, bcr, sep = '-')) %>% 
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
  
  # fips by bcr and neonics
  fips_by_bcr <- bbs_filter %>%
    dplyr::select(fips_bcr, bcr, fips_merge) %>% unique() %>% 
    left_join(neonics_counties_merge, by = 'fips_merge') %>% 
    group_by(bcr) %>% 
    mutate(neonic_mean_cent = neonic_cubert - mean(neonic_cubert)) %>% ungroup() %>% 
    mutate(fips_bcr_no = as.numeric(as.factor(fips_bcr)), bcr_no = as.numeric(as.factor(bcr))) %>% 
    arrange(fips_bcr_no)
  
  # # summary
  # nrow(bbs_filter)
  # length(unique(bbs_filter$bcr))
  # length(unique(bbs_filter$fips_bcr))
  # filter(species_df, aou == aou_focal)
  
  # arrange data in list for stan
  stan_dat <- list(
    ncounts = nrow(bbs_filter),
    ncounties = length(unique(bbs_filter$fips_bcr)),
    nbcrs = length(unique(bbs_filter$bcr)),
    nobservers = length(unique(bbs_filter$obs_route)),
    nyears = length(unique(bbs_filter$year)),
    bird_count = bbs_filter$count,
    county = as.numeric(as.factor(bbs_filter$fips_bcr)),
    bcr = as.numeric(as.factor(bbs_filter$bcr)),
    year = as.numeric(as.factor(bbs_filter$year)),
    observer = as.numeric(as.factor(bbs_filter$obs_route)),
    firstyr = bbs_filter$firstyr,
    county_by_bcr = fips_by_bcr$bcr_no,
    neonic = fips_by_bcr$neonic_mean_cent
  )

  # stan parameters
  stan_pars <- c('gamma', 'eta', 'int_bcr', 'beta_bcr',
                 'dev_int_county', 'dev_beta_county', 'dev_beta_county_hat',
                 'mu_int_bcr', 'mu_beta_bcr', 'sd_int_bcr', 'sd_beta_bcr',
                 'sd_int_county', 'sd_beta_county', 'sd_obs', 'sd_overdisp')
  
  # fit stan model
  fit <- stan(
    file = 'stan/spp-spatial-neonicotinoid.stan',
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
  
  # save and remove file
  file_out <- paste0('stanfit/spatial-neonic-mismatch/fit-neonic-mismatch-', aou_focal, '-', min(years_focal), '-', max(years_focal), '.RData')
  save(fit, file = file_out)
  rm(fit)
  
  # rerun if poor convergence/diagnostics
  if(n_diverg > 30 | rhat_high > 0) {
    
    fit <- stan(
      file = 'stan/spp-spatial-neonicotinoid.stan',
      data = stan_dat,
      par = stan_pars,
      warmup = 2000,
      iter = 3000,
      thin = 2,
      chains = as.numeric(opts$c),
      control = list(adapt_delta = 0.99, stepsize = 0.01)
    )
    
    file_out_redo <- paste0('stanfit/spatial-neonic-mismatch/fit-neonic-mismatch-redo-', aou_focal, '-', min(years_focal), '-', max(years_focal), '.RData')
    save(fit, file = file_out_redo)
    rm(fit)
  }

}

