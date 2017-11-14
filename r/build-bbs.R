

# load libraries
library(dplyr)
library(readr)

# setwd
setwd('~/birds-neonics/')

# read data
dt <- read_csv('data/bbs-allstates-2015.txt')        # BBS bird counts
wx <- read.csv('data/bbs_raw_2015/weather.csv')      # BBS weather info
rt <- read_csv('data/bbs_raw_2015/routes.csv')       # BBS route info
rg <- read_csv('data/bbs_raw_2015/RegionCodes.csv')  # BBS region info
species <- read_csv('data/species-list.csv')         # Barks and Galpern species list

# Note that bbs_allstates-2015.txt is a concatenation of state-specific BBS tables found
# within the folder 'States/'. It can by produced manually, or by using the following
# shell commands on Mac OS X:
# mkdir States
# cd States
# ftp -i ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/Archivefiles/Version2015v1/States/
# mget *.zip
# exit
# unzip "*.zip"
# rm *.zip
# head -1 Alabama.csv > ../bbs-allstates-2015.txt
# tail -n +2 -q *.csv >> ../bbs-allstates-2015.txt

# remove extraneous columns
dt <- select(dt, -matches('count[[:digit:]]+'), -StopTotal, count = SpeciesTotal)
wx <- select(wx, countrynum, statenum, Route, RPID, Year, observer = ObsN, RunType)
rt <- select(rt, countrynum, statenum, Route, Stratum, BCR)

# format state, route, and region IDs (to ensure consistency in # of digits)
wx$statenum <- formatC(wx$statenum, width = 2, flag = "0")
wx$Route <- formatC(wx$Route, width = 3, flag = "0")
rt$statenum <- formatC(rt$statenum, width = 2, flag = "0")
rt$Route <- formatC(rt$Route, width = 3, flag = "0")
rg$RegionCode <- formatC(rg$RegionCode, width = 2, flag = "0")

# limit to contiguous USA (exclude Alaska, statenum == 03)
dt <- filter(dt, countrynum == 840, statenum != "03")
wx <- filter(wx, countrynum == 840, statenum != "03")

# create unique route identifiers (paste country-state-route)
dt$route_uniq <- with(dt, paste0("C", countrynum, "S", statenum, "R", Route))
wx$route_uniq <- with(wx, paste0("C", countrynum, "S", statenum, "R", Route))
rt$route_uniq <- with(rt, paste0("C", countrynum, "S", statenum, "R", Route))

# create unique route-year identifiers (paste route-year)
dt$route_yr <- with(dt, paste0(route_uniq, "YR", Year))
wx$route_yr <- with(wx, paste0(route_uniq, "YR", Year))

# create unique route-year-rpid identifiers (RPID indicates whether survey is within-year replicate)
dt$route_yr_rpid <- with(dt, paste0(route_yr, "RPID", RPID))
wx$route_yr_rpid <- with(wx, paste0(route_yr, "RPID", RPID))

# create unique observer-route identifier
wx$obs_route <- with(wx, paste0(observer, route_uniq))

# join region df to route df
rt_rg <- left_join(rt, rg, by = c('countrynum', 'statenum' = 'RegionCode'))

# create observer firstyr var in wx (1 == firstyr, 0 == not-firstyr)
wx <- wx %>% 
  group_by(observer) %>%
  mutate(firstyr = ifelse(Year == min(Year), 1, 0)) %>% 
  ungroup()

# check that number of first-year runs is sensible
table(wx$firstyr) 

# join weather df to route-region df
wx_rt <- left_join(wx, select(rt_rg, route_uniq, state_prov, Stratum, BCR), by = 'route_uniq')

# put together the full dataset (excluding counts of zero)
bbs_full_nozeros <- left_join(select(dt, route_yr_rpid, Aou, count), wx_rt, by = 'route_yr_rpid') %>%
  setNames(tolower(names(.))) %>%     # set column names to lowercase
  select(aou, route_uniq, year, rpid, observer, firstyr, obs_route,
         runtype, count, state_prov, stratum, bcr) %>%     # rearrange columns
  filter(runtype == 1) %>%          # remove rows with runtype == 0 (BBS category for unacceptable data quality)
  arrange(aou, route_uniq, year)    # arrange rows by route then year

# write to file
save(bbs_full_nozeros, file = 'data/bbs-full-nozeros.RData')



##### Filter to species of interest, and add in species counts of zero
# filter dt to species of interest, and remove unnecessary columns
dt_subset <- filter(dt, Aou %in% species$aou) %>% 
  select(route_yr_rpid, Aou, count)

# filter weather-route df to RunType == 1 (remove runs with unacceptable data quality),
#  and remove unnecessary columns
wx_rt_sub <- filter(wx_rt, RunType == 1) %>% 
  select(route_yr_rpid, route_uniq, Year, route_yr, observer, obs_route, firstyr, state_prov, Stratum, BCR)

# generate all combinations of species and route-year-rpid, where RunType == 1
all_rows <- expand.grid(route_yr_rpid = wx_rt_sub$route_yr_rpid,
                        Aou = species$aou,
                        stringsAsFactors = F)

# put together the full dataset for species of interest, including counts of zero
bbs_use <- all_rows %>% 
  as_tibble() %>% 
  left_join(dt_subset, by = c('Aou', 'route_yr_rpid')) %>%   # join to dt_subset
  left_join(wx_rt, by = 'route_yr_rpid') %>%                 # add back columns of wx_rt
  mutate(count = ifelse(is.na(count), 0, count)) %>%         # change NA bird count to 0
  setNames(tolower(names(.))) %>% 
  select(aou, route_uniq, year, rpid, observer, firstyr, obs_route,
         runtype, count, state_prov, stratum, bcr) %>%   # rearrange columns
  arrange(aou, route_uniq, year)      # arrange rows by route then year
  
# write to file
save(bbs_use, file = 'data/bbs-use.RData')



# ###### simple quality checks
# head(bbs_use)
# nrow(bbs_use)
# length(unique(bbs_use$aou))
# length(unique(bbs_use$stratum))
# length(unique(bbs_use$bcr))
# length(unique(bbs_use$state_prov))
# all(bbs_use$runtype == 1)
# apply(bbs_use, 2, function(x) any(is.na(x)))
# sort(unique(bbs_use$year))
# hist(log(bbs_use$count))
# length(which(bbs_use$count == 0)) / nrow(bbs_use)
# 
# 
# # plot trends over time
# df_trends <- bbs_use %>% 
#   group_by(year) %>% 
#   summarize(n_route = length(unique(route_uniq)),
#             n_observers = length(unique(observer)),
#             sum_counts = sum(count),
#             mean_counts_per_route = sum_counts / n_route)
# 
# plot(n_route ~ year, df_trends, type = 'l')
# plot(n_observers ~ year, df_trends, type = 'l')
# plot(sum_counts ~ year, df_trends, type = 'l')
# plot(mean_counts_per_route ~ year, df_trends, type = 'l')
# 
# 
# # summarize by observer
# df_observer <- bbs_use %>%
#   group_by(observer) %>%
#   summarize_all(function(x) length(unique(x)))
# 
# df_observer_tables <- apply(df_observer[,-1], 2, table)
# 
# df_observer_tables$bcr
# df_observer_tables$stratum
# df_observer_tables$state_prov
# df_observer_tables$year
# df_observer_tables$aou    # should be 28, because includes 0 counts
# 
# 
# # summarize by route
# df_route <- bbs_use %>%
#   group_by(route_uniq) %>%
#   summarize_all(function(x) length(unique(x)))
# 
# df_route_tables <- apply(df_route[,-1], 2, table)
# 
# df_route_tables$bcr         # should all be 1
# df_route_tables$stratum     # should all be 1
# df_route_tables$state_prov  # should all be 1
# df_route_tables$year
# df_route_tables$observer
# 
# 
# # summarize by species (based on counts > 0)
# df_species <- bbs_use %>%
#   filter(count > 0) %>% 
#   group_by(aou) %>%
#   summarize_all(function(x) length(unique(x)))
# 
# df_species_tables <- apply(df_species[,-1], 2, table)
# 
# df_species_tables$bcr
# df_species_tables$stratum
# df_species_tables$state_prov
# df_species_tables$year
# df_species_tables$observer
# 
# 
# # summarize by bcr
# df_bcr <- bbs_use %>%
#   group_by(bcr) %>%
#   summarize_all(function(x) length(unique(x)))
# 
# df_bcr_tables <- apply(df_bcr[,-1], 2, table)
# 
# df_bcr_tables$stratum
# df_bcr_tables$state_prov
# df_bcr_tables$year
# df_bcr_tables$observer
# 
# 
# # check routes where more than one survey in given year
# df_dup_survey <- bbs_use %>% 
#   group_by(route_uniq, year) %>% 
#   summarize(n = n()) %>% 
#   ungroup() %>% 
#   filter(n > 28) %>% 
#   as.data.frame()
# 
# # for given route-year-species, each row should have unique rpid
# i <- sample(1:nrow(df_dup_survey), 1)
# filter(bbs_use, route_uniq == df_dup_survey$route_uniq[i], year == df_dup_survey$year[i], aou == sample(aou, 1))
# 
# 
# # check routes with very high counts
# tail(sort(bbs_use$count), 20)
# dt_raw <- read_csv('data/bbs-allstates-2015.txt')
# filter(dt_raw, SpeciesTotal > 3000) %>% as.data.frame()  # lots of 999 entries


