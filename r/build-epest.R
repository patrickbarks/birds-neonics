
# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(sp)

# set directory
setwd('~/birds-neonics/')

# load polygons for counties within contiguous usa, aea projection
load('data/shape-counties-aea.RData')

# read epest dat
path_epest_2009 <- 'data/epest_raw/1992_2009/'
path_epest_2012 <- 'data/epest_raw/2008_2012/'

files_epest_2009 <- list.files(path_epest_2009)
files_epest_2012 <- list.files(path_epest_2012)

ReadFiles <- function(files, path) read_tsv(paste(path, files, sep = ''))

list_epest_2009 <- lapply(files_epest_2009, ReadFiles, path = path_epest_2009)
list_epest_2012 <- lapply(files_epest_2012, ReadFiles, path = path_epest_2012)

# format county fips codes (unique federal identifiers)
epest_2009 <- do.call(rbind.data.frame, list_epest_2009) %>% 
  mutate(STATE_FIPS_CODE = as.numeric(STATE_FIPS_CODE)) %>% 
  mutate(STATE_FIPS_CODE = formatC(STATE_FIPS_CODE, width = 2, flag = '0')) %>% 
  unite(FIPS, STATE_FIPS_CODE, COUNTY_FIPS_CODE, sep = '')    # merge state and county fips

epest_2012 <- do.call(rbind.data.frame, list_epest_2012) %>% 
  mutate(STATE_FIPS_CODE = as.numeric(STATE_FIPS_CODE)) %>% 
  mutate(STATE_FIPS_CODE = formatC(STATE_FIPS_CODE, width = 2, flag = '0')) %>% 
  unite(FIPS, STATE_FIPS_CODE, COUNTY_FIPS_CODE, sep = '')    # merge state and county fips

# any counties in epest not also in counties shapefile from census.gov?
setdiff(epest_2009$FIPS, counties_aea@data$fips)
setdiff(epest_2012$FIPS, counties_aea@data$fips)

# correct Dade County (name/fips change in 1997)
# https://www.census.gov/geo/reference/county-changes.html
epest_2009$FIPS[epest_2009$FIPS == 12025] <- 12086

# merge epest_2009 and epest_2012 (take 2008-2009 from 2012 file)
epest_merge <- epest_2009 %>% filter(YEAR < 2008) %>%
  rbind(epest_2012) %>% setNames(tolower(names(.)))

# write epest_merge to file
write.csv(epest_merge, 'data/epest-merge.csv', row.names = F)




###################### Get scaled neonic use by county, county-year, etc.

# read merged epest data
epest_merge <- read_csv('data/epest-merge.csv')

# list of neonic compounds
neonics <- c('ACETAMIPRID', 'CLOTHIANIDIN', 'DINOTEFURAN', 'IMIDACLOPRID', 'THIACLOPRID', 'THIAMETHOXAM')

# create df with every unique combination of county, compound, and year
fips_neonic_year_base <- expand.grid(fips = counties_aea@data$fips,
                                     compound = neonics,
                                     year = sort(unique(epest_merge$year)),
                                     stringsAsFactors = F)

# merge epest to county data; scale pesticide use by county area
epest_merge_scaled <- fips_neonic_year_base %>% 
  left_join(epest_merge, by = c('fips', 'compound', 'year')) %>%
  left_join(counties_aea@data, by = 'fips') %>% 
  mutate(kg = ifelse(is.na(kg), 0, kg)) %>% 
  mutate(kg_scaled = kg / (CENSUSAREA / 0.3861019)) %>%  # sq mile to sq km conversion
  as_tibble()

# total neonic use in each county and year, scaled by county area
neonic_county_year <- epest_merge_scaled %>% 
  filter(compound %in% neonics) %>% 
  group_by(fips, year) %>% 
  summarize(neonics = sum(kg_scaled))

write.csv(neonic_county_year, 'data/neonic-county-year.csv', row.names = F)

# average annual total neonic use in each county over the period 2005-2012, scaled by county area
neonic_county_mean <- epest_merge_scaled %>% 
  filter(year >= 2005) %>% group_by(fips) %>% 
  summarize(neonic_mean = mean(kg_scaled))

write.csv(neonic_county_mean, 'data/neonic-county-mean-2005-2012.csv', row.names = F)

# average annual imidacloprid use in each county over the period 2005-2012, scaled by county area
imid_county_mean <- epest_merge_scaled %>% 
  filter(compound == 'IMIDACLOPRID', year >= 2005) %>% group_by(fips) %>% 
  summarize(imid_mean = mean(kg_scaled))

write.csv(imid_county_mean, 'data/imidacloprid-county-mean-2005-2012.csv', row.names = F)

