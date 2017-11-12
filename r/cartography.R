

# load packages
library(dplyr)
library(ggplot2)
library(raster)
library(rgdal)
library(rgeos)
library(rgdal)
library(broom)
library(tidyr)
library(readr)
library(tibble)


# setwd
setwd('~/desktop/bbs/')


# map projections
proj_base <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj_aea <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"  # same as nlcd




############### Create shapefile of American counties within contiguous US

# read counties shapefile
# from https://www.census.gov/geo/maps-data/data/cbf/cbf_counties.html
counties_shp_raw <- readOGR('data/gz_2010_us_050_00_500k', 'gz_2010_us_050_00_500k')

# remove Alaska (02), Hawaii (15), Puerto Rico (72)
counties_shp <- subset(counties_shp_raw, !STATE %in% c('02', '15', '72'))

# add fips code (unique federal identifier for counties)
counties_shp@data$fips <- with(counties_shp@data, paste(STATE, COUNTY, sep = ''))

# re-project counties to AEA for plotting, and remove columns not needed
counties_aea <- counties_shp %>%
  spTransform(CRS(proj_aea)) %>% 
  subset(select = c(STATE, COUNTY, NAME, CENSUSAREA))

# create row and fips ids
counties_aea@data$id <- rownames(counties_aea@data)
counties_aea@data$fips <- with(counties_aea@data, paste(STATE, COUNTY, sep = ''))

# save contiguous counties AEA shapefile
save(counties_aea, file = 'data/shape-counties-aea.RData')




############### Which county does each BBS route fall within, based on coordinates of starting point?

# read bbs routes data
rt <- read_csv('data/bbs_raw_2015/routes.csv') %>% 
  subset(countrynum == '840' & statenum != '03') %>%            # subset to contiguous usa
  mutate(statenum = formatC(statenum, width = 2, flag = '0'),   # format locations vars
         Route = formatC(Route, width = 3, flag = '0'), 
         route_uniq = paste('C', countrynum, 'S', statenum, 'R', Route, sep = '')) # add unique route id

# convert bbs route data to spatial dataframe
rt_sp <- rt
coordinates(rt_sp) <- c('Longitude', 'Latitude')
proj4string(rt_sp) <- proj4string(counties_shp)

# find fips (i.e. county) that each route falls within
rt_sp@data$fips <- over(rt_sp, counties_shp)$fips

# table of number of counties containing X routes (note majority of counties have zero bbs routes)
table(table(rt_sp@data$fips))

# find routes not within any county boundaries
routes_outside <- rt_sp@data$route_uniq[which(is.na(rt_sp@data$fips))]

# Function to find county nearest to route
FindCounty <- function(focal_rt, routes, counties) {
  pDist <- vector()
  route <- subset(routes, route_uniq == focal_rt)
  
  for(j in 1:nrow(counties_shp))
    pDist <- suppressWarnings(append(pDist, gDistance(route, counties[j,])))
  
  return(counties@data$fips[which.min(pDist)])
}

# find counties nearest to the routes that are not within any county polygons
fips_outside <- sapply(routes_outside, FindCounty, routes = rt_sp, counties = counties_shp)
route_fips_outside <- tibble(route_uniq = routes_outside, fips_outside = fips_outside)

# arrange route-county assignments into data frame
route_fips_df <- tibble(route_uniq = rt_sp@data$route_uniq,
                        county_fips = rt_sp@data$fips) %>% 
  left_join(route_fips_outside, by = 'route_uniq') %>% 
  mutate(fips = ifelse(is.na(county_fips), fips_outside, county_fips)) %>% 
  dplyr::select(route_uniq, fips) %>% 
  left_join(dplyr::select(rt, route_uniq, lat = Latitude, lon = Longitude), by = 'route_uniq')

# write to file
write.csv(route_fips_df, 'data/route-county-start.csv', row.names = F)




############### Which counties does each BBS route overlap with, based on full route shapefiles?

# load counties shapefile
load('data/shape-counties-aea.RData')
route_fips_df <- read_csv('data/route-county-start.csv')

# read shapfiles for bbs routes (missing shapes for ~200/4000 routes)
#  and transform to same projection as counties
rt_shp <- readOGR('data/bbs_raw_2015/bbsrte_2012_alb', 'bbsrte_2012_alb') %>% 
  spTransform(CRS(proj4string(counties_shp)))

# format route id
rt_shp@data$rteno <- formatC(rt_shp@data$rteno, width = 5, flag = '0')
rt_shp@data$route_uniq <- paste0('C840', 'S', substr(rt_shp@data$rteno, 1, 2), 'R', substr(rt_shp@data$rteno, 3, 5))

# remove alaska
rt_shp@data$alaska <- F
rt_shp@data$alaska[grep('^03', rt_shp@data$rteno)] <- T
rt_shp <- subset(rt_shp, alaska == F)

# check plot
plot(rt_shp)

# merge so one polygon per route
rt_shp_merge <- gLineMerge(rt_shp, id = rt_shp@data$route_uniq)

# calculate route centroids
rt_centroids <- gCentroid(rt_shp_merge, byid = T)

# arrange route centroids in df
rt_centroids_df <- coordinates(rt_centroids) %>% 
  as.data.frame() %>% 
  setNames(c('lon', 'lat')) %>% 
  dplyr::select(lat, lon) %>% 
  rownames_to_column('route_uniq') %>% 
  as_tibble()

# find counties with which each route overlaps
route_county <- over(rt_shp, counties_shp, returnList = T)

route_county_shape_df <- do.call(rbind.data.frame, route_county) %>% as_tibble() %>% 
  mutate(route_uniq = rep(rt_shp@data$route_uniq, sapply(route_county, nrow))) %>% 
  unique() %>% 
  dplyr::select(route_uniq, fips) %>% 
  left_join(rt_centroids_df, by = 'route_uniq') %>% 
  mutate(type = 'centroid')

# merge route start coords and route shape centroids into one df
route_county_full <- route_fips_df %>% 
  filter(!route_uniq %in% route_county_shape_df$route_uniq) %>%
  mutate(type = 'start') %>% 
  rbind.data.frame(route_county_shape_df) %>% 
  arrange(route_uniq)

# write to file
write.csv(route_county_full, 'data/route-county-full.csv', row.names = F)

# what proportion of routes span X number of counties?
route_county_tab <- route_county_full %>%
  group_by(route_uniq) %>%
  summarize(n = length(unique(fips)))

table(route_county_tab$n) / nrow(route_county_tab)




############### Create grid of quadrats overtop contiguous US, and assign routes to quadrats

# quadrat resolution (km) (150^2 = 22,500 sq. km)
quad_res <- 150

# extent of contiguous USA (conservative estimated with Google Maps)
clip <- as(extent(-130, -60, 15, 55), 'SpatialPolygons')
proj4string(clip) <- CRS(proj_base)

# transform clip to aea projection
clip_aea <- spTransform(clip, CRS(proj_aea))

# transform rt to aea projection
rt_aea <- spTransform(rt_sp, CRS(proj_aea))
rt_aea_coords <- as.data.frame(coordinates(rt_aea))

# raster prep
base_ras <- raster(extent(clip_aea))
res(base_ras) <- quad_res
projection(base_ras) <- proj_aea
base_ras[] <- 1:ncell(base_ras)

# get raster cell for every bbs route in contiguous usa
route_cells <- cellFromXY(base_ras, coordinates(rt_aea))
which(is.na(route_cells))   # any bbs routes not within raster?

# tabulate bbs routes per cell
route_cells_tab <- tabulate(route_cells, ncell(base_ras))
route_cells_tab <- ifelse(route_cells_tab == 0, NA, 1)

# make quadrat raster (1 if contains 1 or more BBS routes, else NA)
quad_ras <- base_ras
quad_ras[] <- route_cells_tab

# convert raster to polygons
quad_poly <- rasterToPolygons(quad_ras)
quad_poly@data$id <- rownames(quad_poly@data)

# plot quadrat polygons
plot(quad_poly)

# find quadrat (polygon id) associated with each route
route_quads <- as.numeric(over(rt_aea, quad_poly)$id)
route_quads_df <- tibble(route_uniq = rt_aea$route_uniq, quad_poly_id = route_quads)

# tabulate routes per quadrat polygon id
poly_route_tab <- tabulate(route_quads, nrow(quad_poly@data))
poly_route_df <- tibble(id = rownames(quad_poly@data), n_routes = poly_route_tab)

# tidy and join quadrat polygons
quad_poly_tidy <- tidy(quad_poly) %>% left_join(poly_route_df, by = 'id')

# plot quadrats and routes (fill with n_routes)
quad_poly_tidy$col <- cut(quad_poly_tidy$n_routes,
                          breaks = c(0, 1, 5, 10, 20, Inf),
                          labels = c('1', '2-5', '6-10', '11-20', '>20'))

ggplot() +
  geom_polygon(data = quad_poly_tidy, aes(long, lat, group = id, fill = col)) +
  geom_path(data = quad_poly_tidy, aes(long, lat, group = id), size = 0.1, col = 'grey50') +
  geom_point(data = rt_aea_coords, aes(x = Longitude, y = Latitude), size = 1, alpha = 0.2) +
  scale_fill_brewer(name = 'Routes per\nquadrat', palette = 'Spectral')

# save quadrat polygons
save(quad_poly, file = 'data/shape-quadrats-aea.RData')

# save route-quadrat dataframe (quadrat associated with each route)
write.csv(route_quads_df, 'data/route-quadrat.csv', row.names = F)




############### Create outline of Contiguous USA for plotting

# shapefile for outline of USA
usa_shp <- raster::getData('GADM', path = 'data/', country = 'USA', level = 0)

# create bbox around contiguous usa (i.e. excluding Alaska, etc.)
clip <- as(extent(-130, -60, 15, 55), 'SpatialPolygons')
proj4string(clip) <- CRS(proj4string(usa_shp))

# clip usa_shp to contiguous usa, and re-project to aea
usa_shp_contig <- gIntersection(usa_shp, clip, byid = T)  # takes ~ 2 mins
usa_shp_contig_aea <- spTransform(usa_shp_contig, CRS(proj_aea))

# simplify shapefile for plotting
usa_shp_contig_aea_simple <- gSimplify(usa_shp_contig_aea, tol = 5)

# plot
plot(usa_shp_contig_aea_simple)

# save to file
save(usa_shp_contig_aea_simple, file = 'data/shape-usa-outline-aea.RData')




############### Create shapefile of Bird Conservation Regions intersecting with contiguous US

# read bcr shapefile
bcr_shp <- readOGR('data/bbs_raw_2015/bcr_shp', 'BCR')

# create row ids
bcr_shp$rowid <- rownames(bcr_shp@data)

# reproject to aea
bcr_shp_aea <- spTransform(bcr_shp, CRS(proj_aea))

# clip to contiguous usa
bcr_shp_contig_aea <- gIntersection(bcr_shp_aea, usa_shp_contig_aea_simple, byid = T)

# tidy
bcr_tidy <- bcr_shp_contig_aea %>% tidy() %>% 
  mutate(rowid = sapply(as.character(group), function(x) strsplit(x, ' ')[[1]][1])) %>% 
  left_join(bcr_shp@data, by = 'rowid') %>% 
  filter(BCRNumber != 39) %>%  # tiny part of Sierras de Baja overlaps with US (contains no routes)
  droplevels()

# bcr centroids
bcr_centroids <- coordinates(bcr_shp_contig_aea) %>%
  as.data.frame() %>%
  rename(x = V1, y = V2) %>% 
  mutate(rowid = sapply(rownames(.), function(x) strsplit(x, ' ')[[1]][1])) %>% 
  left_join(bcr_shp@data, by = 'rowid') %>% 
  filter(BCRNumber != 39) %>% 
  droplevels()

# plot
ggplot(bcr_tidy) +
  geom_path(aes(long, lat, group = group), size = 0.3, alpha = 0.3) +
  geom_text(data = bcr_centroids, aes(x, y, label = BCRNumber), size = 5, color = 'darkblue') +
  theme_bw()

# write bcr_tidy to file
save(bcr_tidy, file = 'data/shape-bcr-aea-tidy.RData')




############### Get proportion cropland surrounding each route

# contiguous us landcover image
nlcd <- raster('data/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img')

# read centroids for bbs routes (for ~200/4000 routes we have only coordinates of starting point)
route_county_full <- read_csv('data/route-county-full.csv')

# convert bbs route data to spatial dataframe
coordinates(route_county_full) <- c('lon', 'lat')
proj4string(route_county_full) <- proj4string(counties_shp)

# transform to same projection as nlcd
route_county_full_aea <- spTransform(route_county_full, CRS(proj4string(nlcd)))

# arrange route centroids into dataframe
route_coords_full <- route_county_full_aea@data %>% 
  mutate(lon = route_county_full_aea@coords[,1],
         lat = route_county_full_aea@coords[,2]) %>% 
  dplyr::select(-fips) %>% unique()

# plot route centroids over landcover raster
plot(nlcd)
plot(route_county_full_aea, add = T, col = 'red')

# Function to get percentage pasture and cropland around each route
GetAgRoute <- function(x, y, route_uniq, buffer_rad = 20000, temp1, temp2) {
  
  focal_crop <- crop(
    nlcd,
    extent(x - buffer_rad,
           x + buffer_rad,
           y - buffer_rad,
           y + buffer_rad),
    filename = temp1
  )
  
  base_raster <- focal_crop
  base_raster[] <- NA
  base_raster[cellFromXY(focal_crop, t(c(x, y)))] <- 1
  
  buff <- buffer(base_raster, width = buffer_rad, filename = temp2)
  focal_crop[which(is.na(buff[]))] <- NA
  
  focal_crop <- getValues(focal_crop)
  
  cells_crop <- length(which(focal_crop %in% c(61, 82:84)))
  cells_past <- length(which(focal_crop == 81))
  cells_tot <- length(which(!focal_crop %in% c(NA, 0, 11)))
  pct_crop <- round(cells_crop / cells_tot * 100, 3)
  pct_past <- round(cells_past / cells_tot * 100, 3)
  
  return(data.frame(pct_crop, pct_past))
}

# create temporary directories to write raster subsets for each route (otherwise will run out of memory)
rasterOptions(overwrite = T)
temp_base <- tmpDir()
temp1 <- file.path(temp_base, 'garbage1')
temp2 <- file.path(temp_base, 'garbage2')

# get percentage cropland arround each route
landcover <- route_coords_full %>%
  group_by(route_uniq) %>% 
  do(GetAgRoute(x = .$lon, y = .$lat, route_uniq = .$route_uniq, buffer_rad = 20000, temp1, temp2)) %>%
  ungroup()

# write to file
write.csv(landcover, 'data/landcover-route-2011-20km.csv', row.names = F)


