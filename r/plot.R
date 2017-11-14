
# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(broom)
library(tibble)
library(sp)
library(ggplot2)
library(grid)
library(gridExtra)

# setwd
setwd('~/birds-neonics/')

# read species list and epest data
species_df <- read_csv('data/species-list.csv')
epest <- read_csv('data/epest-merge.csv')



######### Figure 1
# population abundance over the period 1985-2014, by species

# organize data
summary_trend_allyrs <- read_csv('analysis/post-summary-spp-temporal-strata-trend-allyears.csv')

summary_abundance_85_14 <- read_csv('analysis/post-summary-spp-temporal-strata-abundance.csv') %>% 
  left_join(species_df, by = 'aou')

beta_summary_85_14 <- summary_trend_allyrs %>%
  mutate(sig_change = ifelse(beta_upp90 < 0, '↓', ''),
         sig_change = ifelse(beta_low90 > 0, '↑', sig_change)) %>% 
  as.data.frame()

beta_summary_85_14 %>%
  summarize(sig_down = length(which(beta_upp90 < 0)),
            sig_up = length(which(beta_low90 > 0)))

df_plot_abundance <- summary_abundance_85_14 %>%
  left_join(beta_summary_85_14, by = 'aou') %>% 
  mutate(species_short = reorder(factor(species_short), aou, min))

# create dummy data to add extra space at top of individual panels
dummy_plot <- df_plot_abundance %>% group_by(aou) %>% 
  summarize(plot_min = min(comp_low90),
            plot_max = max(comp_upp90),
            log_min = log(plot_min),
            log_max = log(plot_max),
            log_buffer = (log_max - log_min) * 0.27,
            new_max = exp(log_max + log_buffer),
            x = min(year)) %>% 
  left_join(dplyr::select(species_df, aou, species_short), by = 'aou') %>% 
  mutate(species_short = reorder(factor(species_short), aou, min))

# plotting functions to format axes breaks and labels
BreakFn <- function(limits) {
  log_low <- log(min(limits))
  log_upp <- log(max(limits))
  log_buffer <- 0.08 * (log_upp - log_low)
  new_low <- exp(log_low + log_buffer)
  new_upp <- exp(log_upp - log_buffer)
  out <- c((1:9)/10, (1:9)/1, (1:9)/0.1, (1:10)/0.01)
  return(out[out >= new_low & out <= new_upp])
}

LabelFn <- function(breaks) {
  n_10s <- length(which(breaks %in% c(0.01, 0.1, 1, 10, 100, 1000)))
  
  if(n_10s >= 2) {
    return(ifelse(breaks %in% c(0.01, 0.1, 1, 10, 100), as.character(breaks), ''))
  } else if(any(breaks %in% c(1, 10, 100))) {
    i <- which(breaks %in% c(1, 10, 100))
    j <- ifelse(i == length(breaks), 1, length(breaks))
    out <- as.character(breaks)
    out[-c(i, j)] <- ''
    return(out)
  } else {
    out <- as.character(breaks)
    out[-c(1, length(breaks))] <- ''
    return(out)
  }
}

# create full plot
fig_1 <- ggplot(df_plot_abundance) +
  geom_ribbon(aes(x = year, ymin = comp_low90, ymax = comp_upp90), fill = 'blue', alpha = 0.3) +
  geom_line(aes(x = year, y = comp_med), size = 0.4, alpha = 0.9, col = 'grey30') +
  geom_blank(data = dummy_plot, aes(x, new_max)) +
  scale_y_log10(breaks = BreakFn, labels = LabelFn) +
  geom_text(aes(x = Inf, y = Inf, label = sig_change), hjust = 1.4, vjust = 1.3, size = 4, fontface = 'bold') +
  geom_text(aes(x = 1985, y = Inf, label = species_short), hjust = 0, vjust = 1.6, size = 2.6) +
  facet_wrap(~ species_short, nrow = 7, scales = 'free_y') +
  xlab('Year') +
  ylab('Index of abundance') +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        text = element_text(size = 12),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.2, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm'))) +
  theme(axis.text.x = element_text(size = 9, color = 'grey40'),
        axis.text.y = element_text(size = 7.5, angle = 90, hjust = 0.5, color = 'grey40'),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.1, 'cm'),
        strip.text = element_text(size = 7),
        panel.spacing.x = unit(0.1, 'cm'),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, 'cm'),
        strip.background = element_blank(),
        strip.text.x = element_blank())

# check plot
dev.off()
quartz(height = 6, width = 6)
print(fig_1)

# write to file
ggsave('figures/fig-1.png', fig_1, height = 6, width = 6, units = 'in')





######### Figure 2
# temporal relationship between neonic use and 7-year aggregate population trends

# organize data
load('analysis/post-summary-agg-temporal-strata-trend-7years.RData')   # summary_agg_trend_7years

df_plot_7years <- summary_agg_trend_7years %>%
  mutate(year_med = year_start + 4)

neonics <- c('ACETAMIPRID', 'CLOTHIANIDIN', 'DINOTEFURAN', 'IMIDACLOPRID', 'THIACLOPRID', 'THIAMETHOXAM')
neonics_focal <- c('CLOTHIANIDIN', 'IMIDACLOPRID', 'THIAMETHOXAM')
neonics_other <- c('ACETAMIPRID', 'DINOTEFURAN', 'THIACLOPRID')

neonic_summary_focal <- epest %>% 
  filter(compound %in% neonics_focal) %>% 
  group_by(compound, year) %>% 
  summarize(kg_tot = sum(kg)) %>% 
  ungroup()

neonic_summary_total <- epest %>% 
  filter(compound %in% neonics) %>% 
  group_by(year) %>% 
  summarize(kg_tot = sum(kg)) %>% 
  ungroup() %>% 
  mutate(compound = 'TOTAL') %>% 
  dplyr::select(compound, year, kg_tot)

neonic_summary_other <- epest %>% 
  filter(compound %in% neonics_other) %>% 
  group_by(year) %>% 
  summarize(kg_tot = sum(kg)) %>% 
  ungroup() %>% 
  mutate(compound = 'other neonic.') %>% 
  dplyr::select(compound, year, kg_tot)

neonic_full <- rbind.data.frame(neonic_summary_focal, neonic_summary_other) %>% 
  mutate(compound = factor(tolower(compound), levels = c('clothianidin', 'imidacloprid', 'thiamethoxam', 'other neonic.')))

# upper panel, 7-year aggregate population trends from 1985 to 2014
p1 <- ggplot(df_plot_7years) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_linerange(aes(x = year_med, ymin = mu_theta_low90, ymax = mu_theta_upp90), size = 0.9) +
  geom_linerange(aes(x = year_med, ymin = mu_theta_low99, ymax = mu_theta_upp99), size = 0.3) +
  geom_point(aes(x = year_med, y = mu_theta_med)) +
  xlab(NULL) + ylab(expression(paste('Aggregate trend (% ', year^-1, ')'))) +
  scale_x_continuous(limits = c(1988, 2012)) +
  scale_y_continuous(breaks = seq(-4, 2, 2)) +
  annotate('text', x = -Inf, y = Inf, label = 'a', hjust = -0.8, vjust = 1.5, size = 6, fontface = 'bold') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 13.5),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.title.x = element_text(margin = margin(.2, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')))

# lower panel, total neonic use over time
p2 <- ggplot(neonic_full) +
  geom_line(aes(year, kg_tot/1000000, color = compound), size = 2) +
  xlab('Year') + ylab('Neonicotinoid use (million kg)') +
  scale_x_continuous(limits = c(1988, 2012)) +
  annotate('text', x = -Inf, y = Inf, label = 'b', hjust = -0.8, vjust = 1.5, size = 6, fontface = 'bold') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 13.5),
        axis.title = element_text(size = 14),
        axis.title.x = element_text(margin = margin(0.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, 0.3, 0, 0, unit = 'cm')),
        legend.position = c(0.6, 0.8),
        legend.title = element_blank(),
        legend.background = element_blank())

# arrange panels
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g1$widths = g2$widths

fig_2 <- arrangeGrob(g1, g2, nrow = 2, heights = c(0.47, 0.53))

# plot
dev.off()
quartz(height = 6, width = 5)
grid.arrange(fig_2)

# save to file
ggsave('figures/fig-2.png', fig_2, height = 6, width = 5, units = 'in', dpi = 300)





######### Figure 3

# load relevant data
load('analysis/post-summary-agg-temporal-quadrat.RData') # summary_agg_temporal_quadrat
load('data/shape-counties-aea.RData')                    # counties_aea
load('data/shape-quadrats-aea.RData')                    # quad_poly
load('data/shape-usa-outline-aea.RData')                 # usa_shp_contig_aea_simple
load('data/shape-bcr-aea-tidy.RData')                    # bcr_tidy
load('analysis/post-summary-agg-spatial-neonic.RData')   # summary_agg_neonic
bcr_dat <- read_csv('data/bbs_raw_2015/BCR.csv')         

# read county-specific neonic use, and add color scale based on quantiles
neonic_county_plot <- read_csv('data/neonic-county-mean-2005-2012.csv') %>% 
  right_join(data.frame(fips = counties_aea@data$fips, stringsAsFactors = F), by = 'fips') %>% 
  mutate(neonic_mean = ifelse(is.na(neonic_mean), 0, neonic_mean)) %>% 
  mutate(col = cut(neonic_mean, breaks = quantile(neonic_mean, seq(0, 1, 0.05)))) %>% 
  mutate(col_no = as.numeric(col))

# read quadrat-specific population trends, and add color scale based on quantiles
betas_quadrat_plot <- summary_agg_temporal_quadrat %>% 
  mutate(col = cut(theta_med, breaks = quantile(theta_med, c(seq(0, 1, 0.01))))) %>% 
  mutate(col_no = as.numeric(col)) %>% 
  rename(id = quadrat)

# tidy shapefiles for gg-plotting
us_outline_tidy <- tidy(usa_shp_contig_aea_simple)

quadrats_tidy <- tidy(quad_poly) %>%
  left_join(betas_quadrat_plot, by = 'id')

counties_tidy <- tidy(counties_aea) %>% 
  left_join(counties_aea@data, by = 'id') %>% 
  left_join(neonic_county_plot, by = 'fips')

# get quadrat centroids, and check whether quadrat represented in betas_quadrat_plot
quadrat_centroid <- coordinates(quad_poly) %>%
  as.data.frame() %>%
  rename(x = V1, y = V2) %>%
  mutate(id = rownames(quad_poly@data),
         display = ifelse(is.element(id, betas_quadrat_plot$id), FALSE, TRUE))

# create scale bars for neonic use and population trends, based on quantiles
neonic_scale_labels <- neonic_county_plot$neonic_mean %>% 
  quantile(c(0, 0.25, 0.5, 0.75, 1)) %>%
  formatC(digits = 1, format = 'fg')

trend_scale_labels <- betas_quadrat_plot$theta_med %>% 
  quantile(c(0, 0.25, 0.5, 0.75, 1)) %>%
  formatC(digits = 1, format = 'f')

# format neonicotinoid effect, by species and aggregate, over the period 2005-2014
df_neonic_spp <- filter(summary_agg_neonic, year == 2010)$df_spp[[1]] %>% 
  left_join(species_df, by = 'aou') %>% 
  mutate(species_short = reorder(species_short, theta_med, mean))

df_neonic_agg <- filter(summary_agg_neonic, year == 2010)$alpha_df[[1]]

# panel 1
p1 <- ggplot() +
  geom_blank(data = quadrats_tidy, aes(long, lat, group = id)) +
  geom_blank(data = quadrats_tidy, aes(long + 500, lat, group = id)) +
  geom_polygon(data = counties_tidy, aes(long, lat, group = group, fill = col_no*5)) +
  geom_path(data = bcr_tidy, aes(long, lat, group = group), size = 0.2) +
  coord_cartesian(expand = F) +
  scale_fill_gradient(low = 'white',
                      high = '#a50026',
                      na.value = 'white',
                      limits = c(-2, 102),
                      labels = neonic_scale_labels,
                      name = 'Neonic. use\n(kg per km²\nper year)',
                      breaks = c(0, 25, 50, 75, 100)) +
  theme(line = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.2, 0.5), 'lines'),
        legend.position = c(1, 0),
        legend.justification = c(1, 0.105),
        legend.title = element_text(size = 8.1, color = 'grey30'),
        legend.text = element_text(size = 7, color = 'grey30'),
        legend.key.height = unit(0.8, 'lines'),
        legend.key.width = unit(1.1, 'lines'))

# panel 2
p2 <- ggplot() +
  geom_blank(data = quadrats_tidy, aes(long + 500, lat, group = id)) +
  geom_polygon(data = quadrats_tidy, aes(long, lat, group = id, fill = col_no), alpha = 1) +
  geom_path(data = quadrats_tidy, aes(long, lat, group = id), size = 0.1, col = 'white') +
  geom_path(data = us_outline_tidy, aes(long, lat, group = group), size = 0.2, alpha = 0.4) +
  geom_point(data = filter(quadrat_centroid, display == T), aes(x, y), size = 1.1, shape = 4) +
  coord_cartesian(expand = F) +
  scale_fill_gradient(low = '#053061',
                      high = '#ece7f2',
                      na.value = 'white',
                      limits = c(-2, 102),
                      labels = trend_scale_labels,
                      name = 'Aggregate\ntrend (% per\nyear)',
                      breaks = c(0, 25, 50, 75, 100)) +
  theme(line = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.2, 0.5, 0.5, 0.5), 'lines'),
        legend.position = c(1, 0),
        legend.justification = c(1, 0.105),
        legend.title = element_text(size = 8.1, color = 'grey30'),
        legend.text = element_text(size = 7, color = 'grey30'),
        legend.key.height = unit(0.8, 'lines'),
        legend.key.width = unit(1.1, 'lines'))

# panel 3
p3 <- ggplot(df_neonic_spp) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5, alpha = 0.3) +
  geom_point(aes(y = species_short, x = theta_med), col = 'black', shape = 16, size = 1.7, alpha = 1) +
  geom_errorbarh(aes(y = species_short, x = theta_med, xmin = theta_low90, xmax = theta_upp90), height = 0, size = 0.6) +
  geom_errorbarh(aes(y = species_short, x = theta_med, xmin = theta_low99, xmax = theta_upp99), height = 0, size = 0.2) +
  scale_x_continuous(breaks = seq(-3, 2, 1)) +
  xlab(expression(paste('Neonic. effect (% ', year^-1, ' ', neonic^-1, ')'))) +
  ylab(NULL) +
  theme_bw() +
  theme(axis.title = element_text(size = 9.5),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.2, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .13, 0, 0, unit = 'cm'))) +
  theme(axis.text.x = element_text(size = 8.5, color = 'grey30'),
        axis.text.y = element_text(size = 6.8, color = 'grey30'),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.2, 'lines'),
        plot.margin = unit(c(0.5, 0.5, 0.4, 0.5), 'lines'))

# panel 4
p4 <- ggplot(df_neonic_agg, aes(alpha_geom)) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5, alpha = 0.3) +
  geom_density(fill = 'darkred', alpha = 0.5, size = 0.1) +
  scale_x_continuous(breaks = seq(-1.2, 0.6, 0.6)) +
  scale_y_continuous(breaks = c(0, 1)) +
  xlab(expression(paste('Aggregate neonic. effect'))) +
  ylab('Posterior density') +
  theme_bw() +
  theme(axis.title = element_text(size = 9.5),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.1, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .1, 0, 0, unit = 'cm'))) +
  theme(axis.text.x = element_text(size = 8.5, color = 'grey30'),
        axis.text.y = element_text(size = 8.5, color = 'grey30'),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.2, 'lines'),
        plot.margin = unit(c(0.5, 0.5, 1, 0.5), 'lines'))

# arrange panels into grobs
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)

# standardize dimensions for left panels
g1$widths = g2$widths
g1$heights = g2$heights

# grob for right panel
g_right <- arrangeGrob(p3, p4, nrow = 2, heights = c(0.72, 0.28))

# layout for full spatial figure
lay <- rbind(c(1, 1, 3), c(2, 2, 3))

# arrange full plot
fig_3 <- arrangeGrob(g1, g2, g_right, layout_matrix = lay, widths = c(0.33, 0.33, 0.36))

# write to file
ggsave('figures/fig-3.png', fig_3, height = 5.5, width = 8, units = 'in')





################### Fig S1 (example of model output from species-level spatial neonic model)
# use data for the Bobolink, over the period 2005-2014

# load stanfit object
load('stanfit/spatial-neonic/fit-neonic-county-nested4-04940-2005to2014.RData')

# load bbs
load('data/bbs-use.RData')

# read route-county data
route_counties <- read_csv('data/route-county-start.csv')

# read neonic data
neonics_counties <- read_csv('data/neonic-county-mean-2005-2012.csv') %>% 
  mutate(neonic_cubert = neonic_mean^(1/3)) %>% 
  mutate(neonic_scale = (neonic_cubert - mean(neonic_cubert)) / sd(neonic_cubert))

# read bcr data
bcr_dat <- read_csv('data/bbs_raw_2015/BCR.csv') %>%
  mutate(bcr = as.character(BCR)) %>%
  dplyr::select(bcr, bcr_short)

# filter bbs data for species of interest
bbs_filter <- bbs_use %>% 
  filter(aou == '04940', runtype == 1, year %in% 2005:2014) %>% 
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

# arrange fips (counties) by bcr and neonic use
fips_by_bcr <- bbs_filter %>%
  dplyr::select(fips_bcr, bcr, fips) %>% unique() %>% 
  left_join(neonics_counties, by = 'fips') %>% 
  group_by(bcr) %>% 
  mutate(neonic_mean_cent = neonic_scale - mean(neonic_scale)) %>% ungroup() %>% 
  mutate(fips_bcr_no = as.numeric(as.factor(fips_bcr)), bcr_no = as.numeric(as.factor(bcr))) %>% 
  arrange(fips_bcr_no)

# extract posterior samples from stanfit object
beta_bcr <- rstan::extract(fit, 'beta_bcr')$beta_bcr
dev_beta_county <- rstan::extract(fit, 'dev_beta_strata')$dev_beta_strata
gamma <- rstan::extract(fit, 'g1')$g1

# tidy posterior samples for region-specific slopes
df_beta_bcr <- beta_bcr %>%
  as_tibble() %>%
  setNames(sort(unique(fips_by_bcr$bcr))) %>%
  rownames_to_column(var = 'rep') %>%
  gather(bcr, beta_bcr, -rep)

# tidy posterior samples for county-specific slope deviations
df_county <- dev_beta_county %>%
  as_tibble() %>%
  setNames(sort(unique(fips_by_bcr$fips_bcr))) %>%
  rownames_to_column(var = 'rep') %>%
  gather(fips_bcr, dev_county, -rep) %>%
  left_join(dplyr::select(fips_by_bcr, fips_bcr, bcr), by = 'fips_bcr') %>%
  mutate(bcr = as.character(bcr)) %>%
  left_join(df_beta_bcr, by = c('rep', 'bcr')) %>%
  mutate(county_trend = 100 * (exp(beta_bcr + dev_county) - 1))  %>%
  group_by(fips_bcr) %>%
  summarize(med = quantile(county_trend, 0.500),
            low = quantile(county_trend, 0.025),
            upp = quantile(county_trend, 0.975)) %>%
  left_join(dplyr::select(fips_by_bcr, fips_bcr, bcr, neonic_mean_cent), by = 'fips_bcr') %>%
  mutate(bcr = as.character(bcr)) %>%
  left_join(bcr_dat, by = 'bcr')

# get range of county-specific neonic values for each BCR
df_county_range <- df_county %>%
  group_by(bcr) %>%
  summarize(neonic_min = min(neonic_mean_cent),
            neonic_max = max(neonic_mean_cent))

# function to produce best-fit lines
PredFn <- function(xmin, xmax, alpha, beta) {
  x_pred <- seq(xmin, xmax, length.out = 20)
  pred <- 100 * (exp(alpha + beta * x_pred) - 1)
  data.frame(x_pred, pred)
}

# get best-fit line and 95% CI for each BCR
# reflects the relationship between county-specific population trends and county-specific neonic use, within BCR
df_gamma <- gamma %>%
  as_tibble() %>%
  setNames(sort(unique(fips_by_bcr$bcr))) %>%
  rownames_to_column(var = 'rep') %>%
  gather(bcr, gamma, -rep) %>%
  left_join(df_county_range, by = 'bcr') %>%
  left_join(df_beta_bcr, by = c('bcr', 'rep')) %>%
  group_by(bcr, rep) %>%
  do(PredFn(xmin = .$neonic_min, xmax = .$neonic_max, alpha = .$beta_bcr, beta = .$gamma)) %>%
  ungroup() %>%
  group_by(bcr, x_pred) %>%
  summarize(med = quantile(pred, 0.500),
            low = quantile(pred, 0.025),
            upp = quantile(pred, 0.975)) %>%
  ungroup() %>%
  left_join(bcr_dat, by = 'bcr')

# create plot
fig_s1 <- ggplot() +
  geom_linerange(data = df_county, aes(x = neonic_mean_cent, ymin = low, ymax = upp), col = 'darkblue', alpha = 0.6) +
  geom_ribbon(data = df_gamma, aes(x = x_pred, ymin = low, ymax = upp), alpha = 0.3) +
  geom_line(data = df_gamma, aes(x_pred, med), alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.3) +
  facet_wrap(~ bcr_short, scale = 'free_y', nrow = 5) +
  xlab('Mean neonicotinoid use, 2005-2012 (z-standardized nationally; centered within regions)') +
  ylab(expression(paste('Bobolink population trend, 2005-2014 (% ', year^-1, ')'))) +
  theme_bw() +
  theme(axis.title = element_text(size = 13),
        text = element_text(size = 13),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.2, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .2, 0, 0, unit = 'cm'))) +
  theme(axis.text.x = element_text(size = 10.5, color = 'grey30'),
        axis.text.y = element_text(size = 9.5, color = 'grey30'),
        axis.ticks = element_line(size = 0.2),
        axis.ticks.length = unit(0.2, 'lines'),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'lines'))

# view plot
dev.off()
quartz(height = 6, width = 9.5)
print(fig_s1)

# write to file
ggsave('figures/fig-s1.png', fig_s1, height = 6, width = 9.5, units = 'in', dpi = 300)





######### Figure S2
# aggregate population trends over the period 1985-2014, by species

# organize data
load('analysis/post-summary-agg-temporal-strata-trend-allyears.RData') # summary_agg_trend_allyears

df_plot_allyears_spp <- summary_agg_trend_allyears$df_spp$df_spp %>% 
  left_join(species_df, by = 'aou') %>% 
  mutate(species_short = reorder(species_short, theta_med, mean))

df_plot_allyears_agg <- summary_agg_trend_allyears$df_mu_theta$df_mu_theta

# upper panel, geometric population trend over 1985-2014, by species
p1 <- ggplot(df_plot_allyears_spp) +
  geom_vline(xintercept = 0, alpha = 0.4, linetype = 2) +
  geom_point(aes(x = theta_med, y = species_short)) +
  geom_errorbarh(aes(x = theta_med, xmin = theta_low90, xmax = theta_upp90, y = species_short), height = 0, size = 0.9) +
  geom_errorbarh(aes(x = theta_med, xmin = theta_low99, xmax = theta_upp99, y = species_short), height = 0, size = 0.3) +
  xlab(expression(paste('Population trend, 1985-2014 (% ', year^-1, ')'))) + ylab(NULL) +
  scale_x_continuous(breaks = seq(-6, 2, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10.5),
        plot.margin = unit(c(0.3, 0.5, 0.5, 1), 'lines'))

# lower panel, density plot of aggregate geometric population trend over 1985-2014
p2 <- ggplot(df_plot_allyears_agg, aes(mu_theta)) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.5, alpha = 0.4) +
  geom_density(fill = 'darkred', alpha = 0.5, size = 0.3) +
  xlab(expression(paste('Aggregate population trend, 1985-2014 (% ', year^-1, ')'))) +
  ylab('Posterior density') +
  scale_x_continuous(breaks = seq(0, -2, -1)) +
  scale_y_continuous(breaks = c(0, 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 10.5),
        plot.margin = unit(c(0.5, 0.5, 0.2, 1), 'lines'))

# arrange panels
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

fig_s2 <- arrangeGrob(g1, g2, nrow = 2, heights = c(0.74, 0.26))

# plot
dev.off()
quartz(height = 6.5, width = 3.75)
grid.arrange(fig_s2)

# save to file
ggsave('figures/fig-s2.png', fig_s2, height = 6.5, width = 3.75, units = 'in', dpi = 300)





######### Other spatial plots (neonic, cropland, imidacloprid, neonic mismatch)

# load relevant data
load('analysis/post-summary-agg-spatial-neonic.RData')           # summary_agg_neonic
load('analysis/post-summary-agg-spatial-cropland.RData')         # summary_agg_cropland
load('analysis/post-summary-agg-spatial-imidacloprid.RData')     # summary_agg_imidacloprid
load('analysis/post-summary-agg-spatial-neonic-mismatch.RData')  # summary_agg_neonic_mismatch
bcr_dat <- read_csv('data/bbs_raw_2015/BCR.csv') %>% mutate(bcr = as.character(BCR))

# ggplot theme
tt <-  theme_bw() +
  theme(axis.title = element_text(size = 13),
        panel.grid = element_blank(),
        axis.title.x = element_text(margin = margin(.2, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm'))) +
  theme(axis.text.x = element_text(size = 10, color = 'grey30'),
        axis.text.y = element_text(size = 9, color = 'grey30'),
        axis.ticks = element_line(size = 0.2),
        plot.margin = unit(c(0.5, 0.5, 0.5, 1), 'lines'))

# functionalize
SpatialPlotFn <- function(data, year_index, xlab1, xlab2, breaks1, breaks2, breaks3, geom = T) {
  
  # organize data
  df_spp <- data$df_spp[[year_index]] %>%
    left_join(species_df, by = 'aou') %>%
    mutate(species_short = reorder(species_short, theta_med, mean))
  
  df_bcr <- data$df_bcr[[year_index]] %>%
    left_join(bcr_dat, by = 'bcr') %>%
    mutate(bcr_short = reorder(bcr_short, theta_med, mean))
  
  df_agg <- data$alpha_df[[year_index]]
  
  # if want to depict untransformed data (i.e. exponential rather than geometric)
  if(geom == F) {
    df_spp <- df_spp %>% dplyr::select(-theta_mean, -theta_se, -theta_med, -theta_low90, -theta_upp90, -theta_low99, -theta_upp99) %>% setNames(gsub('_raw', '', names(.)))
    df_bcr <- df_bcr %>% dplyr::select(-theta_mean, -theta_se, -theta_med, -theta_low90, -theta_upp90, -theta_low99, -theta_upp99) %>% setNames(gsub('_raw', '', names(.)))
    df_agg <- df_agg %>% dplyr::select(-alpha_geom) %>% setNames('alpha_geom')
  }
  
  # panel 1, by species
  p1 <- ggplot(df_spp) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.5, alpha = 0.3) +
    geom_point(aes(y = species_short, x = theta_med), col = 'black', shape = 16, size = 2, alpha = 1) +
    geom_errorbarh(aes(y = species_short, x = theta_med, xmin = theta_low90, xmax = theta_upp90), height = 0, size = 0.9) +
    geom_errorbarh(aes(y = species_short, x = theta_med, xmin = theta_low99, xmax = theta_upp99), height = 0, size = 0.3) +
    scale_x_continuous(breaks = breaks1) +
    xlab1 + ylab(NULL) + tt
  
  # panel 2, by bcr
  p2 <- ggplot(df_bcr) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.5, alpha = 0.3) +
    geom_point(aes(y = bcr_short, x = theta_med), col = 'black', shape = 16, size = 2, alpha = 1) +
    geom_errorbarh(aes(y = bcr_short, x = theta_med, xmin = theta_low90, xmax = theta_upp90), height = 0, size = 0.9) +
    geom_errorbarh(aes(y = bcr_short, x = theta_med, xmin = theta_low99, xmax = theta_upp99), height = 0, size = 0.3) +
    scale_x_continuous(breaks = breaks2) +
    xlab1 + ylab(NULL) + tt
  
  # panel 3, density of aggregate effect
  p3 <- ggplot(df_agg, aes(alpha_geom)) +
    geom_vline(xintercept = 0, linetype = 2, size = 0.5, alpha = 0.3) +
    geom_density(fill = 'darkred', alpha = 0.5, size = 0.1) +
    scale_x_continuous(breaks = breaks3) +
    xlab2 + ylab('Posterior density') + tt + theme(plot.margin = unit(c(1.3, 0.5, 0.5, 1), 'lines'))
  
  # spatial payout of panels
  lay <- rbind(c(1, 2), c(3, 3))
  
  # arrange full plot
  plot_full <- arrangeGrob(p1, p2, p3, layout_matrix = lay, heights = c(0.7, 0.3), widths = c(0.47, 0.53))
  
  return(plot_full)
}

# create all possible x-axis labels for spatial effect plots
xlab1_1 <- xlab('Neonicotinoid effect')
xlab1_2 <- xlab('Cropland effect')
xlab1_3 <- xlab('Imidacloprid effect')

xlab2_1 <- xlab(expression(paste('Aggregate neonicotinoid effect (% ', year^-1, ' ', neonic.^-1, ')')))
xlab2_2 <- xlab(expression(paste('Aggregate cropland effect (% ', year^-1, ' ', crop.^-1, ')')))
xlab2_3 <- xlab(expression(paste('Aggregate imidacloprid effect (% ', year^-1, ' ', imid.^-1, ')')))
xlab2_4 <- xlab(expression(paste('Aggregate neonicotinoid effect (individs. ', individ.^-1, ' ', year^-1, ' ', neonic^-1, ')')))
xlab2_5 <- xlab(expression(paste('Aggregate cropland effect (individs. ', individ.^-1, ' ', year^-1, ' ', crop.^-1, ')')))

# generate all spatial effect plots
plot_neonic <- SpatialPlotFn(summary_agg_neonic, 7, xlab1_1, xlab2_1, seq(-2, 2, 1), seq(-4, 2, 2), seq(-1.5, 1, 0.5))
plot_cropland <- SpatialPlotFn(summary_agg_cropland, 1, xlab1_2, xlab2_2, seq(-2, 2, 1), seq(-4, 4, 2), seq(-1, 1, 0.5))
plot_imidacloprid <- SpatialPlotFn(summary_agg_imidacloprid, 1, xlab1_3, xlab2_3, seq(-3, 2, 1), seq(-4, 2, 2), seq(-1.5, 0.5, 0.5))
plot_neonic_mismatch <- SpatialPlotFn(summary_agg_neonic_mismatch, 1, xlab1_1, xlab2_1, seq(-2, 2, 1), seq(-4, 2, 2), seq(-1, 1, 0.5))

plot_neonic_untrans <- SpatialPlotFn(summary_agg_neonic, 7, xlab1_1, xlab2_4, seq(-0.02, 0.02, 0.02), seq(-0.04, 0.02, 0.02), seq(-0.015, 0.01, 0.005), geom = F)
plot_cropland_untrans <- SpatialPlotFn(summary_agg_cropland, 1, xlab1_2, xlab2_5, seq(-0.02, 0.02, 0.01), seq(-0.04, 0.04, 0.02), seq(-0.01, 0.01, 0.005), geom = F)

# check plots
dev.off()
quartz(height = 7, width = 7)
grid.arrange(plot_neonic)

dev.off()
quartz(height = 7, width = 7)
grid.arrange(plot_cropland)

dev.off()
quartz(height = 7, width = 7)
grid.arrange(plot_imidacloprid)

dev.off()
quartz(height = 7, width = 7)
grid.arrange(plot_neonic_mismatch)

dev.off()
quartz(height = 7, width = 7)
grid.arrange(plot_neonic_untrans)

dev.off()
quartz(height = 7, width = 7)
grid.arrange(plot_cropland_untrans)

# write plots to file
ggsave('figures/fig-s3.png', plot_neonic, height = 7, width = 7, units = 'in', dpi = 300)         # Fig S3
ggsave('figures/fig-s4.png', plot_imidacloprid, height = 7, width = 7, units = 'in', dpi = 300)   # Fig S4
ggsave('figures/fig-s8.png', plot_cropland, height = 7, width = 7, units = 'in', dpi = 300)       # Fig S8
ggsave('figures/appendix-spatial-neonic-mismatch.png', plot_neonic_mismatch, height = 7, width = 7, units = 'in', dpi = 300) # Fig A4
ggsave('figures/appendix-untransformed-neonic.png', plot_neonic_untrans, height = 7, width = 7, units = 'in', dpi = 300)     # Fig A1
ggsave('figures/appendix-untransformed-cropland.png', plot_cropland_untrans, height = 7, width = 7, units = 'in', dpi = 300) # Fig A2





######### Figure S5 (aggregate neonicotinoid effect over time)

# load relevant data
load('analysis/post-summary-agg-spatial-neonic.RData')  # summary_agg_neonic

df_agg_neonic_plot <- summary_agg_neonic %>% 
  dplyr::select(year, alpha_df) %>%
  unnest() %>% 
  group_by(year) %>% 
  summarize(med =   quantile(alpha_geom, 0.500),
            low90 = quantile(alpha_geom, 0.050),
            upp90 = quantile(alpha_geom, 0.950),
            low99 = quantile(alpha_geom, 0.005),
            upp99 = quantile(alpha_geom, 0.995))

# create plot
fig_s5 <- ggplot(df_agg_neonic_plot) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4) +
  geom_point(aes(year, med), size = 1.9) +
  geom_linerange(aes(year, ymin = low90, ymax = upp90), size = 1.4) +
  geom_linerange(aes(year, ymin = low99, ymax = upp99), size = 0.4) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_x_continuous(breaks = seq(1980, 2010, 5), labels = c('1975–\n1984', '1980–\n1989', '1985–\n1994', '1990–\n1999', '1995–\n2004', '2000–\n2009', '2005–\n2014')) +
  xlab('Time period') + ylab(expression(paste('Aggregate neonicotinoid effect (% ', year^-1, ' ', neonic^-1, ')'))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), 'lines'),
        axis.title.x = element_text(margin = margin(0.8, 0, 0, 0, unit = 'lines')))

# check plot
dev.off()
quartz(height = 5.5, width = 6.5)
print(fig_s5)

# write to file
ggsave('figures/fig-s5.png', fig_s5, height = 5.5, width = 6.5, units = 'in', dpi = 300)





######### Figure S6 and S7 (neonicotinoid effect over time, by species or bcr)

# load relevant data
load('analysis/post-summary-agg-spatial-neonic.RData')  # summary_agg_neonic

# organize data
df_spp_neonic_plot <- summary_agg_neonic %>% 
  dplyr::select(year, df_spp) %>% unnest() %>% 
  left_join(dplyr::select(species_df, aou, species_short), by = 'aou')

df_bcr_neonic_plot <- summary_agg_neonic %>% 
  dplyr::select(year, df_bcr) %>% unnest() %>% 
  left_join(bcr_dat, by = 'bcr')

# reorder levels of species or bcr by median value in 2005-2014 time period
lev_bcr <- rev(df_bcr_neonic_plot$bcr_short[df_bcr_neonic_plot$year == 2010][order(df_bcr_neonic_plot$theta_med[df_bcr_neonic_plot$year == 2010])])
df_bcr_neonic_plot$bcr_short <- factor(df_bcr_neonic_plot$bcr_short, levels = lev_bcr)

lev_spp <- rev(df_spp_neonic_plot$aou[df_spp_neonic_plot$year == 2010][order(df_spp_neonic_plot$theta_med[df_spp_neonic_plot$year == 2010])])
df_spp_neonic_plot$aou <- factor(df_spp_neonic_plot$aou, levels = lev_spp)

# create dummy plots so that 0 is vertically centered in each panel
plot_bcr_dummy <- df_bcr_neonic_plot %>%
  group_by(bcr_short) %>% 
  summarize(plot_upp = max(abs(c(theta_low90, theta_upp90)))) %>% 
  mutate(plot_low = -plot_upp, x = 2000) %>% 
  gather(type, val, plot_low:plot_upp)

plot_spp_dummy <- df_spp_neonic_plot %>%
  group_by(aou) %>%
  summarize(plot_upp = max(abs(c(theta_low90, theta_upp90)))) %>%
  mutate(plot_low = -plot_upp, x = 2000) %>%
  gather(type, val, plot_low:plot_upp)

# plot neonicotinoid effect over time, by species
fig_s6 <- ggplot(df_spp_neonic_plot) +
  geom_point(aes(year, theta_med), size = 1.8) +
  geom_linerange(aes(year, ymin = theta_low90, ymax = theta_upp90), size = 1.2) +
  geom_linerange(aes(year, ymin = theta_low99, ymax = theta_upp99), size = 0.4) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_blank(data = plot_spp_dummy, aes(x, val)) +
  scale_x_continuous(breaks = seq(1980, 2010, 5), labels = c("1975-\n1984", '', "1985-\n1994", '', "1995-\n2004", '', "2005-\n2014")) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  coord_cartesian(xlim = c(1978, 2012), ylim = c(-3.5, 3.5)) +
  facet_wrap(~ species_short, nrow = 4) +
  xlab('Time period') + ylab(expression(paste('Neonicotinoid effect (% ', year^-1, ' ', neonic^-1, ')'))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 10.5),
        plot.margin = unit(c(0.3, 0.5, 0.7, 0.5), 'lines'),
        axis.title.x = element_text(margin = margin(0.8, 0, 0, 0, unit = 'lines')))

# plot neonicotinoid effect over time, by bcr
fig_s7 <- ggplot(df_bcr_neonic_plot) +
  geom_point(aes(year, theta_med), size = 1.8) +
  geom_linerange(aes(year, ymin = theta_low90, ymax = theta_upp90), size = 1.2) +
  geom_linerange(aes(year, ymin = theta_low99, ymax = theta_upp99), size = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, size = 0.3, alpha = 0.8) +
  geom_blank(data = plot_bcr_dummy, aes(x, val)) +
  scale_x_continuous(breaks = seq(1980, 2010, 5), labels = c('1975-\n1984', '', '1985-\n1994', '', '1995-\n2004', '', '2005-\n2014')) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  coord_cartesian(xlim = c(1978, 2012), ylim = c(-3.5, 3.5)) +
  facet_wrap(~ bcr_short, nrow = 5) +
  xlab('Time period') + ylab(expression(paste('Neonicotinoid effect (% ', year^-1, ' ', neonic^-1, ')'))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 10.5),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 10.5),
        plot.margin = unit(c(0.3, 0.5, 0.7, 0.5), 'lines'),
        axis.title.x = element_text(margin = margin(0.8, 0, 0, 0, unit = 'lines')))

# check plots
dev.off()
quartz(height = 8, width = 12)
print(fig_s6)

dev.off()
quartz(height = 8, width = 12)
print(fig_s7)

# write to file
ggsave('figures/fig-s6.png', fig_s6, height = 8, width = 12, units = 'in', dpi = 300)
ggsave('figures/fig-s7.png', fig_s7, height = 8, width = 12, units = 'in', dpi = 300)





######### Figure S9 (neonic use over time, in USA and Netherlands)

# list of neonic compounds
neonics <- c('ACETAMIPRID', 'CLOTHIANIDIN', 'DINOTEFURAN', 'IMIDACLOPRID', 'THIACLOPRID', 'THIAMETHOXAM')

# neonic use by county and year, scaled by land area
neonic_county_year <- read_csv('data/neonic-county-year.csv') %>% 
  mutate(type = 'American counties')

# scaled nationwide neonic use in netherlands and usa
neonic_netherlands <- read_csv('data/neonic-use-netherlands.csv') %>% 
  group_by(year) %>%
  summarize(neonic_tot = sum(kg, na.rm = T)) %>% 
  mutate(neonic_scaled = neonic_tot / 33690, country = 'Netherlands')

neonic_usa <- read_csv('data/epest-merge.csv') %>%
  filter(compound %in% neonics) %>%
  group_by(year) %>% 
  summarize(neonic_tot = sum(kg, na.rm = T)) %>% 
  mutate(neonic_scaled = neonic_tot / 7663941, country = 'United States')

neonic_full <- rbind.data.frame(neonic_netherlands, neonic_usa)

# as of 2012, how many american counties have scaled neonic use greater than netherlands-wide rate
neonic_county_year %>% 
  filter(year == 2012, neonics >= filter(neonic_netherlands, year == 2012)$neonic_scaled)

# create figure
fig_s9 <- ggplot(neonic_full, aes(year, neonic_scaled, col = country)) +
  geom_line(inherit.aes = F, data = neonic_county_year, aes(year, neonics, group = fips, linetype = type), alpha = 0.1) +
  geom_point(size = 2.5) +
  geom_line(size = 1.5) +
  coord_cartesian(ylim = c(0, 3)) +
  xlab('Year') + ylab(expression(paste('Neonicotinoid use (kg ', km^-2, ')'))) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 13),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, 'cm'),
        legend.margin = margin(0, 0, 0, 0, unit = 'cm'),
        axis.ticks = element_line(size = 0.3),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')),
        legend.position = c(0.04, 0.75),
        legend.justification = c(0, 0))

# check figure
dev.off()
quartz(height = 4.5, width = 5.5)
print(fig_s9)

# write to file
ggsave('figures/fig-s9.png', fig_s9, height = 4.5, width = 5.5, units = 'in')





################### Figure A5
# correlation between posterior median neonic effects estimated in main analysis
# versus effects estimated using spatial mismatch model

# load relevant data
load('analysis/post-summary-agg-spatial-neonic-mismatch.RData')  # summary_agg_neonic_mismatch
load('analysis/post-summary-agg-spatial-neonic.RData')           # summary_agg_neonic

# organize data
df_mismatch_spp <- summary_agg_neonic_mismatch$df_spp[[1]] %>% mutate(type = 'Mismatch')
df_mismatch_bcr <- summary_agg_neonic_mismatch$df_bcr[[1]] %>% mutate(type = 'Mismatch')

df_neonic_spp <- summary_agg_neonic$df_spp[[7]] %>% mutate(type = 'Original')
df_neonic_bcr <- summary_agg_neonic$df_bcr[[7]] %>% mutate(type = 'Original')

df_spp_full <- rbind.data.frame(df_mismatch_spp, df_neonic_spp) %>% 
  dplyr::select(aou, type, theta_med) %>% 
  spread(type, theta_med)

df_bcr_full <- rbind.data.frame(df_mismatch_bcr, df_neonic_bcr) %>% 
  dplyr::select(bcr, type, theta_med) %>% 
  spread(type, theta_med)

# ggplot theme
tt <- theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size = 14),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 15),
        axis.ticks = element_line(size = 0.3),
        axis.title.x = element_text(margin = margin(.3, 0, 0, 0, unit = 'cm')),
        axis.title.y = element_text(margin = margin(0, .3, 0, 0, unit = 'cm')))

# left panel (by species)
p1 <- ggplot(df_spp_full, aes(Original, Mismatch)) +
  geom_point(size = 2) +
  scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2)) +
  xlab('Median neonic. effect (initial analysis)') +
  ylab('Median neonic. effect (amalgamated counties)') +
  ggtitle('By species') +
  tt

# right panel (by BCR)
p2 <- ggplot(df_bcr_full, aes(Original, Mismatch)) +
  geom_point(size = 2) +
  xlab('Median neonic. effect (initial analysis)') +
  ylab(NULL) +
  ggtitle('By BCR') +
  tt

# full plot
p_full <- arrangeGrob(p1, p2, nrow = 1, widths = c(1.08, 1))

# check figure
dev.off()
quartz(height = 6, width = 10)
grid.arrange(p_full)

# write to file
ggsave('figures/appendix-regular-vs-amalgam-counties.png', p_full, height = 6, width = 10, units = 'in')


