#setwd("C:/Users/aspector/Documents/R Projects/TX Housing/data/")
# Make sure you run the zoning-maps-inputs.R file to get some global settings for the 
# dallas zoning data, i.e.
#options(run.main = FALSE)
#source('zoning-maps-inputs.R')
#options(run.main = TRUE)

library(sf)
library(tidyverse)
library(broom)
library(maptools)

# Part 1: Median Income Data -----------------------------------------------
# Get shapes/shapefile for census tracts
tract_shapes_path <- "cb_2017_48_tract_500k_TX_Census_Tract_Shapefile/Cb_2017_48_tract_500k.shp"
tract_shapes <- st_read(tract_shapes_path) 
# Change coordinate system to lat-long
tract_shapes <- st_transform(tract_shapes, '+init=EPSG:4326')

# Get data from census
tract_data_path <- 'ACS_16_5YR_S1901_Dallas_Income_Statistics/ACS_16_5YR_S1901.csv'
# HC01_EST_VC13 is for median income
tract_data <- read_csv(tract_data_path) %>%
  mutate(median_income = HC01_EST_VC13/1000) %>% # Keep median income in thousands of $
  select(c('GEO.id', 'median_income'))

# Join data
all_data <- inner_join(tract_data, tract_shapes, by = c('GEO.id' = 'AFFGEOID'))

# Part 2: Zoning Data -------------------------------------------------------------
dallas_shp_path <- paste('Zoning Shapefiles/', dallas_inputs$path, sep='')
zoning_data <- st_read(dallas_shp_path)
zoning_data <- st_transform(zoning_data, '+init=EPSG:4326') %>%
  tidyr::separate_(dallas_inputs$feature, into = c('zone_code'), sep='-') %>%
  filter(zone_code == 'PD')

# Part 3: Graphing ----------------------------------------------------------------

library(ggmap)
mapImage <- get_map(location = c(lon = -96.8, lat = 32.8),
                    color = "color",
                    source = "google",
                    #maptype = "terrain",
                    zoom = 10.8)

graph <- ggmap(mapImage) +
  geom_sf(data = all_data, aes(fill = median_income), color = NA, inherit.aes = FALSE,
          alpha = 0.6) + 
  scale_fill_gradient(low = "blue", high = "red") +
  geom_sf(data = zoning_data, aes(), fill = NA,
          color = 'black', inherit.aes = FALSE) +
  labs(x = 'Longitude', y = 'Latitude', title = 'Planned Developement Zones and Median Income 
       in Dallas Texas', caption = expression(paste('Planned Development zones overlaid in light yellow. 
Data for median income from 5 year ACS estimates in 2016, Dallas City Open Data.')))

#setwd('C:/Users/aspector/Documents/R Projects/TX Housing')
png('Figures/pd-dev-dallas.png', w=800, h=500)
print(graph)
dev.off()

svg('Figures/pd-dev-dallas.svg', w=8, h=5)
print(graph)
dev.off()

