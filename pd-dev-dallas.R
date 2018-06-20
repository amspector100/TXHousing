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

# Graph 1: The Choropleth -----------------------------------------------------

# Part 1.1: Median Income Data -----------------------------------------------
# Get shapes/shapefile for census tracts
tx_tract_shapes_path <- "cb_2017_48_TX_tract_500k/Cb_2017_48_tract_500k.shp"
tract_shapes <- st_read(tract_shapes_path) 
# Change coordinate system to lat-long
tract_shapes <- st_transform(tract_shapes, '+init=EPSG:4326')

# Get data from census
tract_data_path <- 'ACS_16_5YR_S1901/ACS_16_5YR_S1901.csv'
# HC01_EST_VC13 is for median income
tract_data <- read_csv(tract_data_path) %>%
  mutate(median_income = HC01_EST_VC13/1000) %>% # Keep median income in thousands of $
  select(c('GEO.id', 'median_income'))

# Join data
all_data <- inner_join(tract_data, tract_shapes, by = c('GEO.id' = 'AFFGEOID'))
all_data <- st_sf(all_data)

# Part 1.2: Zoning Data -------------------------------------------------------------
dallas_shp_path <- paste('Zoning Shapefiles/', dallas_inputs$path, sep='')
zoning_data <- st_read(dallas_shp_path)
zoning_data <- st_transform(zoning_data, '+init=EPSG:4326') %>%
  tidyr::separate_(dallas_inputs$feature, into = c('zone_code'), sep='-') %>%
  filter(zone_code == 'PD')

# Part 1.3: Graphing ----------------------------------------------------------------

library(ggmap)
mapImage <- get_map(location = c(lon = -96.8, lat = 32.8),
                    color = "color",
                    source = "google",
                    #maptype = "terrain",
                    zoom = 11)

graph <- ggmap(mapImage) +
  geom_sf(data = all_data, aes(fill = median_income), color = NA, inherit.aes = FALSE,
          alpha = 0.6) + 
  scale_fill_gradient(low = "blue", high = "red") +
  geom_sf(data = zoning_data, aes(), fill = NA,
          color = 'black', inherit.aes = FALSE) +
  labs(x = 'Longitude', y = 'Latitude', title = 'Planned Developement Zones and Median Income 
       in Dallas Texas', caption = expression(paste('Planned Development zones outlined in black. 
Data for median income from 5 year ACS estimates in 2016, Dallas City Open Data.')))

#setwd('C:/Users/aspector/Documents/R Projects/TX Housing')
#png('Figures/pd-dev-dallas.png', w=800, h=500)
#print(graph)
#dev.off()

#svg('Figures/pd-dev-dallas.svg', w=8, h=5)
#print(graph)
#dev.off()

# Graph 2: The Histogram -----------------------------------------------------

# Part 2.1: Get income data from census

# Start with income data

old_cols <- c('HC01_VC75', 'HC01_VC76', 'HC01_VC77', 'HC01_VC78',
              'HC01_VC79', 'HC01_VC80', 'HC01_VC81', 'HC01_VC82',
              'HC01_VC83', 'HC01_VC84')
select_cols <- append(c('GEO.id'), old_cols)
new_cols <- c('Under 10K', '10-15K', '15-25K', '25-35K', 
              '35-50K', '50-75K', '75-100K', '100-150K', 
              '150-200K', '200K+')

data_path <- "ACS_16_5YR_DP03/ACS_16_5YR_DP03.csv"
tract_data <- read_csv(data_path) %>%
  select(select_cols) %>%
  rename_at(vars(old_cols), ~new_cols)

# Get shapes/shapefile for census tracts
tx_tract_shapes_path <- "cb_2017_48_TX_tract_500k/Cb_2017_48_tract_500k.shp"
tract_shapes <- st_read(tract_shapes_path) 
# Change coordinate system to lat-long
tract_shapes <- st_transform(tract_shapes, '+init=EPSG:4326')

census_data <- inner_join(tract_data, tract_shapes, by = c('GEO.id' = 'AFFGEOID')) #I have checked this does not yield duplicates
census_data <- st_sf(census_data) 

# Now we'll divide by the land area (in square meters) to get the density per
# square meter by income group in each tract. 

all_data <- census_data %>% gather('variable', 'value', new_cols) %>% 
  mutate(value = as.numeric(value)/ALAND) %>%
  spread(variable, value)

# Part 2.2: Get zoning data and merge with census tract data ---------------
dallas_shp_path <- paste('Zoning Shapefiles/', dallas_inputs$path, sep='')
zoning_data <- st_read(dallas_shp_path)
zoning_data <- st_transform(zoning_data, '+init=EPSG:4326') %>%
  tidyr::separate_(dallas_inputs$feature, into = c('zone_code'), sep='-') %>%
  filter(zone_code == 'PD')

# Now merge using the intersect feature. This is a bit inaccurate because we're using
# lat-long coords, but because it's just one city which isn't too close to one of
# the poles, it should be okay.

st_crs(zoning_data) <- '+init=EPSG:4326'
st_crs(all_data) <- '+init=EPSG:4326'
merged_data <- sf::st_intersection(all_data, zoning_data)

# Ongoing question - what about lake and river area? This is tough to account for
# but maybe if I can find a (small enough) shapefile of lakes/rivers in Dallas,
# I could use the st_difference function perhaps...

# Calculate area of intersections, in sq meters
merged_data$AREA <- as.numeric(st_area(merged_data))

# Quick sanity check: Since each part of each zoning district ----------------
# should be in a census tract, make sure the sum of the areas of the intersections
# is approximately equal to the sum of the area of the initial zoning data.
# The merged_data area is 96.3% of the zoning_data area, which isn't bad
# because the intersection function isn't great at handling lat/long coords.

ratio <- sum(merged_data$AREA)/sum(st_area(zoning_data))

# End of quick sanity check ------------------------------------------------

# Multiply densities by area of intersections and sum. First rename columns
# because the intersect function changes column names. Then summarize.

intersect_cols <- c('Under.10K', 'X10.15K', 'X15.25K', 'X25.35K', 
                    'X35.50K', 'X50.75K', 'X75.100K', 'X100.150K', 
                    'X150.200K', 'X200K.')

final_PD_data <- as.data.frame(merged_data) %>% 
  rename_at(vars(intersect_cols), ~ new_cols) %>%
  gather('variable', 'value', new_cols) %>%
  mutate(value = as.numeric(value)*AREA) %>%
  spread(variable, value) %>%
  select(new_cols) %>%
  summarise_all(sum) 

final_PD_data <- 100*final_PD_data/sum(final_PD_data)

# Do the same for city-wide income data
final_city_data <- as.data.frame(census_data) %>% 
  select(new_cols) %>%
  mutate_all(as.numeric) %>%
  summarise_all(sum)

final_city_data <- 100*final_city_data/sum(final_city_data)

# Part 3: Graphing ----------------------------------------------------

# Bind data together
final_data <- rbind(final_city_data, final_PD_data)
final_data$Location <- c("Entire City", "Planned Development Zones") 
final_data <- final_data %>% gather(variable, value, new_cols)
  
# Make sure order is right
final_data$variable <- factor(final_data$variable, 
                                 levels =  new_cols)

# Graph
gg <- ggplot(final_data) + 
  geom_bar(aes(x = variable, y = value, group = Location,
               fill = Location), 
           stat = 'identity',
           position = 'dodge') + 
  labs(x = 'Income Bracket', y = 'Percent of Population',
       title = 'Income Distribution of Dallas and its Planned Development Zones',
       caption = 'Data from 2016 ACS Estimates and the City of Dallas GIS Services')

#setwd('C:/Users/aspector/Documents/R Projects/TX Housing')
#png('Figures/pd-dev-dallas-hist.png', w=800, h=500)
#print(gg)
#dev.off()

#svg('Figures/pd-dev-dallas-hist-nonscaled.svg', w=8, h=5)
#print(gg)
#dev.off()
