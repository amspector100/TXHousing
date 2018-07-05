#setwd("C:/Users/aspector/Documents/R Projects/TX Housing/data/")
# Make sure you run the zoning-maps-inputs.R file to get some global settings for the 
# dallas zoning data, i.e.
#options(run.main = FALSE)
#source('zoning-maps-inputs.R')
#options(run.main = TRUE)

library(tidyverse)
library(broom)
library(maptools)
library(sf)

# Both versions need these paths
tx_places_shapes_path <- "tl_2017_48_TX_place/tl_2017_48_place.shp"
ca_places_shapes_path <- "tl_2017_06_CA_place/tl_2017_06_place.shp"
ma_places_shapes_path <- "tl_2017_25_MA_place/tl_2017_25_place.shp"
ny_places_shapes_path <- "tl_2017_36_NY_place/tl_2017_36_place.shp"
il_places_shapes_path <- "tl_2017_17_IL_place/tl_2017_17_place.shp"

# Get shapes/shapefile for census tracts by state
process_tract_shapes <- function(path){
  tract_shapes <- st_read(path)
  tract_shapes <- st_transform(tract_shapes, '+init=EPSG:4326')
  return(tract_shapes)
}

# Version 1: Census tracts

# Final pre-merge columns
premerge_columns <- c('NAME', 'ALAND', 'AWATER', 'AFFGEOID')


# Intersects them with municipality boundaries by state
intersect_tract_shapes <- function(tract_path, place_path, names){
  tract_shapes <- process_tract_shapes(tract_path) %>%
    select(ALAND, AWATER, AFFGEOID)
  places_shapes <- process_tract_shapes(place_path) %>%
    filter(NAME %in% names)
  result <- st_intersection(tract_shapes, places_shapes) %>%
    select(premerge_columns)
  return(st_as_sf(result))
}

# Get metropolitan areas

# Get census tract shapes. This involves downloading tract shapefiles by county
# and then taking the intersection with city shapefiles. 
tx_tract_shapes_path <- "cb_2017_48_TX_tract_500k/cb_2017_48_tract_500k.shp"
tx_tracts <- intersect_tract_shapes(tx_tract_shapes_path,
                                    tx_places_shapes_path,
                                    c('Austin', 'Houston', 'Dallas'))
tx_county_tracts <- process_tract_shapes(tx_tract_shapes_path) %>%
  filter(COUNTYFP %in% c('201', '453', '113')) %>%
  mutate(NAME = case_when(
    COUNTYFP == '201' ~ 'Harris County',
    COUNTYFP == '453' ~ 'Travis County',
    COUNTYFP == '113' ~ 'Dallas County',
    TRUE ~ 'NA'
  )) %>%
  select(ALAND, AWATER, AFFGEOID, NAME)


# This one's a bit different because we're looking at all of LA county at one point
ca_tract_shapes_path <- "cb_2017_06_CA_tract_500k/cb_2017_06_tract_500k.shp"
ca_city_tracts <- intersect_tract_shapes(ca_tract_shapes_path,
                                         ca_places_shapes_path,
                                         c('San Francisco', 'Los Angeles'))
la_county_tracts <- process_tract_shapes(ca_tract_shapes_path) %>%
  filter(COUNTYFP == '037') %>%
  mutate(NAME = 'Los Angeles County') %>%
  select(ALAND, AWATER, AFFGEOID, NAME)

ma_tract_shapes_path <- "cb_2017_25_MA_tract_500k/cb_2017_25_tract_500k.shp"
ma_tracts <- intersect_tract_shapes(ma_tract_shapes_path,
                                    ma_places_shapes_path,
                                    c('Boston'))


ny_tract_shapes_path <- "cb_2017_36_NY_tract_500k/cb_2017_36_tract_500k.shp"
ny_tract_shapes <- process_tract_shapes(ny_tract_shapes_path)
ny_tracts <- intersect_tract_shapes(ny_tract_shapes_path,
                                    ny_places_shapes_path,
                                    c('New York'))

il_tract_shapes_path <- "cb_2017_17_IL_tract_500k/cb_2017_17_tract_500k.shp"
il_tract_shapes <- process_tract_shapes(il_tract_shapes_path)
il_tracts <- intersect_tract_shapes(il_tract_shapes_path,
                                    il_places_shapes_path,
                                    c('Chicago'))

all_tracts <- rbind(tx_tracts, tx_county_tracts, ca_city_tracts, la_county_tracts, 
                    ny_tracts, ma_tracts, il_tracts)

# Get population numbers for each census tract
pop_data_path <- 'Density_Comparison_ACS_16_5YR_DP05/ACS_16_5YR_DP05.csv'
old_cols <- c('GEO.id', 'HC01_VC03')
new_cols <- c('GEO.id', 'population')
raw_pop_data <- read_csv(pop_data_path) %>% 
  select(old_cols) %>%
  rename_at(vars(old_cols), ~ new_cols)


# Merge all the data - also drop unused levels in the name column
all_data <- inner_join(raw_pop_data, all_tracts, by = c('GEO.id' = 'AFFGEOID'))
all_data <- st_sf(all_data)
all_data$NAME <- droplevels(all_data$NAME)

# Calculate areas of each location
all_data <- all_data %>%
  # Get area, in square meters
  mutate(area = as.numeric(st_area(all_data))) %>%
  # If possible, discount the area from water (possible for roughly 2/3 of tracts).
  # This is possible because a lot of the tracts are totally inside the municipality
  # borders.
  mutate(area = ifelse((ALAND + AWATER - area)/(ALAND + AWATER) < 0.01, 
                       ALAND, area*ALAND/(ALAND + AWATER))) %>%
  
  # Get rid of the tracts which barely intersect the city boundaries
  filter(area/(ALAND + AWATER) > 0.02) %>%
  mutate(area = 0.00024711*area) %>%   # Convert area to acres
  mutate(density = as.numeric(population)/area) %>% # 
  mutate(density = ifelse(density > 100, 100, density)) %>% # Group outliers together

  # Filter out areas that are non-residential by excluding less than 
  #0.25 people/1000 square meters, which is roughly 1 person per acre
  filter(density > 1) 
  # Then put locations in the right order

# Quick reality check
reality_check <- as.data.frame(all_data) %>% 
  select(NAME, population, area) %>% 
  group_by(NAME) %>%
  summarize(pop_sum = sum(as.numeric(population)), area_sum = sum(as.numeric(area)))
  

gg <- ggplot(all_data, aes(density)) + 
  facet_wrap(. ~ NAME, ncol = 3) +
  geom_histogram(aes(y = (..count..)/ sapply(PANEL, FUN=function(x) sum(count[PANEL == x]))),
                 binwidth = 5, fill = 'cornflowerblue') + 
  labs(x='Population per Acre', y='Percent of Census Tracts',
       title = 'Population Density Histograms by Location') + 
  scale_y_continuous(labels = scales::percent)


setwd("C:/Users/aspector/Documents/R Projects/TX Housing")
svg('Figures/comparison-densities.svg', w = 8, h = 5)
print(gg)
dev.off()

# Version 2 -- Block level data. See https://www.census.gov/geo/maps-data/data/tiger-data.html

# Note that these files are huge, so only read some layers at a time.
# These are geodatabases, not shapefiles.

#setwd("C:/Users/aspector/Documents/R Projects/TX Housing/data/").

# Paths
tx_block_path <- "ACS_2016_5YR_BG_48_TEXAS.gdb"
ny_block_path <- "ACS_2016_5YR_BG_36_NEW_YORK.gdb"
ca_block_path <- "ACS_2016_5YR_BG_06_CALIFORNIA.gdb"
il_block_path <- "ACS_2016_5YR_BG_17_ILLINOIS.gdb"
ma_block_path <- "ACS_2016_5YR_BG_25_MASSACHUSETTS.gdb"

# Functions for reading in data.
premerge_columns <- c('NAME', 'ALAND', 'AWATER', 'population')

process_block_shapes <- function(block_path, 
                                 data_layer, 
                                 geo_layer, 
                                 feature = 'B00001e1'
                                 ) {
  data <- st_read(block_path, data_layer) %>% 
    select(feature, 'GEOID') %>%
    rename_at(vars(feature), ~'population') 
  geo <- st_read(block_path, geo_layer)
  geodata <- st_as_sf(inner_join(data, geo, by = c('GEOID' = 'GEOID_Data')))
  geodata <- st_transform(geodata, '+init=EPSG:4326')
  return(geodata)
}

intersect_block_shapes <- function(block_path, place_path, data_layer = NA,
                                   geo_layer = NA, feature = 'B00001e1', names = 'New York'){
  block_shapes <- process_block_shapes(block_path, 
                                       data_layer, 
                                       geo_layer, 
                                       feature = feature)
  place_shapes <- process_tract_shapes(place_path) %>%
    filter(NAME %in% names)
  result <- st_intersection(block_shapes, place_shapes) %>%
    select(ALAND, AWATER, population, NAME)
  result$NAME <- droplevels(result$NAME)
  return(result)
}

ma_blocks <- intersect_block_shapes(ma_block_path, ma_places_shapes_path,
                                    data_layer = 'X01_AGE_AND_SEX', 
                                    geo_layer = 'ACS_2016_5YR_BG_25_MASSACHUSETTS',
                                    feature ='B01001e1',
                                    names = 'Boston')

ny_blocks <- intersect_block_shapes(ny_block_path, ny_places_shapes_path,
                                    data_layer = 'X01_AGE_AND_SEX', 
                                    geo_layer = 'ACS_2016_5YR_BG_36_NEW_YORK',
                                    feature ='B01001e1',
                                    names = c('New York'))

il_blocks <- intersect_block_shapes(il_block_path, il_places_shapes_path,
                                    data_layer = 'X01_AGE_AND_SEX', 
                                    geo_layer = 'ACS_2016_5YR_BG_17_ILLINOIS',
                                    feature ='B01001e1',
                                    names = c('Chicago'))

ca_blocks <- intersect_block_shapes(ca_block_path, ca_places_shapes_path,
                                    data_layer = 'X01_AGE_AND_SEX', 
                                    geo_layer = 'ACS_2016_5YR_BG_06_CALIFORNIA',
                                    feature ='B01001e1',
                                    names = c('San Francisco', 'Los Angeles'))

ca_county_blocks <- process_block_shapes(ca_block_path,
                                         data_layer = 'X01_AGE_AND_SEX', 
                                         geo_layer = 'ACS_2016_5YR_BG_06_CALIFORNIA',
                                         feature ='B01001e1') %>% 
  filter(COUNTYFP == '037') %>%
  mutate(NAME = 'Los Angeles County') %>%
  select(premerge_columns)
  
tx_blocks <- intersect_block_shapes(tx_block_path, tx_places_shapes_path,
                                    data_layer = 'X01_AGE_AND_SEX', 
                                    geo_layer = 'ACS_2016_5YR_BG_48_TEXAS',
                                    feature ='B01001e1',
                                    names = c('Dallas', 'Austin', 'Houston'))

tx_county_blocks <- process_block_shapes(tx_block_path,
                                         data_layer = 'X01_AGE_AND_SEX', 
                                         geo_layer = 'ACS_2016_5YR_BG_48_TEXAS',
                                         feature ='B01001e1') %>%
  filter(COUNTYFP %in% c('201', '453', '113')) %>%
  mutate(NAME = case_when(
    COUNTYFP == '201' ~ 'Harris County',
    COUNTYFP == '453' ~ 'Travis County',
    COUNTYFP == '113' ~ 'Dallas County',
    TRUE ~ 'NA'
  )) %>%
  select(premerge_columns) %>%
  mutate(NAME = as.factor(NAME))

binwidth = 10
maximum = 100

all_blocks <- as.data.frame(rbind(tx_blocks, tx_county_blocks, ca_blocks, ca_county_blocks,
                    il_blocks, ma_blocks, ny_blocks)) %>%
  # Get density in persons per acre
  mutate(density = as.numeric(population)/(0.00024711*ALAND)) %>%
  # Group outliers together
  mutate(density = ifelse(density > maximum, maximum, density)) %>%
  # Create bins
  mutate(bin = droplevels(as.factor(binwidth*(density %/% binwidth)))) #%>%
  #mutate(bin = paste(as.character(bin), as.character(bin + binwidth - 0.1), sep = '-'))

# Get total population counts for each area
area_sums <- as.data.frame(all_blocks) %>%
  group_by(NAME) %>%
  summarise_at(c('population'), sum) %>%
  ungroup()

# Get percent of population in each bin
bin_sums <- as.data.frame(all_blocks) %>%
  group_by(NAME, bin) %>%
  summarise(bin_sum = sum(population)) %>%
  filter(!is.na(bin)) %>%
  ungroup()

# Get percents
final_data <- inner_join(area_sums, bin_sums, by = c("NAME" = "NAME")) %>%
  mutate(percent = 100*bin_sum/population) %>%
  mutate(NAME = factor(NAME, 
                       levels = c('Austin', 'Dallas', 'Houston', 
                       'Travis County', 'Dallas County', 'Harris County',
                       'Los Angeles', 'Los Angeles County', 'San Francisco',
                       'Chicago', 'Boston', 'New York')))

# Graph
bar_width = 0.9
graph <- ggplot(final_data, aes(x = bin, y = percent)) +
  geom_col(position=position_dodge(bar_width), fill = 'cornflowerblue') + 
  facet_wrap(. ~ NAME, ncol = 3) + 
  labs(x = 'Population per Acre', y = 'Percent of Population in Municipality',
       caption = 'Data from 2016 5-Year ACS Estimates (Block-level)',
       title = 'Population Density by Municipality')
#print(graph)

#setwd("C:/Users/aspector/Documents/R Projects/TX Housing")
#svg('Figures/comparison-densities-block.svg', w = 8, h = 5)
#print(graph)
#dev.off()


