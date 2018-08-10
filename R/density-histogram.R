library(tidyverse)
library(broom)
library(maptools)
library(sf)

# Source is 2012-2016 ACS block level data.https://www.census.gov/geo/maps-data/data/tiger-data.html

# Both versions need these paths
tx_places_shapes_path <- "data/Census/cb_2017_48_place_500k/cb_2017_48_place_500k.shp"
ca_places_shapes_path <- "data/Census/tl_2017_06_CA_place/tl_2017_06_place.shp"
ma_places_shapes_path <- "data/Census/tl_2017_25_MA_place/tl_2017_25_place.shp"
ny_places_shapes_path <- "data/Census/tl_2017_36_NY_place/tl_2017_36_place.shp"
il_places_shapes_path <- "data/Census/tl_2017_17_IL_place/tl_2017_17_place.shp"

# Reads and transforms
process_shapes <- function(path){
  shapes <- st_read(path)
  shapes <- st_transform(shapes, '+init=EPSG:4326')
  return(shapes)
}

# Note that these files are huge, so only read some layers at a time.
# These are geodatabases, not shapefiles.

# Paths
tx_block_path <- "data/Census/ACS_2016_5YR_BG_48_TEXAS.gdb"
ny_block_path <- "data/Census/ACS_2016_5YR_BG_36_NEW_YORK.gdb"
ca_block_path <- "data/Census/ACS_2016_5YR_BG_06_CALIFORNIA.gdb"
il_block_path <- "data/Census/ACS_2016_5YR_BG_17_ILLINOIS.gdb"
ma_block_path <- "data/Census/ACS_2016_5YR_BG_25_MASSACHUSETTS.gdb"

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
  place_shapes <- process_shapes(place_path) %>%
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

svg('Figures/comparison-densities.svg', w = 8, h = 8)
print(graph)
dev.off()


