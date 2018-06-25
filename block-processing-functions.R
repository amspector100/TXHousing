
library(tidyverse)
library(maptools)
library(ggmap)
library(sf)

# See https://www2.census.gov/geo/tiger/TIGER_DP/2016ACS/Metadata/BG_METADATA_2016.txt
# for metadata

# Process a block shape - this is slightly modified in some other scripts
process_block_shapes <- function(block_path, 
                                 data_layer, 
                                 geo_layer, 
                                 old_cols,
                                 new_cols) {
  data <- st_read(block_path, data_layer, quietly = TRUE) %>% 
    select(feature, 'GEOID') %>%
    rename_at(vars(old_cols), new_cols)
  geo <- st_read(block_path, geo_layer)
  geodata <- st_as_sf(inner_join(data, geo, by = c('GEOID' = 'GEOID_Data')))
  geodata <- st_transform(geodata, '+init=EPSG:4326') 
  return(geodata)
}