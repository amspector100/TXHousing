#setwd("C:/Users/aspector/Documents/R Projects/TX Housing")


library(tidyverse)
library(broom)
library(maptools)
library(sf)

# Get shapes/shapefile for census tracts by state
process_tract_shapes <- function(path){
  tract_shapes <- st_read(path)
  tract_shapes <- st_transform(tract_shapes, '+init=EPSG:4326')
  return(tract_shapes)
}

tx_places_shapes_path <- "data/tl_2017_48_TX_place/tl_2017_48_place.shp"
austin_boundary <- process_tract_shapes(tx_places_shapes_path) %>%
  filter(NAME == 'Austin') %>%
  select(NAME)

# Already in lat/long coords - won't retransform because it's large (245K rows, 65 cols)
austin_path <- 'data/Zoning Shapefiles/Austin Land Database 2016/geo_export_ade599f6-7cc2-4e34-b51e-3aaccd6a2dc1.shp'
land_use_data <- st_read(austin_path)

st_crs(land_use_data) <- '+init=EPSG:4326'
st_crs(austin_boundary) <- '+init=EPSG:4326'
# I assume this takes a long time
subsetted_data <- st_as_sf(st_intersection(austin_boundary, land_use_data))


