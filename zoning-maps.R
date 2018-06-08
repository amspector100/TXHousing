#setwd("C:/Users/aspector/Documents/R Projects/TX Housing/data/Zoning Shapefiles")

# This will prevent the zoning-maps-inputs script from doing extra work
options(run.main = FALSE)
source('zoning-maps-inputs.R')
options(run.main = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(rgdal)
library(maptools)

coarse_format_data <- function(raw_data, sf, mf, or, feature, 
                               separator = '-', path = NA, read_shapefile = FALSE){
  
  # -----------------------------------------------------------
  # Formats shapefile data. Assumes that the SpatialPolygons object will have an 
  # @data slot. The 'feature' should be the zoning data column name. 
  # -----------------------------------------------------------
  
  # If specified, read in the shapefile
  if (read_shapefile){
    raw_data <- readShapePoly(path) # This operation is expensive
  }
  
  # Else, just manipulate the data
  zoning_data <- raw_data@data %>%
    tidyr::separate_(feature, into = c('zone_code'), sep=separator) %>% 
    select('zone_code') %>% 
    mutate(zone_code = case_when(
      zone_code %in% sf ~ 'SF',
      zone_code %in% mf ~ 'MF',
      zone_code %in% or ~ 'OR',
      TRUE ~ 'NR'
    ))
    
  return(zoning_data)
}

#raw_dallas_data <- readShapePoly(dallas_inputs$path)
#raw_dallas_data@proj4string <- dallas_inputs$proj4string
#transformed_dallas_data <- spTransform(raw_dallas_data, lat_long)


test_inputs <- north_texas_inputs

print('Test one, manipulating data')

raw_data <- readShapePoly(test_inputs$path, delete_null_obj = TRUE)
raw_data@proj4string <- test_inputs$proj4string
transformed_data <- spTransform(raw_data, lat_long)

print('Test two, graphing data')

library(ggmap)
mapImage <- get_map(location = c(lon = -96.8, lat = 32.7),
                    color = "color",
                    source = "google",
                    #maptype = "terrain",
                    zoom = 9)

cleaned_data <- broom::tidy(transformed_data)

library(RColorBrewer)
colors <- brewer.pal(9, "BuGn")

ggmap(mapImage) +
  geom_polygon(aes(x = long, y = lat, group = group),
               data = cleaned_data,
               color = colors[6], 
               fill = colors[5], 
               alpha = 1) +
  labs(x = 'longitude', y = 'latitude')

print('Test three, cleaning data')

formatted_data <- coarse_format_data(transformed_data,
                   test_inputs$sf,
                   test_inputs$mf,
                   test_inputs$or,
                   test_inputs$feature, 
                   separator = test_inputs$separator)

print(formatted_data[1:100, ])






  