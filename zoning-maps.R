#setwd("C:/Users/aspector/Documents/R Projects/TX Housing/data/Zoning Shapefiles")
#setwd("C:/Users/amspe/Documents/R/MI2018/TXHousing/data/Zoning Shapefiles")

# This will prevent the zoning-maps-inputs script from doing extra work
# Make sure the working directory is set for this
#options(run.main = FALSE)
#source('zoning-maps-inputs.R')
#options(run.main = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(rgdal)
library(sf)
library(maptools)
library(tibble)
library(sf)


coarse_format_data <- function(data, sf, mf, or, feature, 
                               separator = '-'){
  
  # -----------------------------------------------------------
  # Formats shapefile data. Assumes that the raw_data is of the sf class.
  # The 'feature' should be the zoning data column name. 
  # -----------------------------------------------------------

  # Not converting to a data frame will preserve the geometry
  data <- data %>%
    select_(feature) %>%
    tidyr::separate_(feature, into = c('zone_code'), sep=separator) %>% 
    select('zone_code') %>% 
    mutate(zone_code = case_when(
      zone_code %in% sf ~ 'SF',
      zone_code %in% mf ~ 'MF',
      zone_code %in% or ~ 'OR',
      TRUE ~ 'NR'
    ))
    
  return(data)
}

#raw_dallas_data <- readShapePoly(dallas_inputs$path)
#raw_dallas_data@proj4string <- dallas_inputs$proj4string
#transformed_dallas_data <- spTransform(raw_dallas_data, lat_long)

# ------------------------------- DALLAS ------------------------------------------------

# Step 1: Load and format data

dallas_raw_data <- st_read(north_texas_inputs$path) # ~0.6 GB, 45 sec

# Before graphing, simplify a bit. This does not eliminate zones or change their zoning
# category, but instead simplifies the edges of the zones to make graphing easier/possible.
# Moreover, this function doesn't work super well on lat/long coordinates, so before using it,
# we'll transform to another coordinate system (North Central Texas system), and then transform back
dallas_raw_data <- st_transform(dallas_raw_data, '+init=EPSG:2276') # Takes about 20 seconds
dallas_simplified_data <- st_simplify(dallas_raw_data, preserveTopology = FALSE, dTolerance = 3000) # Takes about 1-2 min, depending on dTolerance value

dallas_formatted_data <- coarse_format_data(dallas_simplified_data, 
                                            north_texas_inputs$sf,
                                            north_texas_inputs$mf,
                                            north_texas_inputs$or,
                                            north_texas_inputs$feature,
                                            north_texas_inputs$separator)

dallas_formatted_data <- st_transform(dallas_formatted_data, '+init=EPSG:4326') # back to lat/long


# Step 2: Graph

library(ggmap)
mapImage <- get_map(location = c(lon = -96.8, lat = 32.7),
                    color = "color",
                    source = "google",
                    #maptype = "terrain",
                    zoom = 9)

# This call will take a while just because it's a ton of data/shapes (450K shapes)
ggmap(mapImage) + 
  geom_sf(data = dallas_formatted_data, aes(fill = zone_code), 
          alpha = 0.8, inherit.aes = FALSE) +
  labs(x = 'longitude', y = 'latitude')

# ----------------------------- AUSTIN --------------------------------------------------
