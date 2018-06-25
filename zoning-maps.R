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
library(ggmap)
library(maptools)
library(tibble)
library(sf)


# Replace zone codes with sf, mf, or, etc.
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
    )) %>%
    lwgeom::st_make_valid()
    
  return(data)
}

# Full processing of data given an input list
fine_process_data <- function(input){
  # ---------------------------------------------------------
  # Input should be a list having
  # attributes as described in the zoning-maps-inputs script. 
  # ---------------------------------------------------------
  data <- st_read(input$path)
  data <- st_transform(data, '+init=EPSG:4326')
  processed_data <- coarse_format_data(data,input$sf, input$mf, input$or,
                                       input$feature, input$separator) %>%
    mutate(place = input$name)
  return(processed_data)
}


# Call fine_process_data repeatedly on a variety of inputs
process_all_data <- function(inputs){
  
  # ---------------------------------------------------------
  # Inputs should be a vector of lists, with each list having
  # attributes as described in the zoning-maps-inputs script. 
  # ---------------------------------------------------------
  
  if (length(inputs) < 2) {
    warning('Inputs must be greater than length 2; else use fine_process_data')
  }
  
  # Work with first input manually
  initial_input <- inputs[[1]]
  result <- fine_process_data(initial_input)
  
  # Rowbind the data from the rest of the inputs
  for (input in inputs[2:length(inputs)]){
    print(paste('Processing data from', input$name))
    result <- rbind(result, fine_process_data(input))
  }
  
  return(result)
  
}


# ----------------------------- AUSTIN --------------------------------------------------

inputs <- list(austin_inputs,
               round_rock_inputs,
               pflugerville_inputs, 
               georgetown_inputs,
               cedar_park_inputs,
               hutto_inputs) #Add austin inputs

all_data <- process_all_data(inputs)

colors <- c('red', 'lightskyblue2', 'pink', 'blue')

gg <- 
  ggmap::ggmap(mapImage) + 
  geom_sf(data = all_data, aes(fill = zone_code), alpha = 0.8, color = NA,
          inherit.aes = FALSE) +
  labs(x = 'Longitude', y = 'Latitude', title = 'Base Zones around Austin Texas',
       caption = 'Data from GIS Services of Austin, Round Rock, Georgetown,
       Pflugerville, Cedar Park, and Hetto') +
  scale_fill_manual(values = colors)


mapImage <- ggmap::get_map(location = c(lon = -97.7431, lat = 30.2672),
                           color = "bw",
                           source = "google",
                           zoom = 10)

#setwd("C:/Users/aspector/Documents/R Projects/TX Housing")

#png('Figures/austin-base-zones.png', w = 800, h = 800)
#print(gg)
#dev.off()

#svg('Figures/austin-base-zones.svg', w = 8, h = 8)
#print(gg)
#dev.off()



# ------------------------------- DALLAS -Now in Python ------------------------------

# Step 1: Load and format data

#dallas_raw_data <- st_read(north_texas_inputs$path, type = 0006) # ~0.6 GB, 45 sec

# Before graphing, simplify a bit. This does not eliminate zones or change their zoning
# category, but instead simplifies the edges of the zones to make graphing easier/possible.
# Moreover, this function doesn't work super well on lat/long coordinates, so before using it,
# we'll transform to another coordinate system (North Central Texas system), and then transform back
#dallas_raw_data <- st_transform(dallas_raw_data, '+init=EPSG:2276') # Takes about 20 seconds
#dallas_simplified_data <- st_simplify(dallas_raw_data, preserveTopology = TRUE, dTolerance = 3000) # Takes about 1-2 min, depending on dTolerance value#

#dallas_formatted_data <- coarse_format_data(dallas_raw_data, 
                                            #north_texas_inputs$sf,
                                            #north_texas_inputs$mf,
                                            #north_texas_inputs$or,
                                            #north_texas_inputs$feature,
                                            #north_texas_inputs$separator) 

#dallas_formatted_data <- dallas_formatted_data %>% filter(zone_code != 'NR')
#dallas_formatted_data$area <- as.numeric(st_area(dallas_formatted_data))
#dallas_formatted_data <- dallas_formatted_data %>% filter(area > 1000)

#dallas_formatted_data$area <- st_area(dallas_formatted_data)
#dallas_formatted_data <- dallas_formatted_data %>% arrange(desc(area))

#time0 <- proc.time()
#sample <- dallas_formatted_data %>%
#  sf::st_cast("MULTIPOLYGON") %>%
#  lwgeom::st_make_valid() %>%
#  dplyr::group_by(zone_code) %>%
#  dplyr::summarise(geometry = st_union(geometry)) %>%
#  dplyr::ungroup() 

#sample <- st_as_sf(sample)
#sample <- st_simplify(sample, preserveTopology = FALSE, dTolerance = 3000)
#sample <- st_transform(sample, '+init=EPSG:4326') # back to lat/long


# Step 2: Graph
#graph <- ggmap(mapImage) + 
  #geom_sf(data = sample, aes(fill = zone_code), color= NA, inherit.aes = FALSE) +
  #labs(x = 'longitude', y = 'latitude')
#png('All-Zones-Dallas.png', w = 800, h = 500)
#print(graph)
#dev.off()
#save(graph, 'All-Zones-Dallas.RData')
#print(proc.time() - time0)

#library(ggmap)
#mapImage <- get_map(location = c(lon = -96.8, lat = 32.8),
#                    color = "bw",
#                    source = "google",
#                    #maptype = "terrain",
#                    zoom = 9)