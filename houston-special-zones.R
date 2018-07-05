library(sf)
library(tidyverse)
#setwd("C:/Users/aspector/Documents/R Projects/TX Housing")


# For special minimum lot size
houston_spec_min_lot_path <- "data/Zoning Shapefiles/Houston_Special_Minimum_Lot_Size/Minimum_Lot_Size.shp"
houston_spec_min_lot <- st_read(houston_spec_min_lot_path)
sum(st_area(houston_spec_min_lot)) #Output is 6.87 square miles or 17790000 square meters

