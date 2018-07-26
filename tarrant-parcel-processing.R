# This is a very quick script to process the Tarrant County Parcel data (bug in Fiona prevents me from reading it in Python)

library(tidyverse)
library(sf)

tarrant_county_parcel_path <- "data/parcels/tarrant_county_parcels_2018/TADData.gdb" 
outfile_path <- "data/parcels/processed_tarrant_county_parcels_2018/TADData.shp"

tarrant_data <- st_read(tarrant_county_parcel_path, 'TADParcels') %>%
  dplyr::select(c('Property_Class', 'Land_SqFt', 'TAXID'))

st_write(tarrant_data, outfile_path, delete_layer = TRUE) # Overwrites other data
