# Compiled by Asher Spector - these are not officially a part of the shapefiles

# Some notes on methodology, for consistency. I classify:
# (a) Areas which allow both multifamily and single family 
# dwelling as multifamily.
# (b) Areas which allow mixed-use development (i.e. commercial + residential)
# as some type of residential area in most cases. 
# (c) If multifamily housing is permitted but with some restrictions, 
# the area is classified as multifamily housing. 

library(rgdal)
library(maptools)

# Inputs list -----------------------------------------------------------
# Name is the name of the municipality
# Path is the path to the shapefile
# sf is single family codes
# mf is multifamily codes
# or stands for other residential codes
# separator is the separator (str) by which to split codes into their component parts
# feature is the feature which actually lists the zoning codes per area
# -----------------------------------------------------------------------

# Basic projection inputs ----------------------------------------------- 

# For NAD83 - NO (HARN) - in ftUS
texas_north_central <- CRS('+init=EPSG:2276') #Dallas
texas_central <- CRS('+init=EPSG:2277') #Austin
texas_south_central <- CRS('+init=EPSG:2278') #Houston
lat_long <- CRS('+init=EPSG:4326') # standard lat and long


# AUSTIN -----------------------------------------------------------------

austin_inputs <- list(
  name = 'Austin',
  path = "austin_zoning/geo_export_8d172e2c-89a0-4c2f-a91d-337da740f107.shp",
  sf = 'SF', # The codes which correspond to single family
  mf = 'MF', # Codes which correspond to multifamily
  or = c('MH', 'RR', 'LA'), # Codes which correspond to other residential
  separator = '-',
  feature = 'zoning_zty',
  proj4string = lat_long
)

# DALLAS -----------------------------------------------------------------

north_texas_inputs <- list(
  name = 'North Texas', 
  path = '2015_North_Texas_Land_Use/2015_Land_Use.shp'
  
)


dallas_inputs <- list(
  name = 'Dallas',
  path = 'dallas_zoning/BaseZoning.shp',
  sf = c('A(A)', 'D(A)', 'TH', 'R'),
  mf = c('CH', 'MF'),
  or = c('MH(A)'),
  separator = '-',
  feature = 'LONG_ZONE_',
  proj4string = texas_north_central
)

denton_inputs <- list(
  name = 'Denton',
  path ='denton_zoning/Current_zoning.shp',
  feature = 'ZONING',
  sf = c('DR-1', 'NR-1', 'NR-2', 'NR-3', 'NR-4', 'NR-5', 'NR-6', 'RD-5'),
  mf = c('DR-2', 'MF-1', 'NRMU', 'NRMU-12', 'RCR-1', 'RCR-2', 'RCC-D', 'RCC-N'),
  or = c('PD', 'RD-5X'),
  separator = 'no_separator', 
  proj4string = texas_north_central
)


if (getOption('run.main', default=TRUE)) {
  
  print('hello')
  

}
